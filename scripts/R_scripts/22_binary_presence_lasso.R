#!/usr/bin/env Rscript
# =============================================================================
# 22_binary_presence_lasso.R
#
# Does presence/absence of a genus predict GHI better than its continuous
# abundance? Same target, same rules, only the feature representation changes.
#
# RUN:  Rscript 22_binary_presence_lasso.R   (needs 05_maaslin2.R to have run,
#       for its run-to-arm-to-GHI mapping; reads raw counts independently)
# TIME: a few minutes.
#
# -----------------------------------------------------------------------------
# WHY BINARIZE FROM RAW COUNTS, NOT THE ALREADY-BATCH-CORRECTED MATRIX
# -----------------------------------------------------------------------------
# The idea (from GMWI2, Chang et al. 2024) is that presence/absence sidesteps
# sequencing depth and platform differences BY ITSELF, without needing
# continuous-value batch correction at all. Binarizing the CLR matrix - which
# has already been through MMUPHin and zero-replacement - would test something
# else entirely (a discretized version of an already-corrected value) and
# defeat the point of the comparison. This script goes back to raw counts.
#
# Presence is defined as relative abundance >= 0.001% within that run,
# matching GMWI2's own threshold - not simply "count > 0", which would call a
# single stray or index-hopped read "present."
#
# -----------------------------------------------------------------------------
# ARM-LEVEL AGGREGATION: DETECTION FREQUENCY, NOT A SECOND THRESHOLD
# -----------------------------------------------------------------------------
# An arm contains multiple sequencing runs. Rather than force a second,
# arbitrary threshold ("present in the arm" = present in how many runs?), each
# arm gets the FRACTION of its runs where the genus was detected - a
# continuous value from 0 to 1. This loses no information the binary
# definition doesn't already lose, and avoids compounding one threshold
# decision with another.
#
# -----------------------------------------------------------------------------
# EVERYTHING ELSE IS IDENTICAL TO 06_lasso_selection.R
# -----------------------------------------------------------------------------
# Same target (GHI, within-study centered), same 84-genus universe, same
# arms, same 1000-iteration/50%-subsample stability selection, same 70%
# registered threshold. Only X changes: detection frequency instead of CLR
# abundance. A genuinely controlled, one-variable comparison.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

PRESENCE_THRESHOLD <- 0.00001   # 0.001% relative abundance, per GMWI2
N_BOOT         <- 1000
STABILITY_THR  <- 0.70
SUBSAMPLE_FRAC <- 0.50
set.seed(42)
UNASSIGNABLE <- c("g__uncultured","g__unclassified","g__uncultured_bacterium",
                  "g__metagenome","g__Incertae_Sedis","g__unidentified")

log_msg("[1/6] Loading the run-to-arm-to-GHI mapping (same as 06_lasso_selection.R)")
inp <- file.path(R_OUT_DIR, "gdi", "gdi_input.rds")
if (!file.exists(inp)) stop("Not found: ", inp, "\n  Run 05_maaslin2.R first.")
d <- readRDS(inp)
clr_ref <- d$clr; meta <- d$meta
run_universe <- colnames(clr_ref)
genus_universe <- rownames(clr_ref)
log_msg("  reference: %d genera x %d runs (same universe as the continuous analysis)",
        length(genus_universe), length(run_universe))

log_msg("[2/6] Loading raw counts and computing per-run relative abundance")
raw <- fread(file.path(R_OUT_DIR, "genus", "genus_counts.tsv"))
gcol <- names(raw)[1]
raw_mat <- as.matrix(raw[, -1, with = FALSE]); rownames(raw_mat) <- raw[[gcol]]
common_runs <- intersect(run_universe, colnames(raw_mat))
common_genera <- intersect(genus_universe, rownames(raw_mat))
raw_mat <- raw_mat[common_genera, common_runs, drop = FALSE]
log_msg("  %d genera x %d runs matched to the reference universe", nrow(raw_mat), ncol(raw_mat))

lib_size <- colSums(raw_mat)
rel_abund <- sweep(raw_mat, 2, lib_size, "/")
presence <- (rel_abund >= PRESENCE_THRESHOLD) * 1
log_msg("  presence threshold: relative abundance >= %.4f%%", PRESENCE_THRESHOLD * 100)
log_msg("  overall detection rate across the matrix: %.1f%%", 100 * mean(presence))

log_msg("[3/6] Collapsing runs to study x arm (detection frequency, not CLR mean)")
if (!"Treatment_Group" %in% names(meta)) {
  cw <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = "Crosswalk_16S"),
                      stringsAsFactors = FALSE)
  meta$Treatment_Group <- cw$Treatment_Group[match(rownames(meta), trimws(cw$Run))]
}
meta <- meta[common_runs, , drop = FALSE]
arm <- paste(meta$Study_ID, meta$Treatment_Group, sep = "||")
ok <- !is.na(arm) & !is.na(meta$GHI)
presence <- presence[, ok, drop = FALSE]; meta <- meta[ok, , drop = FALSE]; arm <- arm[ok]

ua <- unique(arm)
X <- t(sapply(ua, function(a) rowMeans(presence[, arm == a, drop = FALSE], na.rm = TRUE)))
y <- vapply(ua, function(a) mean(meta$GHI[arm == a]), numeric(1))
study <- vapply(ua, function(a) meta$Study_ID[arm == a][1], character(1))

drop <- !stats::complete.cases(X) | is.na(y)
if (any(drop)) { X <- X[!drop, , drop = FALSE]; y <- y[!drop]; study <- study[!drop] }
log_msg("  %d arms across %d studies (matches the continuous analysis if binarization",
        nrow(X), length(unique(study)))
log_msg("  changed nothing about which arms qualify)")

# drop zero-variance columns (genus detected in 0% or 100% of every arm - no
# information for LASSO to use)
keep_x <- apply(X, 2, function(v) stats::sd(v) > 0)
X <- X[, keep_x, drop = FALSE]
log_msg("  %d genera retained after removing zero-variance detection patterns", ncol(X))

log_msg("[4/6] Within-study centering (same reasoning as 06_lasso_selection.R)")
Xc <- X; yc <- y
for (s in unique(study)) { i <- study == s
  Xc[i, ] <- sweep(X[i, , drop = FALSE], 2, colMeans(X[i, , drop = FALSE]), "-")
  yc[i] <- y[i] - mean(y[i]) }

log_msg("[5/6] Stability selection, %d resamples, %.0f%% study subsampling",
        N_BOOT, 100 * SUBSAMPLE_FRAC)
studies <- unique(study); ntake <- max(3, floor(length(studies) * SUBSAMPLE_FRAC))
cnt <- setNames(numeric(ncol(Xc)), colnames(Xc)); nok <- 0
for (b in seq_len(N_BOOT)) {
  tk <- sample(studies, ntake); idx <- which(study %in% tk)
  if (length(idx) < 10 || length(unique(yc[idx])) < 5) next
  f <- tryCatch(cv.glmnet(Xc[idx, , drop = FALSE], yc[idx], alpha = 1,
                          nfolds = min(5, length(idx)), standardize = TRUE, intercept = FALSE),
               error = function(e) NULL)
  if (is.null(f)) next
  co <- as.matrix(coef(f, s = f$lambda.1se)); co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  sel <- rownames(co)[co[,1] != 0]
  cnt[sel] <- cnt[sel] + 1; nok <- nok + 1
}
freq <- if (nok > 0) cnt / nok else cnt
freq <- freq[!(names(freq) %in% UNASSIGNABLE)]
log_msg("  %d of %d resamples usable", nok, N_BOOT)
log_msg("  genera at the registered %.0f%% threshold: %d", 100 * STABILITY_THR, sum(freq >= STABILITY_THR))
top <- sort(freq, decreasing = TRUE)[1:min(10, length(freq))]
log_msg("  strongest:")
for (i in seq_along(top)) log_msg("    %-32s %3.0f%%", gsub("^g__", "", names(top)[i]), 100 * top[i])

log_msg("[6/6] In-sample fit and writing")
full <- tryCatch({ f <- cv.glmnet(Xc, yc, alpha = 1, standardize = TRUE, intercept = FALSE)
                   pr <- as.numeric(predict(f, newx = Xc, s = f$lambda.1se))
                   r2 <- if (sd(pr) == 0) 0 else max(0, 1 - sum((yc - pr)^2) / sum((yc - mean(yc))^2))
                   list(r2 = r2) }, error = function(e) list(r2 = NA_real_))
log_msg("  R2 = %.3f  (continuous-CLR comparison: 0.081 on the same target)", full$r2)

out_dir <- file.path(R_OUT_DIR, "gdi")
out <- c("BINARY PRESENCE/ABSENCE LASSO (vs continuous CLR abundance)",
  strrep("=", 74), "",
  sprintf("Presence threshold: relative abundance >= %.4f%% (GMWI2 convention)", PRESENCE_THRESHOLD*100),
  sprintf("Arms: %d, studies: %d, genera: %d", nrow(X), length(unique(study)), ncol(X)),
  sprintf("In-sample R2: %.3f  (continuous-CLR comparison on the same target: 0.081)", full$r2),
  sprintf("Genera at %.0f%% stability threshold: %d", 100*STABILITY_THR, sum(freq >= STABILITY_THR)),
  "", "Strongest candidates:")
for (i in seq_along(top)) out <- c(out, sprintf("  %-32s %3.0f%%", gsub("^g__","",names(top)[i]), 100*top[i]))
writeLines(out, file.path(out_dir, "binary_presence_lasso.txt"))
log_msg("summary -> %s", file.path(out_dir, "binary_presence_lasso.txt"))
log_msg("=== DONE ===")
