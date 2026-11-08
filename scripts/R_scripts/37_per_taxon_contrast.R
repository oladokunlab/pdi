#!/usr/bin/env Rscript
# =============================================================================
# 37_per_taxon_contrast.R
#
# Per-genus within-study CONTRAST, pooled across all eligible studies -
# not the within-study CORRELATION used in 18_taxon_metaanalysis.R.
#
# RUN:  Rscript 37_per_taxon_contrast.R
#
# -----------------------------------------------------------------------------
# WHY THIS EXISTS
# -----------------------------------------------------------------------------
# 18_taxon_metaanalysis.R computes, per genus, a within-study Spearman
# correlation between that genus's abundance and the GHI. A correlation needs
# at least three arms in one study to mean anything - two points always sit
# on a line, so a "correlation" from two arms is mathematically forced to
# +1, -1 or undefined, not a real measurement. Requiring three arms cut 41
# eligible studies down to 8, because most studies in this corpus are simple
# two-arm designs (one control, one treatment).
#
# That three-arm requirement was never asked for as a design choice - it is
# a mathematical necessity for correlation specifically, not for every kind
# of within-study test. A CONTRAST does not have that requirement. It only
# needs a study to have a Control arm and one arm of the type being
# contrasted - exactly the design already used for GHI validation
# (04e_mixed_model_validation.R) and the guild ratios (33_guild_ratios.R).
# This script applies that same contrast-and-pool design to each of the 84
# genera individually, so far more of the 44-study corpus can contribute.
#
# -----------------------------------------------------------------------------
# WHAT CHANGES AND WHAT DOES NOT
# -----------------------------------------------------------------------------
# Same genus abundances (CLR, batch-corrected). Same random-effects pooling
# (DerSimonian-Laird, weighted by number of group means) as every other
# meta-analysis in this project. What changes is the per-study statistic:
# a paired arm-level DIFFERENCE instead of a CORRELATION, and consequently a
# much larger and different set of eligible studies per genus.
#
# This is a different question from 18's, not a corrected version of it.
# 18 asks: within a study, does this genus track the GHI continuously.
# This script asks: does this genus's level shift between control and
# treated arms, on average, across as much of the corpus as will support it.
# Both are legitimate; they can validly disagree for a given genus, and
# where they do, that disagreement is informative rather than a bug in
# either one.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

CONTROL   <- "Control"
ADDITIVE  <- "Non_Disruptive_Additive"
CHALLENGE <- c("Pathogen_Challenge", "Antibiotic_Challenge", "Dietary_Feed_Challenge")

# same 22 genera evaluated in 18_taxon_metaanalysis.R, for direct comparison.
# Set to NULL to run all 84 genera instead.
CANDIDATE_GENERA <- NULL

log_msg("[1/5] Loading genus abundances and metadata")
clr <- as.data.frame(fread(file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv")))
gname_col <- names(clr)[1]
genera <- clr[[gname_col]]
X <- as.matrix(clr[, -1, drop = FALSE])   # genera x runs
rownames(X) <- genera
runs <- colnames(X)

md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2study <- setNames(trimws(as.character(md$Study_ID)), rid)
r2arm   <- setNames(trimws(as.character(md$Treatment_Group)), rid)

# Challenge_Type does not exist on the raw QIIME2 metadata file - it is never
# meant to. GDI_metadata_qiime2.tsv only carries Intervention_Class and
# Challenge_Subtype, the granular fields. The simplified five-category
# Challenge_Type (Control, Non_Disruptive_Additive, Pathogen_Challenge,
# Antibiotic_Challenge, Dietary_Feed_Challenge) is built once, upstream,
# during GHI construction, and lives on ghi_scores_pca.tsv. Every other
# script that uses Challenge_Type (08_loso_validate.R, 33_guild_ratios.R,
# 29_ghi_diagnostics.R) reads it from there, keyed by Study_ID +
# Treatment_Group, not from the raw QIIME2 file. This script originally
# tried to read it from METADATA_TSV directly, which was the actual bug -
# nothing was renamed or lost, the join was just pointed at the wrong file.
ghi <- as.data.frame(fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")))
arm2cls <- setNames(trimws(as.character(ghi$Challenge_Type)),
                     paste(trimws(as.character(ghi$Study_ID)),
                           trimws(as.character(ghi$Treatment_Group)), sep = "||"))
r2cls <- unname(arm2cls[paste(r2study, r2arm, sep = "||")])
names(r2cls) <- rid
log_msg("  Challenge_Type joined from ghi_scores_pca.tsv: %d of %d runs matched",
        sum(!is.na(r2cls)), length(r2cls))

if (is.null(CANDIDATE_GENERA)) {
  CANDIDATE_GENERA <- rownames(X)
  log_msg("  testing all %d genera", length(CANDIDATE_GENERA))
} else {
  log_msg("  testing %d candidate genera (same set as 18_taxon_metaanalysis.R)", length(CANDIDATE_GENERA))
}

log_msg("[2/5] Collapsing runs to arms (mean CLR per Study_ID x Treatment_Group)")
common <- intersect(runs, names(r2study))
arm_key <- paste(r2study[common], r2arm[common], sep = "||")
uniq_arms <- unique(arm_key)
A <- matrix(NA_real_, nrow = nrow(X), ncol = length(uniq_arms),
            dimnames = list(rownames(X), uniq_arms))
for (a in uniq_arms) {
  cols <- common[arm_key == a]
  A[, a] <- rowMeans(X[, cols, drop = FALSE], na.rm = TRUE)
}
arm_study <- sub("\\|\\|.*$", "", uniq_arms)
arm_cls   <- r2cls[sapply(strsplit(uniq_arms, "\\|\\|"), function(z)
                    common[which(arm_key == paste(z, collapse = "||"))[1]])]
names(arm_cls) <- uniq_arms
log_msg("  %d arms across %d studies", length(uniq_arms), length(unique(arm_study)))

# =============================================================================
# per-study, per-genus contrast: mean(treated arms) - mean(control arms)
# =============================================================================
compute_contrasts <- function(genus, treated_classes) {
  vals <- A[genus, ]
  out <- data.frame(study = character(0), diff = numeric(0), n_ctrl = integer(0), n_trt = integer(0))
  for (s in unique(arm_study)) {
    idx <- arm_study == s
    ctrl_idx <- idx & arm_cls[uniq_arms] == CONTROL
    trt_idx  <- idx & arm_cls[uniq_arms] %in% treated_classes
    ctrl_idx[is.na(ctrl_idx)] <- FALSE; trt_idx[is.na(trt_idx)] <- FALSE
    nc <- sum(ctrl_idx & !is.na(vals)); nt <- sum(trt_idx & !is.na(vals))
    if (nc >= 1 && nt >= 1) {
      d <- mean(vals[trt_idx], na.rm = TRUE) - mean(vals[ctrl_idx], na.rm = TRUE)
      out <- rbind(out, data.frame(study = s, diff = d, n_ctrl = nc, n_trt = nt))
    }
  }
  out
}

# random-effects pooling, DerSimonian-Laird, weight = n_ctrl + n_trt
# (number of group means behind the estimate) - matches 04e and 33's method.
pool_re <- function(df) {
  k <- nrow(df)
  if (k < 2) return(list(k = k, est = NA, p = NA, i2 = NA))
  y  <- df$diff
  w0 <- df$n_ctrl + df$n_trt
  mu_f <- sum(w0 * y) / sum(w0)
  Q <- sum(w0 * (y - mu_f)^2)
  C <- sum(w0) - sum(w0^2) / sum(w0)
  tau2 <- max(0, (Q - (k - 1)) / C)
  wi <- 1 / (1 / w0 + tau2)
  mu <- sum(wi * y) / sum(wi)
  se <- sqrt(1 / sum(wi))
  z <- mu / se
  p <- 2 * pnorm(-abs(z))
  i2 <- max(0, (Q - (k - 1)) / Q) * 100
  list(k = k, est = mu, p = p, i2 = i2)
}

log_msg("[3/5] Additive vs control, all genera")
res_add <- lapply(CANDIDATE_GENERA, function(g) pool_re(compute_contrasts(g, ADDITIVE)))
names(res_add) <- CANDIDATE_GENERA

log_msg("[4/5] Challenge vs control, all genera")
res_chal <- lapply(CANDIDATE_GENERA, function(g) pool_re(compute_contrasts(g, CHALLENGE)))
names(res_chal) <- CANDIDATE_GENERA

log_msg("[5/5] Combined (additive OR challenge) vs control, all genera")
res_comb <- lapply(CANDIDATE_GENERA, function(g) pool_re(compute_contrasts(g, c(ADDITIVE, CHALLENGE))))
names(res_comb) <- CANDIDATE_GENERA

# =============================================================================
build_table <- function(res_list, label) {
  data.frame(
    genus = names(res_list),
    contrast = label,
    k   = vapply(res_list, function(x) x$k, numeric(1)),
    estimate = vapply(res_list, function(x) x$est, numeric(1)),
    p   = vapply(res_list, function(x) x$p, numeric(1)),
    I2  = vapply(res_list, function(x) x$i2, numeric(1))
  )
}
tab <- rbind(build_table(res_add, "ADDITIVE"),
             build_table(res_chal, "CHALLENGE"),
             build_table(res_comb, "COMBINED"))
tab <- tab[order(tab$contrast, tab$p), ]

out_dir <- file.path(R_OUT_DIR, "gdi")
fwrite(tab, file.path(out_dir, "taxon_contrast_meta.tsv"), sep = "\t")

log_msg("")
log_msg("=========================================================")
log_msg("  k (studies contributing) BY CONTRAST - compare to 18's fixed k=8")
log_msg("=========================================================")
for (lab in c("ADDITIVE", "CHALLENGE", "COMBINED")) {
  ks <- tab$k[tab$contrast == lab]
  log_msg("  %-10s median k = %d, range %d-%d across %d genera",
          lab, median(ks, na.rm=TRUE), min(ks, na.rm=TRUE), max(ks, na.rm=TRUE), length(ks))
}

log_msg("")
log_msg("=========================================================")
log_msg("  TOP RESULTS BY CONTRAST (lowest p, k >= 5)")
log_msg("=========================================================")
for (lab in c("ADDITIVE", "CHALLENGE", "COMBINED")) {
  sub <- tab[tab$contrast == lab & tab$k >= 5, ]
  sub <- sub[order(sub$p), ]
  log_msg("  --- %s ---", lab)
  log_msg("  %-24s %5s %10s %10s %8s", "genus", "k", "estimate", "p", "I2")
  for (i in seq_len(min(8, nrow(sub)))) {
    r <- sub[i, ]
    log_msg("  %-24s %5d %10.3f %10.4f %7.0f%%%s", r$genus, r$k, r$estimate, r$p, r$I2,
            if (!is.na(r$p) && r$p <= 0.05) "  *" else "")
  }
  log_msg("")
}

log_msg("=== DONE ===")
log_msg("Full table -> %s", file.path(out_dir, "taxon_contrast_meta.tsv"))
log_msg("")
log_msg("READING: k here can exceed 20-30 per contrast where 18_taxon_metaanalysis.R")
log_msg("was fixed at k=8 for every genus. The two scripts answer different")
log_msg("questions - continuous within-study tracking (18) versus a group-level")
log_msg("shift (this script) - and are not expected to rank genera identically.")
