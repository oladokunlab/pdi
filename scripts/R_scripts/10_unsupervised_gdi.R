#!/usr/bin/env Rscript
# =============================================================================
# 10_unsupervised_gdi.R
#
# Build microbiome indices WITHOUT using the GHI, then ask whether they
# correlate with it.
#
# RUN:  Rscript 10_unsupervised_gdi.R
# TIME: 2-5 minutes.
#
# -----------------------------------------------------------------------------
# WHY THIS IS A DIFFERENT QUESTION
# -----------------------------------------------------------------------------
# Every GDI attempt so far selected taxa USING the GHI and then validated
# against the GHI. Even with leave-one-study-out and per-fold re-selection, the
# construction is bound to that target, and when the target changed the index
# vanished.
#
# The indices below use no outcome at all. Nothing is fitted, nothing is
# selected, no threshold is applied. They are computed from the community data
# alone, on definitions that predate this project. The correlation with the GHI
# is therefore genuine external validation: the GHI had no opportunity to shape
# what was built.
#
# That also means a null result here is more informative than a failed feature
# selection. If none of these track the GHI, the finding is that community
# structure does not carry gut health information in this corpus - not that a
# particular estimator failed to find it.
#
# -----------------------------------------------------------------------------
# THE FOUR MEASURES
# -----------------------------------------------------------------------------
# 1. DEVIATION FROM CONTROL. Within each study the control arms define what
#    normal looks like; every arm is scored by its Aitchison distance from that
#    centroid. This is the most literal reading of dysbiosis - departure from a
#    healthy reference - and is how the canine and feline indices are framed,
#    minus the training step.
#
#    Its weakness is direction: a beneficial additive also moves the community
#    away from control and would score as deviant. Measure 2 addresses that.
#
# 2. SIGNED POSITION ON THE CONTROL-TO-CHALLENGE AXIS. In studies containing
#    both controls and challenged arms, those two centroids define a direction.
#    Every arm is projected onto it, so moving towards the challenged end
#    scores positive and towards the control end negative. Direction is taken
#    from the study design, not from the GHI.
#
# 3. ALPHA DIVERSITY. Shannon and observed richness. "Reduced diversity
#    indicates dysbiosis" is the oldest heuristic in the field and is worth
#    testing rather than assuming.
#
# 4. FIRST PRINCIPAL COMPONENT of the within-study centred community matrix -
#    whatever the dominant axis of variation is, does it track gut health?
#
# -----------------------------------------------------------------------------
# HOW TO READ IT
# -----------------------------------------------------------------------------
# Expected signs, if a measure tracks gut health:
#   deviation from control   NEGATIVE with GHI (more deviation, worse health)
#   signed challenge axis    NEGATIVE with GHI (towards challenged, worse)
#   alpha diversity          POSITIVE with GHI (more diverse, better) - though
#                            this is contested and the test is the point
#   PC1                      either sign; only the magnitude is meaningful
#
# Correlations are reported both pooled and within study. The within-study
# figure is the one that matters: pooled correlation can be produced by
# between-study differences alone, which is what created the R2 = 0.000
# diagnosis in the supervised model.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

HEALTHY   <- "Control"
CHALLENGE <- c("Pathogen_Challenge", "Antibiotic_Challenge", "Dietary_Feed_Challenge")
ADDITIVE  <- "Non_Disruptive_Additive"

# =============================================================================
# 1 - LOAD THE COMMUNITY DATA, WITHOUT ANY OUTCOME
# =============================================================================
log_msg("[1/5] Loading community data (no GHI involved)")
clr_f <- file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv")
cnt_f <- file.path(R_OUT_DIR, "genus", "genus_counts.tsv")
if (!file.exists(clr_f)) stop("Not found: ", clr_f, "\n  Run 03_batch_correction.R first.")
clr <- fread(clr_f)
gcol <- names(clr)[1]
M <- as.matrix(clr[, -1, with = FALSE]); rownames(M) <- clr[[gcol]]
log_msg("  CLR matrix: %d genera x %d runs", nrow(M), ncol(M))

md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
runid <- trimws(as.character(md[[1]]))
r2s <- setNames(trimws(as.character(md$Study_ID)), runid)
r2t <- setNames(trimws(as.character(md$Treatment_Group)), runid)

runs <- colnames(M)
keep <- runs %in% names(r2s)
M <- M[, keep, drop = FALSE]; runs <- colnames(M)
log_msg("  %d runs mapped to a study and treatment group", length(runs))

# ---- collapse runs to treatment arms -----------------------------------------
arm_key <- paste(r2s[runs], r2t[runs], sep = "||")
uarm <- unique(arm_key)
A <- vapply(uarm, function(a) rowMeans(M[, arm_key == a, drop = FALSE]), numeric(nrow(M)))
A <- t(A)                                  # arms x genera
arm_study <- sub("\\|\\|.*$", "", rownames(A))
log_msg("  collapsed to %d arms across %d studies", nrow(A), length(unique(arm_study)))

# ---- intervention class per arm, from the biomarker sheet --------------------
bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
bk <- paste(trimws(bio$Study_ID), trimws(bio$Treatment_Group), sep = "||")
cls <- unname(tapply(trimws(as.character(bio$Challenge_Type)), bk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1]
})[rownames(A)])
log_msg("  class known for %d of %d arms (%d control, %d challenged, %d additive)",
        sum(!is.na(cls)), nrow(A), sum(cls %in% HEALTHY),
        sum(cls %in% CHALLENGE), sum(cls %in% ADDITIVE))

# =============================================================================
# 2 - THE FOUR MEASURES
# =============================================================================
log_msg("[2/5] Computing unsupervised indices")

n <- nrow(A)
dev_ctrl <- rep(NA_real_, n)   # 1. distance from the study's control centroid
sig_axis <- rep(NA_real_, n)   # 2. signed position on control->challenge axis

for (s in unique(arm_study)) {
  i  <- which(arm_study == s)
  ic <- i[cls[i] %in% HEALTHY]
  if (!length(ic)) next
  cen <- colMeans(A[ic, , drop = FALSE])
  dev_ctrl[i] <- sqrt(rowSums(sweep(A[i, , drop = FALSE], 2, cen, "-")^2))

  ix <- i[cls[i] %in% CHALLENGE]
  if (length(ix)) {
    axis <- colMeans(A[ix, , drop = FALSE]) - cen
    L <- sqrt(sum(axis^2))
    if (is.finite(L) && L > 0)
      sig_axis[i] <- as.numeric(sweep(A[i, , drop = FALSE], 2, cen, "-") %*% (axis / L))
  }
}
log_msg("  deviation from control: %d arms (studies with a control arm)", sum(!is.na(dev_ctrl)))
log_msg("  signed challenge axis : %d arms (studies with control AND challenge)",
        sum(!is.na(sig_axis)))

# ---- 3. alpha diversity, from counts -----------------------------------------
shannon <- rep(NA_real_, n); observed <- rep(NA_real_, n)
if (file.exists(cnt_f)) {
  cn <- fread(cnt_f)
  Cm <- as.matrix(cn[, -1, with = FALSE]); rownames(Cm) <- cn[[names(cn)[1]]]
  Cm <- Cm[, colnames(Cm) %in% runs, drop = FALSE]
  ck <- paste(r2s[colnames(Cm)], r2t[colnames(Cm)], sep = "||")
  for (j in seq_len(n)) {
    sel <- which(ck == rownames(A)[j])
    if (!length(sel)) next
    v <- rowSums(Cm[, sel, drop = FALSE])
    p <- v / sum(v); p <- p[p > 0]
    shannon[j]  <- -sum(p * log(p))
    observed[j] <- sum(v > 0)
  }
  log_msg("  alpha diversity: %d arms", sum(!is.na(shannon)))
} else log_msg("  counts table not found - alpha diversity skipped")

# ---- 4. PC1 of the within-study centred matrix -------------------------------
Ac <- A
for (s in unique(arm_study)) {
  i <- arm_study == s
  if (sum(i) >= 2) Ac[i, ] <- sweep(A[i, , drop = FALSE], 2, colMeans(A[i, , drop = FALSE]), "-")
  else Ac[i, ] <- NA
}
ok_pc <- stats::complete.cases(Ac)
pc1 <- rep(NA_real_, n)
if (sum(ok_pc) >= 10) {
  pr <- stats::prcomp(Ac[ok_pc, , drop = FALSE], center = TRUE, scale. = FALSE)
  pc1[ok_pc] <- pr$x[, 1]
  log_msg("  PC1: %d arms, explains %.1f%% of within-study variance",
          sum(ok_pc), 100 * (pr$sdev^2 / sum(pr$sdev^2))[1])
}

# =============================================================================
# 3 - ATTACH THE GHI, FOR CORRELATION ONLY
# =============================================================================
log_msg("[3/5] Attaching GHI scores (used for evaluation only, not construction)")
gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
gk <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
pull <- function(col) {
  if (!col %in% names(gh)) return(rep(NA_real_, n))
  v <- tapply(gh[[col]], gk, function(z) mean(z, na.rm = TRUE))
  unname(v[rownames(A)])
}
GHI <- list(Combined = pull("GHI_Combined_raw"),
            Morphology = pull("GHI_Morphology_raw"),
            Biomarker  = pull("GHI_Biomarker_raw"))
for (nm in names(GHI)) log_msg("  GHI_%-11s available for %d arms", nm, sum(!is.na(GHI[[nm]])))

# =============================================================================
# 4 - CORRELATIONS
# =============================================================================
log_msg("[4/5] Correlation with the GHI")

within_cor <- function(x, y, g) {
  xc <- x; yc <- y
  for (s in unique(g)) {
    i <- g == s & !is.na(x) & !is.na(y)
    if (sum(i) >= 2) { xc[i] <- x[i] - mean(x[i]); yc[i] <- y[i] - mean(y[i]) }
    else { xc[g == s] <- NA; yc[g == s] <- NA }
  }
  ok <- !is.na(xc) & !is.na(yc) & is.finite(xc) & is.finite(yc)
  if (sum(ok) < 10 || stats::sd(xc[ok]) == 0) return(c(NA_real_, NA_real_, sum(ok)))
  ct <- suppressWarnings(stats::cor.test(xc[ok], yc[ok], method = "spearman"))
  c(unname(ct$estimate), ct$p.value, sum(ok))
}

MEAS <- list("deviation from control" = dev_ctrl,
             "signed challenge axis"  = sig_axis,
             "Shannon diversity"      = shannon,
             "observed richness"      = observed,
             "community PC1"          = pc1)
EXPECT <- c("negative", "negative", "positive", "positive", "either")

res <- list()
log_msg("")
log_msg("  %-24s %-11s %9s %9s %6s  %s", "measure", "GHI", "pooled", "within", "n", "expected")
for (mi in seq_along(MEAS)) {
  mn <- names(MEAS)[mi]; x <- MEAS[[mi]]
  if (all(is.na(x))) next
  for (gn in names(GHI)) {
    y <- GHI[[gn]]
    ok <- !is.na(x) & !is.na(y)
    pooled <- if (sum(ok) >= 10) suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman")) else NA_real_
    w <- within_cor(x, y, arm_study)
    res[[paste(mn, gn)]] <- list(meas = mn, ghi = gn, pooled = pooled,
                                 within = w[1], p = w[2], n = w[3])
    log_msg("  %-24s %-11s %+9.3f %+9.3f %6d  %s%s", mn, gn,
            ifelse(is.na(pooled), NA, pooled), ifelse(is.na(w[1]), NA, w[1]), w[3],
            EXPECT[mi],
            if (!is.na(w[2]) && w[2] < 0.05) "  *" else "")
  }
}
log_msg("")
log_msg("  * marks within-study p < 0.05. The within-study column is the one that")
log_msg("  matters; a pooled correlation can arise from between-study differences")
log_msg("  alone.")

# =============================================================================
# 5 - READING
# =============================================================================
log_msg("[5/5] Reading")
sig <- Filter(function(z) !is.na(z$p) && z$p < 0.05 && abs(z$within) >= 0.2, res)
# direction check: deviation and signed axis should be negative, diversity positive
right_dir <- Filter(function(z) {
  e <- EXPECT[match(z$meas, names(MEAS))]
  (e == "negative" && z$within < 0) || (e == "positive" && z$within > 0) || e == "either"
}, sig)

verdict <- if (!length(sig)) {
  paste("No unsupervised community measure correlates with the GHI within study.",
        "Community structure does not carry gut health information in this corpus,",
        "which is a stronger statement than a feature-selection failure because",
        "nothing was fitted or selected here.")
} else if (!length(right_dir)) {
  paste("Correlations reach significance but run opposite to the expected",
        "direction, which argues for an artefact rather than a signal. Inspect",
        "before interpreting.")
} else {
  paste(length(right_dir), "measure/target combination(s) correlate with the GHI",
        "within study in the expected direction. This is unsupervised, so it is",
        "external validation rather than a fitted result, and is worth reporting",
        "and following up.")
}
log_msg("  %s", verdict)

out <- c("UNSUPERVISED MICROBIOME INDICES AGAINST THE GHI", strrep("=", 74), "",
  "Indices computed from community data alone. No outcome was used in their",
  "construction, nothing was fitted, and no threshold was applied, so the",
  "correlations below are external validation rather than model fit.", "",
  sprintf("Arms: %d across %d studies", nrow(A), length(unique(arm_study))), "",
  sprintf("  %-24s %-11s %9s %9s %6s", "measure", "GHI", "pooled", "within", "n"))
for (z in res)
  out <- c(out, sprintf("  %-24s %-11s %+9.3f %+9.3f %6d", z$meas, z$ghi,
                        ifelse(is.na(z$pooled), NA, z$pooled),
                        ifelse(is.na(z$within), NA, z$within), z$n))
out <- c(out, "", verdict, "",
  "Expected directions: deviation from control and signed challenge axis should",
  "correlate NEGATIVELY with the GHI; alpha diversity POSITIVELY. The community",
  "PC1 has no expected sign and only its magnitude is interpretable.")
writeLines(out, file.path(R_OUT_DIR, "gdi", "unsupervised_gdi.txt"))

fwrite(data.table(arm = rownames(A), study = arm_study, class = cls,
                  deviation = dev_ctrl, signed_axis = sig_axis,
                  shannon = shannon, observed = observed, pc1 = pc1,
                  GHI_Combined = GHI$Combined),
       file.path(R_OUT_DIR, "gdi", "unsupervised_scores.tsv"), sep = "\t")
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "unsupervised_gdi.txt"))
log_msg("scores  -> %s", file.path(R_OUT_DIR, "gdi", "unsupervised_scores.tsv"))
log_msg("=== DONE ===")
