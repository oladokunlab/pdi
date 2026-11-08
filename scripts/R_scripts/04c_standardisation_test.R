#!/usr/bin/env Rscript
# =============================================================================
# 04c_standardisation_test.R
#
# Diagnostic, not a pipeline stage. Answers one question:
#
#   Is the weak control-vs-challenge discrimination caused by the
#   standardisation method, or is it inherent to the data?
#
# RUN:  Rscript 04c_standardisation_test.R
# TIME: under a minute. Reads the biomarker workbook only - no sequence data,
#       no cluster job, nothing else in the pipeline is touched or changed.
#
# -----------------------------------------------------------------------------
# WHY THIS TEST EXISTS
# -----------------------------------------------------------------------------
# Each row in the biomarker sheet is a TREATMENT-GROUP MEAN extracted from a
# published table or figure, not an individual bird - papers almost never
# publish bird-level data. With roughly 11 rows per study spread across 3-4
# treatment groups, z-scoring within a study computes a standard deviation
# from only 3-4 numbers. That is unstable enough to plausibly erode a real
# biological signal.
#
# Supporting evidence that the signal exists in the raw data: on untransformed
# villus height, control groups exceed challenged groups in 14 of 19 studies
# that contain both arms. If the z-scored version performs worse than that,
# the standardisation is destroying information rather than revealing it.
#
# -----------------------------------------------------------------------------
# THE THREE METHODS COMPARED
# -----------------------------------------------------------------------------
# RAW         Values as reported, oriented only (multiplied by +1/-1 so higher
#             always means healthier). No standardisation at all. Cannot be
#             pooled across studies - units differ - but IS valid for the
#             within-study comparison this script performs.
#
# ZSCORE      (value - study mean) / study SD, computed within each study.
#             The method currently used in the pipeline.
#
# PERCENTILE  Each value expressed as its percentile rank within its own
#             study, scaled 0-1. Rank-based, so it is unaffected by outliers
#             and by a standard deviation estimated from 3-4 points. This is
#             the method used by the Gut Microbiome Index (Oh et al. 2022),
#             which reached cross-study AUC 0.71 on 1,663 samples, and was
#             the approach the literature review specifically recommended for
#             block-wise missing multi-study data.
#
# -----------------------------------------------------------------------------
# HOW PERFORMANCE IS MEASURED
# -----------------------------------------------------------------------------
# Within-study concordance: inside a single study, does the untreated control
# group score healthier than that same study's challenged groups? Reported as
# the proportion of eligible studies where that holds.
#
# This is a paired, within-study comparison, so it is valid for all three
# methods including RAW. Pooled AUC is deliberately NOT the headline metric:
# within-study standardisation centres every study on zero, so a control in a
# beneficial-additive study scores below its study mean while a control in a
# challenge study scores above it. Pooling cancels the signal and reports ~0.5
# regardless of whether the index works.
#
# The null is 50% (a coin flip), the same null as AUC 0.5, so the numbers are
# interpretable on the same scale.
#
# Beneficial additives (probiotics, phytogenics) are excluded from the
# contrast entirely rather than counted as challenges - a bird on a probiotic
# is not diseased, and treating it as such previously inverted the result.
# =============================================================================

.locate_and_source_config <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  d <- if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
  source(file.path(d, "00_r_config.R"))
}
.locate_and_source_config()

HEALTHY_LEVELS   <- c("Control")
CHALLENGE_LEVELS <- c("Pathogen_Challenge", "Antibiotic_Challenge", "Dietary_Feed_Challenge")

DOMAINS <- list(
  Morphology     = c(VH_um = +1, CD_um = -1),
  Barrier        = c(ZO1 = +1, Ocln = +1, Cldn1 = +1),
  Inflammation   = c(IL6 = -1, TNFa = -1),
  SCFA           = c(Ace = +1, Prop = +1, Buty = +1),
  # Split mucosal / systemic to match 04b_ghi_pca.R. Combining them scores a
  # domain the index does not actually use (see 04g: combined ran r = -0.247
  # against Morphology). Systemic +1 orients for PCA only - it is not a claim
  # that higher circulating immunoglobulin is healthier.
  Immunoglobulin_Mucosal  = c(sIgA = +1, IgA = +1),
  Immunoglobulin_Systemic = c(IgG  = +1, IgM = +1)
)

log_msg("[1/4] Loading biomarker sheet")
bio <- as.data.frame(read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
log_msg("  %d rows x %d columns", nrow(bio), ncol(bio))

if (!"Challenge_Type" %in% names(bio)) stop("Challenge_Type column required for this test.")

# How many treatment groups does each study actually contribute? This is the
# quantity that determines whether a within-study SD is estimable at all.
grp_per_study <- tapply(seq_len(nrow(bio)), bio$Study_ID, length)
log_msg("  rows per study: median %.0f, range %d-%d",
        stats::median(grp_per_study), min(grp_per_study), max(grp_per_study))
log_msg("  (a within-study SD computed from this few points is the thing under test)")

# ---- the three standardisation methods ------------------------------------
std_raw <- function(v, study) as.numeric(v)

std_zscore <- function(v, study) {
  v <- as.numeric(v); out <- rep(NA_real_, length(v))
  for (s in unique(study)) {
    k <- study == s & !is.na(v)
    if (sum(k) < 2) next
    sdv <- stats::sd(v[k])
    if (!is.finite(sdv) || sdv == 0) next
    out[k] <- (v[k] - mean(v[k])) / sdv
  }
  out
}

std_percentile <- function(v, study) {
  v <- as.numeric(v); out <- rep(NA_real_, length(v))
  for (s in unique(study)) {
    k <- study == s & !is.na(v)
    n <- sum(k)
    if (n < 2) next
    # rank-based: robust to outliers and to an SD estimated from few points
    out[k] <- (rank(v[k], ties.method = "average") - 0.5) / n
  }
  out
}

METHODS <- list(RAW = std_raw, ZSCORE = std_zscore, PERCENTILE = std_percentile)

# ---- within-study concordance ---------------------------------------------
concordance <- function(score, study, is_ctrl, is_chal) {
  ok <- !is.na(score)
  n_eval <- 0; n_ok <- 0
  for (s in unique(study[ok])) {
    k <- ok & study == s
    cv <- score[k & is_ctrl]; hv <- score[k & is_chal]
    if (!length(cv) || !length(hv)) next
    n_eval <- n_eval + 1
    if (mean(cv) > mean(hv)) n_ok <- n_ok + 1
  }
  list(n_eval = n_eval, n_ok = n_ok,
       prop = if (n_eval) n_ok / n_eval else NA_real_)
}

# Binomial test against a 50% null, so a result can be called better than
# chance rather than merely looking better.
binom_p <- function(k, n) if (n > 0) stats::binom.test(k, n, 0.5)$p.value else NA_real_

is_ctrl <- bio$Challenge_Type %in% HEALTHY_LEVELS
is_chal <- bio$Challenge_Type %in% CHALLENGE_LEVELS
log_msg("  contrast: %d control rows vs %d challenged rows (%d additive/other excluded)",
        sum(is_ctrl), sum(is_chal), sum(!is_ctrl & !is_chal))

# =============================================================================
log_msg("[2/4] Per-marker concordance under each standardisation")
# Tested per marker first, before any domain aggregation, so that a weak
# result can be attributed to the standardisation rather than to how domains
# were combined.
# =============================================================================
marker_rows <- list()
for (dom in names(DOMAINS)) {
  for (m in names(DOMAINS[[dom]])) {
    if (!m %in% names(bio)) next
    dir <- DOMAINS[[dom]][[m]]
    res <- c()
    for (mn in names(METHODS)) {
      sc <- METHODS[[mn]](bio[[m]], bio$Study_ID) * dir
      cc <- concordance(sc, bio$Study_ID, is_ctrl, is_chal)
      res[mn] <- if (cc$n_eval >= 3) cc$prop else NA_real_
      if (mn == "RAW") res["n"] <- cc$n_eval
    }
    marker_rows[[m]] <- c(domain = dom, marker = m, res)
  }
}

log_msg("  %-8s %-16s %5s  %7s %7s %7s", "domain", "marker", "n", "RAW", "ZSCORE", "PCTILE")
for (r in marker_rows) {
  fmt <- function(v) if (is.na(v)) "    -" else sprintf("%6.0f%%", 100 * as.numeric(v))
  log_msg("  %-8s %-16s %5s  %7s %7s %7s",
          substr(r["domain"], 1, 8), r["marker"], r["n"],
          fmt(r["RAW"]), fmt(r["ZSCORE"]), fmt(r["PERCENTILE"]))
}

# =============================================================================
log_msg("[3/4] Domain-level and index-level concordance under each method")
# =============================================================================
build_index <- function(std_fn) {
  dom_scores <- list()
  for (dom in names(DOMAINS)) {
    cols <- intersect(names(DOMAINS[[dom]]), names(bio))
    if (!length(cols)) next
    M <- sapply(cols, function(m) std_fn(bio[[m]], bio$Study_ID) * DOMAINS[[dom]][[m]])
    if (is.null(dim(M))) M <- matrix(M, ncol = 1)
    n_present <- rowSums(!is.na(M))
    s <- ifelse(n_present >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
    s[!is.finite(s)] <- NA_real_
    dom_scores[[dom]] <- s
  }
  D <- do.call(cbind, dom_scores)
  morph <- dom_scores[["Morphology"]]
  bio_cols <- setdiff(names(dom_scores), "Morphology")
  B <- if (length(bio_cols)) do.call(cbind, dom_scores[bio_cols]) else NULL
  agg <- function(M) {
    if (is.null(M)) return(rep(NA_real_, nrow(bio)))
    if (is.null(dim(M))) M <- matrix(M, ncol = 1)
    n <- rowSums(!is.na(M)); s <- ifelse(n >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
    s[!is.finite(s)] <- NA_real_; s
  }
  list(domains = dom_scores,
       Morphology = morph, Biomarker = agg(B), Combined = agg(D))
}

results <- list()
for (mn in names(METHODS)) results[[mn]] <- build_index(METHODS[[mn]])

log_msg("")
log_msg("  --- domain level ---")
log_msg("  %-18s %5s  %7s %7s %7s", "domain", "n", "RAW", "ZSCORE", "PCTILE")
for (dom in names(DOMAINS)) {
  row <- c(); nn <- NA
  for (mn in names(METHODS)) {
    s <- results[[mn]]$domains[[dom]]
    if (is.null(s)) { row[mn] <- NA_real_; next }
    cc <- concordance(s, bio$Study_ID, is_ctrl, is_chal)
    row[mn] <- if (cc$n_eval >= 3) cc$prop else NA_real_
    if (mn == "RAW") nn <- cc$n_eval
  }
  fmt <- function(v) if (is.na(v)) "    -" else sprintf("%6.0f%%", 100 * v)
  log_msg("  %-18s %5s  %7s %7s %7s", dom, nn,
          fmt(row["RAW"]), fmt(row["ZSCORE"]), fmt(row["PERCENTILE"]))
}

log_msg("")
log_msg("  --- index level (p = binomial test vs a 50%% coin-flip null) ---")
summary_rows <- list()
for (v in c("Morphology", "Biomarker", "Combined")) {
  parts <- c()
  for (mn in names(METHODS)) {
    cc <- concordance(results[[mn]][[v]], bio$Study_ID, is_ctrl, is_chal)
    parts[[mn]] <- cc
  }
  log_msg("  GHI_%s", v)
  for (mn in names(METHODS)) {
    cc <- parts[[mn]]
    if (cc$n_eval == 0) { log_msg("    %-11s not evaluable", mn); next }
    p <- binom_p(cc$n_ok, cc$n_eval)
    log_msg("    %-11s %2d/%2d studies = %3.0f%%   p=%.3f%s",
            mn, cc$n_ok, cc$n_eval, 100 * cc$prop, p,
            if (!is.na(p) && p < 0.05) "  *significant" else "")
  }
  summary_rows[[v]] <- parts
}

# =============================================================================
log_msg("[4/4] Writing summary")
# =============================================================================
out <- c("=== Standardisation method comparison ===", "",
  "Question: is weak control-vs-challenge discrimination caused by the",
  "standardisation method, or inherent to the data?", "",
  "Metric: within-study concordance - the proportion of studies in which the",
  "untreated control scores healthier than that study's challenged arms.",
  "Null expectation 50%. Beneficial additives excluded from the contrast.", "",
  sprintf("Contrast: %d control rows vs %d challenged rows (%d additive/other excluded)",
          sum(is_ctrl), sum(is_chal), sum(!is_ctrl & !is_chal)),
  sprintf("Rows per study: median %.0f (range %d-%d)",
          stats::median(grp_per_study), min(grp_per_study), max(grp_per_study)),
  "", "--- index level ---")
for (v in names(summary_rows)) {
  out <- c(out, sprintf("  GHI_%s", v))
  for (mn in names(METHODS)) {
    cc <- summary_rows[[v]][[mn]]
    out <- c(out, if (cc$n_eval == 0) sprintf("    %-11s not evaluable", mn) else
      sprintf("    %-11s %2d/%2d = %3.0f%%   p=%.3f", mn, cc$n_ok, cc$n_eval,
              100 * cc$prop, binom_p(cc$n_ok, cc$n_eval)))
  }
}
out <- c(out, "", "--- how to read this ---",
  "If PERCENTILE clearly beats ZSCORE: the standard deviation estimated from",
  "  3-4 treatment-group means is unstable, and rank-based normalisation",
  "  recovers signal that z-scoring was destroying. Switch the pipeline to",
  "  percentile normalisation (precedent: Gut Microbiome Index, Oh 2022).",
  "",
  "If RAW clearly beats both: standardisation itself is the problem. Note RAW",
  "  is only usable for within-study comparison - raw units are not poolable",
  "  across studies - so this would point toward a within-study or paired",
  "  modelling strategy rather than a single pooled index.",
  "",
  "If all three are similar and near 50%: the limitation is the data, not the",
  "  method. Group-level means from heterogeneous published studies may simply",
  "  not carry enough signal. That is a legitimate finding, consistent with",
  "  Pietruska & Hauck (2023), where cross-study AUC fell to 0.42-0.68 in a",
  "  comparable 17-trial broiler meta-analysis. It would argue for treating",
  "  the GHI as a regression target rather than a standalone diagnostic, and",
  "  for resting the diagnostic claim on the qPCR validation instead.")

dir.create(file.path(R_OUT_DIR, "ghi"), recursive = TRUE, showWarnings = FALSE)
f <- file.path(R_OUT_DIR, "ghi", "standardisation_comparison.txt")
writeLines(out, f)
log_msg("  -> %s", f)
log_msg("=== DONE ===")
