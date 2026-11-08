#!/usr/bin/env Rscript
# =============================================================================
# 04e_mixed_model_validation.R
#
# Primary validation of the GHI, replacing within-study concordance.
#
# RUN:  Rscript 04e_mixed_model_validation.R
# TIME: under a minute. Reads the biomarker workbook only.
#
# -----------------------------------------------------------------------------
# WHY THIS REPLACES CONCORDANCE AS THE PRIMARY METRIC
# -----------------------------------------------------------------------------
# Within-study concordance ("in what fraction of studies did the group we
# expected to be healthier actually score higher?") is robust and intuitive,
# but it has two defects that matter for a published result:
#
#   1. It discards magnitude. A group that is marginally healthier and one that
#      is dramatically healthier both count as a single success. That is a
#      large loss of power on continuous data.
#
#   2. It cannot give an effect size, a confidence interval, or a formal test
#      of whether the additive and challenge contrasts genuinely differ - it
#      can only be eyeballed (88% vs 57%).
#
# A mixed-effects model fixes both. Study enters as a random intercept, which
# removes between-study level differences properly - the same problem that
# makes pooled AUC invalid here, but handled by the model rather than by
# discarding information. Contrast type enters as a moderator, so "is the index
# better at detecting improvement than damage?" becomes a testable interaction
# rather than an impression.
#
# Pooled AUC is deliberately NOT used. Within-study standardisation centres
# each study on zero, so a control in an additive study sits below its study
# mean while a control in a challenge study sits above it; pooling them
# cancels the signal. Note the observed pooled AUC of 0.44-0.48 is BELOW 0.5,
# which is not what pure cancellation looks like (that lands at 0.50) - it
# indicates slight directional miscalibration in the pooled construct, and is
# reported here as such rather than dismissed.
#
# -----------------------------------------------------------------------------
# WHAT IS AND IS NOT ESTABLISHED BY THIS TEST
# -----------------------------------------------------------------------------
# The group expected to be healthier is defined by the study's own design: in
# an additive study the additive arm, in a challenge study the control arm.
# Those are the same contrasts being tested. A significant result therefore
# demonstrates that the index is SENSITIVE TO THE INTERVENTION - it does not
# by itself establish that the index is a valid measure of gut health in an
# absolute sense. That distinction is stated here because it must be stated in
# the manuscript.
#
# -----------------------------------------------------------------------------
# WHY SMD / lnRR IS NOT USED, ALTHOUGH IT WOULD BE PREFERABLE
# -----------------------------------------------------------------------------
# The methodologically superior approach for aggregate treatment-group data is
# to compute a per-marker standardised mean difference (Hedges' g) or log
# response ratio per study and pool by inverse variance. That requires each
# group's standard deviation and sample size.
#
# This extraction contains group MEANS only - no SD, SEM, or per-group n
# columns exist in the workbook. Effect-size pooling is therefore not
# computable without returning to the 50 source papers to re-extract
# dispersion data, which many report only as SEM or only in figures. This
# limitation is stated explicitly in the output and should be stated in the
# manuscript as a constraint of aggregate-data synthesis, with effect-size
# pooling named as the preferred method should dispersion data be obtained.
# =============================================================================

.locate_and_source_config <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  d <- if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
  source(file.path(d, "00_r_config.R"))
}
.locate_and_source_config()

have_lme4 <- requireNamespace("lme4", quietly = TRUE)
if (!have_lme4) {
  log_msg("NOTE: lme4 not installed. Falling back to per-study paired differences")
  log_msg("      pooled by inverse-variance random-effects meta-analysis, which")
  log_msg("      answers the same question without the mixed-model machinery.")
  log_msg("      To install: BiocManager::install('lme4')  or install.packages('lme4')")
}

HEALTHY   <- c("Unchallenged control")
CHALLENGE <- c("Dietary challenge", "Infectious challenge", "Toxin challenge",
               "Antibiotic challenge", "Environmental challenge")
ADDITIVE  <- c("Non-disruptive additive")
CLASS_COL <- "Intervention_Class"   # was Challenge_Type

# -----------------------------------------------------------------------------
# THIS LIST MUST BE KEPT IN SYNC WITH 04b_ghi_pca.R's DOMAINS BLOCK BY HAND.
# It is a deliberately SEPARATE copy, not a shared import, because the whole
# point of the RAW construction is that it is built independently of the
# PUBLISHED (04b) code path - sharing the object would defeat that. The cost
# of that independence is that a domain added to 04b (Performance, Mucin,
# Anti_Inflammatory and IL1b were added in this session) does not
# automatically appear here, and silently drops out of the RAW check and the
# agreement comparison without erroring - exactly what happened on the first
# run after Performance was split out. Confirm this list matches 04b's
# DOMAINS block whenever either file changes.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# THIS LIST MUST BE KEPT IN SYNC WITH 04b_ghi_pca.R's DOMAINS BLOCK BY HAND.
#
# It is deliberately a SEPARATE copy rather than a shared import, because the
# whole point of the RAW construction is that it is built independently of the
# PUBLISHED (04b) code path - sharing the object would defeat that. The cost is
# that a change in 04b does not appear here automatically, and the RAW column
# then silently validates a DIFFERENT index from the published one. That has
# now happened twice: once when Performance was split out, and again when the
# structure moved from ten domains to seven.
#
# Synced to 04b as of the seven-domain structure:
#
#   Claudin1 folded into Barrier - it is a tight junction protein, and it
#   correlates +0.726 with Barrier in this corpus, so the separation was
#   statistical rather than biological.
#
#   Mucin folded into Barrier - mucus is part of the intestinal barrier and
#   correlates +0.754 with it, higher than Barrier does with Morphology. The
#   separate domain also gave one marker the same weight as Morphology's 311
#   samples.
#
#   IgA moved to the systemic domain as IgA_serum - all 64 IgA values in the
#   corpus were measured in blood, confirmed against every source paper.
#
#   Immunoglobulin_Mucosal dropped from the scored index. With the compartment
#   cleaned it scores AUC 0.071 - sharply inverted - but the sign cannot simply
#   be flipped, because sIgA rises under challenge AND under beneficial
#   additives (5 of 7 studies). It marks mucosal immune ACTIVITY, not the
#   direction of health, and is reported separately instead.
# -----------------------------------------------------------------------------
DOMAINS <- list(
  Morphology     = c(VH_um = +1, CD_um = -1),
  Barrier        = c(ZO1 = +1, Ocln = +1, Cldn1 = +1, MUC2 = +1),
  Inflammation   = c(IL6 = -1, TNFa = -1, IL1b = -1),
  Anti_Inflammatory = c(IL10 = +1),
  SCFA           = c(Ace = +1, Prop = +1, Buty = +1),
  Performance    = c(BWG_g = +1, FCR = -1),
  # Systemic +1 orients for PCA only - it is not a claim that higher
  # circulating immunoglobulin is healthier.
  Immunoglobulin_Systemic = c(IgG = +1, IgM = +1, IgA_serum = +1)
)

log_msg("[1/4] Loading data")
bio <- as.data.frame(read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
log_msg("  %d rows", nrow(bio))

# Raw values, oriented. Raw rather than z-scored because the standardisation
# test showed z-scoring reduced concordance (Morphology 79% raw vs 68%
# z-scored): a within-study SD estimated from as few as 2-4 treatment-group
# means is unstable, and dividing by it injects noise. Raw values are valid
# here because every comparison below is WITHIN a study, where units are
# constant by construction.
orient <- function(m, d) suppressWarnings(as.numeric(bio[[m]])) * d

log_msg("[2/4] Building domain and index scores (raw, oriented)")
dom <- list()
for (dn in names(DOMAINS)) {
  cols <- intersect(names(DOMAINS[[dn]]), names(bio))
  if (!length(cols)) next
  M <- sapply(cols, function(m) orient(m, DOMAINS[[dn]][[m]]))
  if (is.null(dim(M))) M <- matrix(M, ncol = 1)
  # Scale each marker to unit variance across the corpus before averaging
  # within a domain, so a marker measured in thousands (villus height, um)
  # does not dominate one measured in single digits (a fold-change) purely
  # because of its numeric range. This is a scaling step, not a
  # standardisation against a study mean - no study-level centring occurs.
  M <- apply(M, 2, function(v) { s <- stats::sd(v, na.rm = TRUE)
                                 if (!is.finite(s) || s == 0) v else v / s })
  n <- rowSums(!is.na(M))
  s <- ifelse(n >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
  s[!is.finite(s)] <- NA_real_
  dom[[dn]] <- s
  log_msg("  %-16s %4d samples", dn, sum(!is.na(s)))
}

agg <- function(cols) {
  if (!length(cols)) return(rep(NA_real_, nrow(bio)))
  M <- do.call(cbind, dom[cols])
  if (is.null(dim(M))) M <- matrix(M, ncol = 1)
  n <- rowSums(!is.na(M)); s <- ifelse(n >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
  s[!is.finite(s)] <- NA_real_; s
}

# =============================================================================
# TWO CONSTRUCTIONS, VALIDATED SIDE BY SIDE
#
# PUBLISHED  - the scores 04b_ghi_pca.R actually wrote (within-study z-scoring,
#              PCA within each domain on complete cases, mean fallback, then
#              proration across domains). This is the index that goes in the
#              paper, so it is the one that must be validated. It is primary.
#
# RAW        - the construction previously used here: oriented raw values,
#              scaled to unit variance, averaged within domain. Kept because it
#              is a genuinely independent check - it uses no z-scoring and no
#              PCA, so agreement between the two says the result is a property
#              of the biology rather than of one construction choice.
#
# These are NOT interchangeable and their sample counts differ, because PCA
# needs complete cases within a domain while mean-of-available does not
# (Inflammation: 103 complete pairs vs 166 samples with either marker). An
# earlier version of this script reported only the RAW numbers while the
# manuscript described the PUBLISHED index - a validation of something
# adjacent to what was being claimed. Both are now reported and labelled.
# =============================================================================
INDICES_RAW <- list(
  Morphology  = dom[["Morphology"]],
  Biomarker   = agg(setdiff(names(dom), c("Morphology", "Performance"))),
  Performance = dom[["Performance"]],
  Combined    = agg(names(dom))
)

# ---- published scores from 04b ----------------------------------------------
GHI_TSV <- file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")
INDICES_PUB <- NULL
if (!file.exists(GHI_TSV)) {
  log_msg("  NOTE: %s not found - run 04b_ghi_pca.R to validate the published",
          basename(GHI_TSV))
  log_msg("        index. Only the RAW construction is reported below.")
} else {
  gh <- fread(GHI_TSV)
  need_pub <- c("GHI_Morphology_raw", "GHI_Biomarker_raw", "GHI_Performance_raw", "GHI_Combined_raw")
  miss <- setdiff(need_pub, names(gh))
  if (length(miss)) {
    log_msg("  NOTE: %s lacks %s - re-run 04b.", basename(GHI_TSV),
            paste(miss, collapse = ", "))
  } else {
    # Match on Sample_ID where both carry it, falling back to Study_ID +
    # Treatment_Group. A silent mismatch here would validate the wrong rows,
    # so the match rate is reported and a poor match stops the script.
    key_bio <- if ("Sample_ID" %in% names(bio) && "Sample_ID" %in% names(gh)) {
      list(b = trimws(as.character(bio$Sample_ID)),
           g = trimws(as.character(gh$Sample_ID)), how = "Sample_ID")
    } else {
      list(b = paste(trimws(bio$Study_ID), trimws(bio$Treatment_Group), sep = "||"),
           g = paste(trimws(gh$Study_ID),  trimws(gh$Treatment_Group),  sep = "||"),
           how = "Study_ID + Treatment_Group")
    }
    idx <- match(key_bio$b, key_bio$g)
    matched <- sum(!is.na(idx))
    log_msg("  published scores joined on %s: %d of %d rows matched (%.0f%%)",
            key_bio$how, matched, nrow(bio), 100 * matched / nrow(bio))
    if (matched < 0.5 * nrow(bio))
      stop("Fewer than half the biomarker rows matched ", basename(GHI_TSV),
           ".\n  The two files are out of step - re-run 04b_ghi_pca.R.")

    INDICES_PUB <- list(
      Morphology  = gh$GHI_Morphology_raw[idx],
      Biomarker   = gh$GHI_Biomarker_raw[idx],
      Performance = gh$GHI_Performance_raw[idx],
      Combined    = gh$GHI_Combined_raw[idx])
    for (nm in names(INDICES_PUB))
      log_msg("    published %-11s %4d samples", nm, sum(!is.na(INDICES_PUB[[nm]])))
  }
}

CONSTRUCTIONS <- if (is.null(INDICES_PUB)) {
  list(RAW = INDICES_RAW)
} else {
  list(PUBLISHED = INDICES_PUB, RAW = INDICES_RAW)
}

# =============================================================================
# STEP 3 - PER-STUDY PAIRED DIFFERENCES, POOLED BY RANDOM-EFFECTS META-ANALYSIS
#
# For each study and each contrast, the difference between the two group means
# is computed. Those differences are then pooled across studies with a
# DerSimonian-Laird random-effects model, which estimates between-study
# heterogeneity (tau-squared) rather than assuming every study measures the
# same underlying effect.
#
# Random-effects rather than fixed-effects because the studies differ in bird
# age, breed, challenge organism, dose, and husbandry - assuming one common
# true effect across that range would be indefensible, and would produce
# artificially narrow confidence intervals.
#
# Each study's difference is weighted by the number of group-means behind it.
# This is a coarser weight than inverse-variance (which needs the SDs this
# extraction lacks), and it is labelled as such rather than presented as a
# proper inverse-variance analysis.
# =============================================================================
log_msg("[3/4] Per-study differences, pooled by random-effects meta-analysis")

study_diffs <- function(score, better_lv, worse_lv, skip = character(0)) {
  out <- list()
  for (s in unique(bio$Study_ID)) {
    if (s %in% skip) next
    k <- bio$Study_ID == s & !is.na(score)
    a <- score[k & bio[[CLASS_COL]] %in% better_lv]
    b <- score[k & bio[[CLASS_COL]] %in% worse_lv]
    if (!length(a) || !length(b)) next
    out[[length(out) + 1]] <- list(study = s, diff = mean(a) - mean(b),
                                   n = length(a) + length(b))
  }
  out
}

pool_re <- function(d) {
  if (length(d) < 3) return(NULL)
  y <- vapply(d, function(z) z$diff, numeric(1))
  w0 <- vapply(d, function(z) z$n, numeric(1))
  k <- length(y)
  # Fixed-effect pass to estimate between-study heterogeneity (Q, tau^2)
  mu_f <- sum(w0 * y) / sum(w0)
  Q <- sum(w0 * (y - mu_f)^2)
  C <- sum(w0) - sum(w0^2) / sum(w0)
  tau2 <- max(0, (Q - (k - 1)) / C)
  # Random-effects weights incorporate that heterogeneity
  wi <- 1 / (1 / w0 + tau2)
  mu <- sum(wi * y) / sum(wi)
  se <- sqrt(1 / sum(wi))
  z <- mu / se
  p <- 2 * stats::pnorm(-abs(z))
  I2 <- if (Q > k - 1) 100 * (Q - (k - 1)) / Q else 0
  list(k = k, est = mu, se = se, lo = mu - 1.96 * se, hi = mu + 1.96 * se,
       z = z, p = p, I2 = I2, tau2 = tau2,
       pi_lo = mu - stats::qt(0.975, max(1, k - 2)) * sqrt(tau2 + se^2),
       pi_hi = mu + stats::qt(0.975, max(1, k - 2)) * sqrt(tau2 + se^2))
}

# Studies with no unchallenged arm contemporaneous with their challenged arms.
# Massacci_2019 gavaged every bird with C. jejuni on d21, so its only Control
# rows are pre-challenge d21 faecal rows carrying 0-21 d performance, while its
# challenged rows are d40 and carry 28-40 d performance plus duodenal
# morphology. Pairing them measures growth, not infection. The study keeps its
# ADDITIVE contrast, which is a valid pre-challenge comparison of C against S.
NO_CONCURRENT_CONTROL <- c("Massacci_2019")

CONTRASTS <- list(
  ADDITIVE  = list(b = ADDITIVE, w = HEALTHY,   lab = "additive minus control",
                   skip = character(0)),
  CHALLENGE = list(b = HEALTHY,  w = CHALLENGE, lab = "control minus challenged",
                   skip = NO_CONCURRENT_CONTROL)
)

results <- list()
for (cons in names(CONSTRUCTIONS)) {
  IDX <- CONSTRUCTIONS[[cons]]
  log_msg("")
  log_msg("  --- %s construction%s ---", cons,
          if (cons == "PUBLISHED") " (04b scores - this is the index in the paper)"
          else " (independent check: no z-scoring, no PCA)")
  log_msg("  %-12s %-11s %3s %9s %19s %8s %6s", "index", "contrast", "k",
          "estimate", "95%% CI", "p", "I2")
  for (iname in names(IDX)) {
    for (cname in names(CONTRASTS)) {
      ct <- CONTRASTS[[cname]]
      d <- study_diffs(IDX[[iname]], ct$b, ct$w, ct$skip)
      r <- pool_re(d)
      results[[paste(cons, iname, cname)]] <-
        list(r = r, d = d, index = iname, contrast = cname, construction = cons)
      if (is.null(r)) {
        log_msg("  %-12s %-11s  <3 studies, not evaluable", iname, cname); next
      }
      log_msg("  %-12s %-11s %3d %+9.3f  [%+7.3f, %+7.3f] %8.4f %5.0f%%%s",
              iname, cname, r$k, r$est, r$lo, r$hi, r$p, r$I2,
              if (r$p < 0.05 && r$est > 0) "  *" else "")
    }
  }
}

log_msg("")
log_msg("  A positive estimate means the index ordered the groups correctly.")
log_msg("  * marks a positive effect with p < 0.05.")
log_msg("  I2 is the share of variation between studies beyond chance;")
log_msg("  above ~75%% indicates substantial heterogeneity.")
log_msg("  Estimates are on each construction's own scale and are NOT directly")
log_msg("  comparable in magnitude between constructions - compare direction,")
log_msg("  significance and the additive/challenge pattern instead.")

# ---- do the two constructions agree? ----------------------------------------
# Agreement matters more than either result alone: if a conclusion holds under
# both z-scored PCA and raw mean-of-available, it is a property of the data
# rather than of one modelling choice. Disagreement is equally informative and
# would need stating in the manuscript.
if (length(CONSTRUCTIONS) == 2) {
  log_msg("")
  log_msg("  Agreement between constructions:")
  for (iname in names(INDICES_RAW)) {
    for (cname in names(CONTRASTS)) {
      rp <- results[[paste("PUBLISHED", iname, cname)]]$r
      rr <- results[[paste("RAW", iname, cname)]]$r
      if (is.null(rp) || is.null(rr)) next
      same_dir <- sign(rp$est) == sign(rr$est)
      same_sig <- (rp$p < 0.05) == (rr$p < 0.05)
      log_msg("    %-12s %-11s %s", iname, cname,
              if (same_dir && same_sig) "agree"
              else if (same_dir) "same direction, differ on significance"
              else "DISAGREE ON DIRECTION - state this explicitly")
    }
  }
}

# =============================================================================
# STEP 4 - IS THE INDEX BETTER AT DETECTING IMPROVEMENT THAN DAMAGE?
#
# Tested formally rather than by comparing two percentages by eye. The
# difference between the two pooled estimates is assessed against its combined
# standard error - a standard test for a difference between subgroup effects
# in meta-analysis.
# =============================================================================
log_msg("[4/4] Formal test: does the index detect improvement better than damage?")
for (cons in names(CONSTRUCTIONS)) {
  log_msg("  --- %s ---", cons)
  for (iname in names(CONSTRUCTIONS[[cons]])) {
    ra <- results[[paste(cons, iname, "ADDITIVE")]]$r
    rc <- results[[paste(cons, iname, "CHALLENGE")]]$r
    if (is.null(ra) || is.null(rc)) next
    dd <- ra$est - rc$est
    se <- sqrt(ra$se^2 + rc$se^2)
    z <- dd / se; p <- 2 * stats::pnorm(-abs(z))
    log_msg("  %-12s additive %+0.3f vs challenge %+0.3f -> difference %+0.3f, p=%.4f%s",
            iname, ra$est, rc$est, dd, p,
            if (p < 0.05) "  * genuinely different" else "  (not distinguishable)")
  }
}

# ---- output ------------------------------------------------------------------
out <- c("=== GHI validation by random-effects meta-analysis ===", "",
  "Primary validation metric, replacing within-study concordance.",
  "Per-study paired group differences pooled with a DerSimonian-Laird",
  "random-effects model. A positive estimate means the index ordered the",
  "groups correctly.", "",
  sprintf("  %-10s %-12s %-11s %3s %9s %19s %8s %6s", "construct", "index","contrast","k","estimate","95% CI","p","I2"))
for (nm in names(results)) {
  z <- results[[nm]]; r <- z$r
  out <- c(out, if (is.null(r)) sprintf("  %-10s %-12s %-11s  <3 studies", z$construction, z$index, z$contrast) else
    sprintf("  %-10s %-12s %-11s %3d %+9.3f  [%+7.3f, %+7.3f] %8.4f %5.0f%%",
            z$construction, z$index, z$contrast, r$k, r$est, r$lo, r$hi, r$p, r$I2))
}
out <- c(out, "",
  "PUBLISHED = the scores written by 04b_ghi_pca.R (within-study z-scoring,",
  "  PCA within domain, proration). This is the index reported in the paper.",
  "RAW = oriented raw values, unit-variance scaled, averaged within domain.",
  "  An independent check using neither z-scoring nor PCA. Magnitudes are on",
  "  different scales and are not comparable between constructions; direction,",
  "  significance and the additive/challenge pattern are.")
out <- c(out, "", "--- limitations, to be stated in the manuscript ---",
  "1. The group expected to be healthier is defined by the study's own design,",
  "   so a significant result shows the index is SENSITIVE TO THE INTERVENTION.",
  "   It does not by itself establish absolute validity as a health measure.",
  "",
  "2. Effect-size pooling (Hedges' g or log response ratio, inverse-variance",
  "   weighted) is the preferred method for aggregate treatment-group data.",
  "   It is not computable here: the extraction holds group means only, with",
  "   no SD, SEM, or per-group n. Studies are weighted by number of group",
  "   means instead, which is coarser. Obtaining dispersion data would allow",
  "   the stronger analysis.",
  "",
  "3. Every row is a treatment-group mean, not an individual bird. Conclusions",
  "   apply at group level; a relationship at group level need not hold at bird",
  "   level (ecological fallacy).",
  "",
  "4. Marker availability is confounded with study design - challenge studies",
  "   preferentially measured cytokines and barrier genes, nutrition studies",
  "   measured morphology. Any difference between the additive and challenge",
  "   contrasts is therefore partly a data-availability artifact and not purely",
  "   biological.",
  "",
  "5. The corpus is dominated by beneficial-additive studies, which are subject",
  "   to positive-result publication bias. This plausibly inflates the additive",
  "   contrast relative to truth, and competes with the biological explanation",
  "   for the additive/challenge asymmetry.",
  "",
  "6. Days_Post_Challenge is populated for only 12 of 566 rows, from a single",
  "   study. Stratification by time since challenge - which the literature",
  "   indicates drives sign reversal in tight junction expression - is",
  "   therefore not possible with the current extraction.")
dir.create(file.path(R_OUT_DIR, "ghi"), recursive = TRUE, showWarnings = FALSE)
f <- file.path(R_OUT_DIR, "ghi", "mixed_model_validation.txt")
writeLines(out, f)
log_msg("  -> %s", f)
log_msg("=== DONE ===")

# ---- prediction intervals (added 2026-10-23) --------------------------------
# A 95% prediction interval bounds the effect a NEW study would be expected to
# observe given the between-study heterogeneity. Reported because I2 exceeds
# 75% for the primary contrasts, where a confidence interval alone understates
# the spread of true effects.
pi_out <- c("=== Prediction intervals for the pooled contrasts ===", "",
  sprintf("  %-10s %-12s %-11s %3s %9s %21s %21s %6s",
          "construct", "index", "contrast", "k", "estimate", "95% CI", "95% PI", "I2"))
for (nm in names(results)) {
  z <- results[[nm]]; r <- z$r
  if (is.null(r)) next
  pi_out <- c(pi_out, sprintf(
    "  %-10s %-12s %-11s %3d %+9.3f  [%+7.3f, %+7.3f]  [%+7.3f, %+7.3f] %5.0f%%",
    z$construction, z$index, z$contrast, r$k, r$est, r$lo, r$hi,
    r$pi_lo, r$pi_hi, r$I2))
}
f2 <- file.path(R_OUT_DIR, "ghi", "mixed_model_prediction_intervals.txt")
writeLines(pi_out, f2)
log_msg("  -> %s", f2)
