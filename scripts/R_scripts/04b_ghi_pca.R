#!/usr/bin/env Rscript
# =============================================================================
# 04b_ghi_pca.R
#
# Constructs the Gut Health Index by PCA, as specified in the registered
# protocol, with 0-100 rescaling and ROC-derived interpretive cutoffs.
#
# HOW TO RUN: Rscript 04b_ghi_pca.R  (or open in RStudio and Source)
# Requires 00_r_config.R alongside it. Run after the QIIME 2 pipeline.
#
# -----------------------------------------------------------------------------
# RELATIONSHIP TO THE REGISTERED PROTOCOL
# -----------------------------------------------------------------------------
# PROSPERO and the Methods specify: "a composite Gut Health Index derived by
# principal component analysis of the standardised biomarker matrix." PCA
# requires complete cases - every marker present for every sample it scores.
#
# In this corpus that is achievable WITHIN a biological domain but not ACROSS
# domains. Complete-case counts, measured directly from the data:
#
#     Morphology (VH + CD)            311
#     Barrier (ZO-1, occludin, CLDN1) 111
#     Inflammation (IL-6 + TNF-a)     103
#     SCFA (acetate/propionate/butyrate) 74
#     Immunoglobulin (all four)         0   [now split mucosal / systemic]
#     All non-morphology markers        0
#     All markers                       0
#
# The zeros are not a data-entry problem: which biomarkers a study measured is
# determined by what that study was designed to test (challenge studies assay
# cytokines, nutrition studies assay morphology), so no sample carries the full
# panel. Complete-case PCA across domains is therefore not merely underpowered,
# it is undefined.
#
# This script applies PCA where the protocol's method is computable - within
# each domain - and aggregates domain scores by mean where it is not. The
# deviation is confined to the cross-domain step, is reported explicitly in the
# output summary, and each domain's PC1 variance-explained and loadings are
# printed so the composite can be inspected rather than taken on trust.
#
# A domain with too few complete cases for PCA falls back to the mean of its
# available z-scored markers, and is labelled as such in the output.
#
# -----------------------------------------------------------------------------
# WHY VH:CD IS EXCLUDED FROM THE MORPHOLOGY PCA
# -----------------------------------------------------------------------------
# VH:CD is the arithmetic ratio of villus height to crypt depth. Entering it
# alongside its own numerator and denominator gives PCA a variable that is a
# deterministic function of the other two, which inflates the shared component
# and distorts the loadings. VH and CD alone are used (n=311, also more samples
# than the three-marker set at n=295); VH:CD is retained in the output as a
# separately reported descriptive measure.
# =============================================================================

.locate_and_source_config <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    this_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else if (requireNamespace("rstudioapi", quietly = TRUE) &&
             isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    this_dir <- if (!is.null(ctx) && nzchar(ctx$path)) dirname(normalizePath(ctx$path)) else normalizePath(".")
  } else {
    this_dir <- normalizePath(".")
  }
  source(file.path(this_dir, "00_r_config.R"))
}
.locate_and_source_config()

# =============================================================================
# CONFIGURATION
# =============================================================================

# Minimum complete cases for a domain to be scored by PCA rather than by the
# mean-of-available fallback. Below roughly 30, PC1 loadings are unstable
# enough that the fallback is the more honest choice.
# A domain's PC1 must correlate at least this strongly with the mean of its
# oriented markers, otherwise PC1 is measuring something other than the
# construct and the mean is used instead. See the guard in step 4.
#
# Settable from the environment ONLY so that the pre-correction behaviour can
# be reproduced for comparison:
#
#   PC1_ALIGN_MIN=0 Rscript 04b_ghi_pca.R
#
# At zero the guard never fires and the Morphology domain reverts to computing
# the SUM of villus height and crypt depth rather than their contrast - that
# is, intestinal size rather than gut health. Measured against the published
# scores, the domain then correlates +1.000 with the size axis and +0.017 with
# the intended health contrast, and the morphology additive effect falls from
# +0.346 to +0.070.
#
# This setting exists to document the difference, not to be used. Any result
# produced under it describes the earlier, incorrect index.
PC1_ALIGN_MIN <- {
  e <- Sys.getenv("PC1_ALIGN_MIN", "")
  v <- if (nzchar(e)) suppressWarnings(as.numeric(e)) else 0.30
  if (!is.finite(v)) 0.30 else v
}
if (PC1_ALIGN_MIN < 0.30)
  log_msg("  *** PC1_ALIGN_MIN = %.2f - alignment guard weakened or disabled.", PC1_ALIGN_MIN)
if (PC1_ALIGN_MIN <= 0)
  log_msg("  *** Morphology will revert to the size axis. Results are not the corrected index.")

MIN_N_FOR_PCA <- 30

# PC1 must explain at least this share of variance for the domain's markers to
# be considered as sharing a common signal. Reported either way, never used to
# silently drop a domain.
PC1_VARIANCE_FLOOR <- 0.35

# Claudin-1 direction switch. +1 = higher is healthier (current default,
# supported by this corpus); -1 = higher marks damage (supported by the
# poultry challenge literature). Set from the environment so both can be run
# without editing the file:  CLDN1_DIR=-1 Rscript 04b_ghi_pca.R
CLDN1_DIR <- {
  e <- Sys.getenv("CLDN1_DIR", "")
  if (nzchar(e)) as.numeric(e) else +1
}
if (!CLDN1_DIR %in% c(-1, 1)) stop("CLDN1_DIR must be +1 or -1, got ", CLDN1_DIR)

# -----------------------------------------------------------------------------
# HOW CLAUDIN-1 IS HANDLED
#
# Three configurations, set by CLDN1_MODE:
#
#   "drop"     excluded from the index entirely            (default)
#   "separate" scored as its own domain
#   "merge"    scored inside Barrier alongside ZO-1 and occludin
#
# DROP is the default, on a prior interpretability argument rather than on
# performance. The direction of effect for claudin-1 is contested: the poultry
# challenge literature reports it RISING during active damage while ZO-1 and
# occludin fall, whereas in this corpus it correlates +0.726 with Barrier
# (n = 111), the strongest pairing in the domain matrix. Most arms here are
# additive trials or late-sampled challenges, where the compensatory-repair
# phase has every junction transcript moving together, so the disagreement is
# most plausibly one of sampling window rather than of biology.
#
# A marker whose sign cannot be defended should not carry weight in an index
# whose every other marker has a stable, justifiable direction. Excluding it
# removes the only contested direction in the panel, and removes the awkward
# case of a domain built on a single marker.
#
# The alternatives are retained so all three can be reported:
#   CLDN1_MODE=separate Rscript 04b_ghi_pca.R
#   CLDN1_MODE=merge    Rscript 04b_ghi_pca.R
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ALTERNATIVE CONSTRUCTIONS
#
# MORPH_MODE controls how intestinal morphology enters the index.
#   "markers"  villus height and crypt depth as two oriented markers (default)
#   "ratio"    the study-reported VH:CD ratio as a single marker
#
# The ratio is the measure the field itself reports and is size-invariant by
# construction, so it sidesteps the problem that villus height and crypt depth
# both scale with body mass. It is available on 295 rows against 311 carrying
# both raw measurements, so it costs a little coverage.
#
# AGGREGATION controls how markers become an index.
#   "domain"   markers grouped into domains, each domain scored, domain scores
#              averaged (default). Every facet carries equal weight regardless
#              of how many markers a study measured for it.
#   "flat"     every oriented marker averaged directly, no domains. Simpler,
#              but a sample carrying four barrier markers and one morphology
#              marker becomes barrier-dominated, which is the imbalance the
#              domain structure exists to prevent.
#
# Both default to the registered construction. The alternatives exist so the
# choice can be reported as a sensitivity analysis rather than asserted.
# -----------------------------------------------------------------------------
MORPH_MODE  <- { e <- tolower(Sys.getenv("MORPH_MODE", ""));  if (e == "ratio") "ratio" else "markers" }
AGGREGATION <- { e <- tolower(Sys.getenv("AGGREGATION", "")); if (e == "flat")  "flat"  else "domain" }

CLDN1_MODE <- {
  e <- tolower(Sys.getenv("CLDN1_MODE", ""))
  if (e %in% c("drop", "separate", "merge")) e else "merge"
}
MERGE_CLDN1 <- CLDN1_MODE == "merge"

# -----------------------------------------------------------------------------
# MUC2: INSIDE BARRIER BY DEFAULT
#
# MUC2 was scored as its own domain on the grounds that a secreted mucus gel is
# a different structure from a tight-junction protein. That distinction is real,
# but it does not survive the same two arguments already applied to claudin-1
# in this script:
#
#   - Mucus IS part of the intestinal barrier. The poultry literature reports
#     it that way, grouping mucin 2 with ZO-1 and occludin as barrier outcomes
#     while holding cytokines separate. The separation was a judgment about
#     structure, not about function.
#
#   - A one-marker domain carried the same weight in the domain average as
#     Morphology, which rests on 311 samples. That is the identical
#     over-weighting objection used to fold claudin-1 in, left unapplied here.
#
# The corpus agrees: Mucin correlates +0.754 with Barrier, higher than Barrier
# correlates with Morphology (+0.472) and far higher than Mucin with Morphology
# (+0.316).
#
# Tested before adopting: merging changed the additive contrast from +0.764 to
# +0.755 and the challenge contrast from +0.813 to +0.794, both differences far
# inside the between-study heterogeneity. The merge is adopted because it
# removes an arbitrary split, not because it improves anything.
#
# MERGE_MUC2=no restores the separate domain so both can be reported.
# -----------------------------------------------------------------------------
MERGE_MUC2 <- { e <- tolower(Sys.getenv("MERGE_MUC2", "")); e != "no" }

# sIgA is excluded from the scored index by default - see the note on the
# Immunoglobulin_Mucosal domain below. SCORE_SIGA=yes puts it back.
SCORE_SIGA <- { e <- tolower(Sys.getenv("SCORE_SIGA", "")); e == "yes" }

# Direction of effect: +1 means a higher raw value indicates BETTER gut health.
# Every marker is multiplied by this before analysis so that, after orientation,
# higher always means healthier. A wrong entry here does not raise an error -
# it silently inverts that marker - so these are listed explicitly.
DOMAINS <- list(
  Morphology = if (MORPH_MODE == "ratio") c(
    VH_CD = +1    # the study-reported ratio; size cancels in the quotient
  ) else c(
    VH_um = +1,   # taller villi = greater absorptive surface
    CD_um = -1    # deeper crypts = higher epithelial turnover / damage response
  ),
  # -------------------------------------------------------------------------
  # CLAUDIN-1: INSIDE BARRIER BY DEFAULT
  #
  # It was originally scored as its own domain, on the grounds that the poultry
  # literature reports claudin-1 RISING under challenge while ZO-1 and occludin
  # fall - graded Eimeria upregulates CLDN1 with occludin unchanged, necrotic
  # enteritis strains raise CLDN1 while lowering OCLN and TJP1, AFB1 raises
  # jejunal CLDN1. On that evidence, averaging the three would cancel signal.
  #
  # THIS CORPUS DOES NOT REPRODUCE THAT PATTERN. In 04g, Claudin1 correlated
  # +0.726 with Barrier (n = 111), the strongest pairing in the whole matrix.
  # If claudin-1 were opposing the other two here, that would be negative. The
  # likely reason is sampling window: most arms are additive trials or
  # late-sampled challenges, where the compensatory-repair phase has all
  # junction transcripts moving together.
  #
  # Three further reasons to merge:
  #   - claudin-1 IS a tight junction protein; the separation was statistical,
  #     not biological
  #   - a one-marker domain carried the same weight in the domain average as
  #     Morphology, which rests on 323 samples - that over-weights one marker
  #   - six domains are easier to justify than seven built on one exception
  #
  # The original evidence for separating it came from concordance analysis,
  # which was later replaced as an invalid metric. It is not evidence any more.
  #
  # MERGE_CLDN1=no restores the separate domain, so both can be reported.
  # -------------------------------------------------------------------------
  Barrier = c(
    c(ZO1 = +1, Ocln = +1),                       # tight junction integrity
    if (MERGE_CLDN1) c(Cldn1 = CLDN1_DIR) else NULL,
    if (MERGE_MUC2)  c(MUC2  = +1)        else NULL   # secreted mucus layer
  ),
  # IL1b joins IL6/TNFa - same direction (pro-inflammatory, higher = worse),
  # so no mixed-orientation risk. Previously extracted but never scored;
  # 139 rows / 21 studies in the source workbook.
  Inflammation = c(
    IL6  = -1,    # pro-inflammatory
    TNFa = -1,
    IL1b = -1
  ),
  # IL10 is anti-inflammatory - the opposite direction from IL6/TNFa/IL1b, and
  # kept as its own domain rather than folded in with a flipped sign, for the
  # same reason claudin-1 and systemic Ig are separated: a marker whose
  # biological role differs from its domain-mates should not be allowed to
  # cancel their signal inside one PCA. Previously extracted but never
  # scored; 114 rows / 19 studies.
  Anti_Inflammatory = c(
    IL10 = +1    # anti-inflammatory cytokine; higher = better-regulated response
  ),
  SCFA = c(
    Ace  = +1,    # microbial fermentation products
    Prop = +1,
    Buty = +1
  ),
  # MUC2: mucin, the structural material of the protective mucus layer.
  # Scored as its own domain rather than folded into Barrier, on the same
  # logic as claudin-1 - it is a distinct barrier component (secreted gel
  # layer, not a tight-junction protein) and its own domain keeps that
  # distinction visible rather than assumed. Previously extracted but never
  # scored; 92 rows / 15 studies.
  Mucin = if (MERGE_MUC2) NULL else c(
    MUC2 = +1    # more mucin = a more intact protective barrier
  ),
  # PERFORMANCE: body weight gain and feed conversion ratio. These were in
  # the original registered biomarker panel and were extracted (BWG_g: 125
  # rows / 16 studies; FCR: 238 rows / 31 studies - among the best-covered
  # markers in the whole panel) but never entered any domain. FCR is oriented
  # -1 because a LOWER feed conversion ratio means the bird converted feed to
  # body mass more efficiently, i.e. healthier.
  Performance = c(
    BWG_g = +1,
    FCR   = -1
  ),
  # The immunoglobulin panel is split into mucosal and systemic rather than
  # scored as one domain. This is the same treatment given to claudin-1 above
  # and for the same reason: markers whose direction of effect is not stable
  # across contexts contaminate a domain whose other markers are stable.
  #
  # The split is driven by a measurement, not by argument. In the domain
  # correlation diagnostic (04g), the combined Immunoglobulin domain ran
  # NEGATIVE against Morphology (r = -0.247, n = 68) while sitting at zero
  # against Barrier (+0.015, n = 36) and Claudin1 (+0.028, n = 26) - it neither
  # joined the coherent biomarker block nor tracked morphology. Every other
  # inter-domain correlation on usable n was positive. That is the signature of
  # markers pulling in opposite directions inside one domain.
  #
  # -------------------------------------------------------------------------
  # MUCOSAL sIgA: MEASURED AND REPORTED, BUT NOT SCORED
  #
  # This domain was "sIgA and IgA, higher is better". Two corrections and one
  # test removed the basis for that.
  #
  # First, the IgA in it was serum. All 64 IgA values in the corpus were
  # measured in blood, confirmed against every source paper. They now sit in
  # the systemic domain. A further 15 values in the sIgA column itself were
  # serum IgA from three studies that measured no secretory IgA at all
  # (Jansseune_2024, Han_2023, Chen_2024), leaving 55 genuinely mucosal rows
  # across 8 studies - all jejunal or ileal homogenate.
  #
  # Second, with the compartment clean the domain scores AUC 0.071. That is
  # not weak, it is inverted: sIgA is one of the sharpest discriminators in
  # the corpus, pointing the opposite way to the +1 sign.
  #
  # Third - and this is why the sign is not simply flipped - sIgA rises under
  # BOTH conditions. Challenged arms: 2 studies higher, 0 lower. Additive
  # arms: 5 studies higher, 2 lower. Flipping would score the challenge
  # contrast at 0.929 while marking probiotic-supplemented birds as unhealthy,
  # because their sIgA rose too. Neither sign is right.
  #
  # The mechanism is documented: secretory IgA reflects a mature, well
  # regulated mucosal barrier in a healthy bird AND a defensive response to
  # infection in a sick one. The standing recommendation is to interpret it
  # alongside performance, lesions and inflammatory markers rather than alone.
  # A composite index scores each domain alone, which is exactly what sIgA
  # cannot support.
  #
  # It is therefore excluded from the scored index and reported separately.
  # The cost is nil: combined AUC is 0.646 with it and 0.651 without.
  #
  # SCORE_SIGA=yes restores it, so both can be reported.
  # -------------------------------------------------------------------------
  Immunoglobulin_Mucosal = if (SCORE_SIGA) c(
    sIgA = +1
  ) else NULL,
  # SYSTEMIC (IgG, IgM): circulating immunoglobulin rises in RESPONSE to
  # antigenic exposure. In a challenge study, elevated systemic Ig marks the
  # insult, not health. It is scored as its own domain and read separately
  # rather than assigned a sign the corpus does not support.
  #
  # The +1 below is retained ONLY so the domain is oriented consistently for
  # PCA; it is not a claim that higher is healthier. Do not fold this domain
  # into a composite without stating that assumption explicitly, and run the
  # downstream taxon selection both ways as a sensitivity analysis.
  Immunoglobulin_Systemic = c(
    IgG  = +1,    # REVIEW: direction unresolved - see note above
    IgM  = +1,
    # Serum IgA, moved here from the mucosal domain by fix_iga_compartment.R.
    # All eleven studies contributing IgA measured it in blood serum, each
    # confirmed by a direct methods quote; none described it as secretory.
    # Four of them measured mucosal sIgA separately, in tissue homogenates,
    # which stays in the mucosal domain. Leaving serum IgA there put roughly a
    # third of the "mucosal" domain in the wrong compartment and plausibly
    # explains why it discriminated worse than the systemic domain
    # (AUC 0.588 vs 0.685).
    IgA_serum = +1
  )
)

# domains switched off by MERGE_* flags come through as NULL; drop them
DOMAINS <- Filter(function(z) !is.null(z) && length(z) > 0, DOMAINS)

# When claudin-1 is not merged into Barrier it is restored as its own domain,
# so the two configurations differ only in where Cldn1 sits.
if (CLDN1_MODE == "separate") {
  DOMAINS <- append(DOMAINS,
                    list(Claudin1 = c(Cldn1 = CLDN1_DIR)),
                    after = which(names(DOMAINS) == "Barrier"))
}
log_msg("  claudin-1: %s", switch(CLDN1_MODE,
        drop     = "EXCLUDED from the index (direction contested)",
        separate = "scored as its own domain",
        merge    = "inside the Barrier domain"))
log_msg("  domains: %d (%s)", length(DOMAINS), paste(names(DOMAINS), collapse = ", "))


MORPHOLOGY_DOMAIN <- "Morphology"
# "Biomarker" was originally every non-morphology domain, defined when every
# domain in that bucket was a direct physiological measurement (barrier,
# inflammation, immunoglobulins, SCFA, mucin). Performance (BWG, FCR) is a
# zootechnical OUTCOME, not a biomarker - a downstream result of gut health
# among other things, not a measurement of a gut process itself. Left inside
# "Biomarker" it was a mislabelled category, not just an odd inclusion.
# Split out as its own third top-level index.
PERFORMANCE_DOMAIN <- "Performance"

# =============================================================================
# STEP 1 - LOAD
# =============================================================================
log_msg("[1/7] Loading biomarker sheet '%s'", BIOMARKER_SHEET)
bio <- as.data.frame(read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
log_msg("  %d rows x %d columns", nrow(bio), ncol(bio))
for (req in c("Study_ID", "Sample_ID")) {
  if (!req %in% names(bio)) stop(sprintf("Required column '%s' not found", req))
}

# =============================================================================
# STEP 2 - WITHIN-STUDY Z-SCORING
#
# Absolute biomarker values are not comparable across studies: different assay
# kits, units, reference ranges and bird ages mean a value that is high in one
# study may be ordinary in another. Standardising within study expresses each
# sample relative to its own study's samples, which is what makes values
# poolable at all, and removes the need for cross-study unit conversion.
#
# A marker cannot be standardised within a study that has fewer than two
# samples carrying it, or where all its values are identical (zero variance).
# Those yield NA rather than a fabricated number, and the loss is reported.
# =============================================================================
log_msg("[2/7] Within-study z-scoring")

zscore_within_study <- function(values, study) {
  out <- rep(NA_real_, length(values))
  for (s in unique(study)) {
    sel <- study == s & !is.na(values)
    v <- values[sel]
    if (length(v) < 2) next
    sdv <- stats::sd(v)
    if (!is.finite(sdv) || sdv == 0) next
    out[sel] <- (v - mean(v)) / sdv
  }
  out
}

all_markers <- unlist(lapply(DOMAINS, names), use.names = FALSE)
missing_from_sheet <- setdiff(all_markers, names(bio))
if (length(missing_from_sheet) > 0) {
  log_msg("  NOTE: %d configured marker(s) absent from the sheet, skipped: %s",
          length(missing_from_sheet), paste(missing_from_sheet, collapse = ", "))
}
usable <- intersect(all_markers, names(bio))

for (m in usable) {
  raw <- suppressWarnings(as.numeric(bio[[m]]))
  bio[[paste0("Z_", m)]] <- zscore_within_study(raw, bio$Study_ID)
  n_raw <- sum(!is.na(raw)); n_z <- sum(!is.na(bio[[paste0("Z_", m)]]))
  if (n_raw > n_z) {
    log_msg("  %-6s %d raw -> %d z-scored (%d in studies with <2 values or zero variance)",
            m, n_raw, n_z, n_raw - n_z)
  }
}

# =============================================================================
# STEP 3 - ORIENT SO HIGHER ALWAYS MEANS HEALTHIER
# =============================================================================
log_msg("[3/7] Orienting markers by direction of effect")
n_flipped <- 0
for (dom in names(DOMAINS)) {
  for (m in names(DOMAINS[[dom]])) {
    if (!m %in% usable) next
    d <- DOMAINS[[dom]][[m]]
    bio[[paste0("O_", m)]] <- bio[[paste0("Z_", m)]] * d
    if (d < 0) n_flipped <- n_flipped + 1
  }
}
log_msg("  %d marker(s) inverted (higher raw value = worse gut health)", n_flipped)

# =============================================================================
# STEP 4 - DOMAIN SCORES BY PCA (protocol method), FALLING BACK TO MEAN
#
# PC1 is taken as the domain score. Its sign is arbitrary in PCA output, so it
# is fixed here to correlate positively with the mean of the domain's oriented
# markers - without this, a domain could silently invert between runs.
# =============================================================================
log_msg("[4/7] Domain scores (PCA where complete cases >= %d, else mean-of-available)",
        MIN_N_FOR_PCA)

domain_info <- list()
domain_cols <- character(0)

for (dom in names(DOMAINS)) {
  cols <- paste0("O_", intersect(names(DOMAINS[[dom]]), usable))
  if (length(cols) == 0) {
    log_msg("  %-16s SKIPPED - no markers available", dom); next
  }
  M <- as.matrix(bio[, cols, drop = FALSE])
  complete <- stats::complete.cases(M)
  n_complete <- sum(complete)
  scol <- paste0("GHI_dom_", dom)
  bio[[scol]] <- NA_real_
  bio[[paste0("GHI_nmark_", dom)]] <- rowSums(!is.na(M))

  if (n_complete >= MIN_N_FOR_PCA && length(cols) >= 2) {
    p <- stats::prcomp(M[complete, , drop = FALSE], center = TRUE, scale. = FALSE)
    pc1 <- p$x[, 1]
    varexp <- (p$sdev^2 / sum(p$sdev^2))[1]
    # fix the arbitrary PC1 sign against the mean of the oriented markers
    ref <- rowMeans(M[complete, , drop = FALSE])
    if (stats::cor(pc1, ref) < 0) { pc1 <- -pc1; p$rotation[, 1] <- -p$rotation[, 1] }

    # -------------------------------------------------------------------------
    # IS PC1 MEASURING THE INTENDED CONSTRUCT?
    #
    # PCA returns the direction of greatest variance, which is not necessarily
    # the direction that means gut health. Where a domain mixes orientation
    # signs, the two can come apart entirely.
    #
    # Morphology is the case in point, and the only domain affected: villus
    # height carries +1 and crypt depth -1, while every other domain has all
    # its markers pointing the same way. Raw villus height and crypt depth
    # correlate POSITIVELY (+0.44 within study) because both scale with bird
    # size. After orienting crypt depth by -1 the pair correlates NEGATIVELY,
    # and PCA on a negatively correlated pair extracts their contrast - which
    # in oriented terms is O_VH - O_CD, and in raw terms is z(VH) + z(CD).
    #
    # That is the SUM of the two markers: the size axis. A bird with taller
    # villi and deeper crypts scores higher, which is growth, not health. The
    # intended quantity is z(VH) - z(CD), the VH:CD contrast.
    #
    # Measured against the published scores, the domain correlated +1.000 with
    # the size axis and +0.017 with the health contrast.
    #
    # The sign-fixing line above cannot catch this. When the oriented pair is
    # negatively correlated, PC1 and the mean of the oriented markers are
    # mathematically orthogonal, so the correlation it tests sits at zero and
    # the flip never fires.
    #
    # The guard below tests the same correlation for MAGNITUDE rather than
    # sign. If PC1 is near-orthogonal to the mean of the oriented markers, PC1
    # is not measuring the construct and the mean is used instead - which for
    # Morphology is exactly the VH:CD contrast the field reports.
    # -------------------------------------------------------------------------
    align <- suppressWarnings(stats::cor(pc1, ref))
    if (!is.finite(align) || abs(align) < PC1_ALIGN_MIN) {
      log_msg("  %-16s PC1 rejected: alignment with the oriented-marker mean is %+.3f",
              dom, ifelse(is.finite(align), align, NA))
      log_msg("  %-16s   PC1 is capturing variance orthogonal to the intended", "")
      log_msg("  %-16s   construct; using the mean of oriented markers instead.", "")
      pc1 <- ref
      p$rotation[, 1] <- rep(1 / length(cols), length(cols))
      varexp <- NA_real_
      method_used <- "MEAN (PC1 rejected)"
    } else method_used <- "PCA"

    bio[[scol]][complete] <- pc1

    # -------------------------------------------------------------------------
    # CALIBRATED FILL FOR INCOMPLETE ROWS
    #
    # PC1 is only defined for rows carrying every marker of the domain. Leaving
    # the rest as NA discards usable data, and unevenly: Inflammation lost 63
    # of 166 rows that carried a marker, while Barrier lost 5 and SCFA lost 1.
    # That asymmetry is an artefact of which domains happened to clear the
    # complete-case threshold for PCA, not a biological difference - a study
    # measuring one cytokine is no less informative about inflammation than a
    # study measuring one immunoglobulin is about mucosal immunity.
    #
    # Incomplete rows are therefore scored by the mean of whichever oriented
    # markers they carry, which is exactly how sparse domains are scored
    # already. The two cannot simply be concatenated, because a mean of
    # oriented z-scores and a PC1 score have different centres and spreads;
    # mixing them raw would shift a sample up or down for a reason that is
    # arithmetic rather than biological.
    #
    # So the fill is CALIBRATED. On the complete rows, where both quantities
    # exist, the linear map from mean-of-markers onto PC1 is estimated by
    # matching centre and spread. That map is then applied to the incomplete
    # rows. The result sits on the PC1 scale by construction, and on complete
    # rows nothing changes - PC1 is used as before.
    #
    # The two constants are exported in the index definition so that
    # score_new_ghi.R reproduces the same behaviour on new data.
    # -------------------------------------------------------------------------
    fill_off <- NA_real_; fill_sc <- NA_real_; n_filled <- 0L
    mean_all <- ifelse(rowSums(!is.na(M)) >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
    ref_c <- mean_all[complete]
    if (stats::sd(ref_c, na.rm = TRUE) > 0 && is.finite(stats::sd(ref_c, na.rm = TRUE))) {
      fill_sc  <- stats::sd(pc1) / stats::sd(ref_c)
      fill_off <- mean(pc1) - fill_sc * mean(ref_c)
      need <- !complete & !is.na(mean_all)
      if (any(need)) {
        bio[[scol]][need] <- fill_off + fill_sc * mean_all[need]
        n_filled <- sum(need)
      }
    }
    if (n_filled > 0)
      log_msg("  %-16s   + %d incomplete row(s) filled by calibrated mean (scale %.3f)",
              "", n_filled, fill_sc)

    domain_info[[dom]] <- list(method = method_used, n = n_complete, varexp = varexp,
                               loadings = p$rotation[, 1], markers = cols,
                               fill_offset = fill_off, fill_scale = fill_sc,
                               n_filled = n_filled)
    # varexp is NA when PC1 was rejected by the alignment guard, since the score
    # is then the mean of the oriented markers and no component was retained.
    if (is.na(varexp)) {
      log_msg("  %-16s MEAN  n=%-4d equal weight on each oriented marker",
              dom, n_complete)
      log_msg("        markers: %s",
              paste(sub("^O_", "", cols), collapse = "  "))
    } else {
      flag <- if (varexp < PC1_VARIANCE_FLOOR) "  <-- BELOW variance floor, inspect loadings" else ""
      log_msg("  %-16s PCA   n=%-4d PC1 explains %.1f%%%s", dom, n_complete, 100*varexp, flag)
      log_msg("        loadings: %s",
              paste(sprintf("%s=%+.2f", sub("^O_", "", names(p$rotation[,1])), p$rotation[,1]),
                    collapse = "  "))
    }
  } else {
    score <- ifelse(rowSums(!is.na(M)) >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
    score[!is.finite(score)] <- NA_real_
    bio[[scol]] <- score
    domain_info[[dom]] <- list(method = "mean-of-available", n = sum(!is.na(score)),
                               varexp = NA_real_, loadings = NULL, markers = cols)
    log_msg("  %-16s MEAN  n=%-4d (only %d complete cases - too few for stable PCA)",
            dom, sum(!is.na(score)), n_complete)
  }
  domain_cols <- c(domain_cols, scol)
}

# =============================================================================
# STEP 5 - COMBINE DOMAINS INTO THE THREE INDEX VERSIONS
#
# Cross-domain PCA is not computable here (zero samples carry every marker),
# so domain scores are combined by mean of those available. Arithmetic rather
# than geometric mean, because these scores are z-score/PC1 based and centred
# on zero - a geometric mean requires strictly positive input.
# =============================================================================
log_msg("[5/7] Combining domains into the three index versions")

combine <- function(dt, cols) {
  if (length(cols) == 0) return(list(score = rep(NA_real_, nrow(dt)), n = rep(0L, nrow(dt))))
  M <- as.matrix(dt[, cols, drop = FALSE])
  n <- rowSums(!is.na(M))
  s <- ifelse(n >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
  s[!is.finite(s)] <- NA_real_
  list(score = s, n = n)
}

morph_col <- paste0("GHI_dom_", MORPHOLOGY_DOMAIN)

if (AGGREGATION == "flat") {
  # ---------------------------------------------------------------------------
  # FLAT AGGREGATION
  #
  # Every oriented marker averaged directly, with no domain layer. The domain
  # structure is what gives each facet of gut health equal weight regardless of
  # how many markers a study happened to measure for it; removing it means a
  # sample carrying four barrier markers and one morphology marker is scored
  # mostly on barrier. Reported as an alternative construction, not as the
  # registered one.
  # ---------------------------------------------------------------------------
  o_all   <- grep("^O_", names(bio), value = TRUE)
  o_morph <- intersect(paste0("O_", names(DOMAINS[[MORPHOLOGY_DOMAIN]])), o_all)
  o_perf  <- intersect(paste0("O_", names(DOMAINS[[PERFORMANCE_DOMAIN]])), o_all)
  o_bio   <- setdiff(o_all, c(o_morph, o_perf))
  flat <- function(cols) {
    if (!length(cols)) return(list(score = rep(NA_real_, nrow(bio)), n = 0))
    M <- as.matrix(bio[, cols, drop = FALSE])
    n <- rowSums(!is.na(M))
    s <- ifelse(n >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
    s[!is.finite(s)] <- NA_real_
    list(score = s, n = n)
  }
  log_msg("  AGGREGATION = flat: %d oriented markers averaged directly, no domain layer",
          length(o_all))
  r_morph <- flat(o_morph)
  r_bio   <- flat(o_bio)
  r_perf  <- flat(o_perf)
  r_comb  <- flat(o_all)
} else {
  perf_col <- paste0("GHI_dom_", PERFORMANCE_DOMAIN)
  r_morph <- combine(bio, intersect(morph_col, domain_cols))
  r_perf  <- combine(bio, intersect(perf_col, domain_cols))
  r_bio   <- combine(bio, setdiff(domain_cols, c(morph_col, perf_col)))
  r_comb  <- combine(bio, domain_cols)
}

bio$GHI_Morphology_raw <- r_morph$score
bio$GHI_Biomarker_raw  <- r_bio$score
bio$GHI_Performance_raw <- r_perf$score
bio$GHI_Combined_raw   <- r_comb$score
bio$GHI_n_domains      <- r_comb$n

for (v in c("Morphology","Biomarker","Performance","Combined")) {
  log_msg("  GHI_%-12s %4d samples", v, sum(!is.na(bio[[paste0("GHI_", v, "_raw")]])))
}

# =============================================================================
# STEP 6 - RESCALE TO 0-100
#
# The raw scores are z-score based and centred on zero, which is fine
# statistically but not interpretable to a non-statistician. Rescaling is a
# linear transform of the observed range and changes nothing about the
# ordering or the statistics - it only makes the number readable.
#
# CONVENTION: 100 = best observed gut health, 0 = worst observed. Note this is
# the inverse of the dysbiosis-style scale used for the GDI, where a high
# number denotes more dysbiosis; keep the two straight when reporting.
#
# The min and max are corpus-specific, so a score is only comparable to other
# scores from this same corpus. That is a real limitation for external
# application and is stated in the output summary.
# =============================================================================
log_msg("[6/7] Rescaling to 0-100 (100 = best observed)")

rescale_100 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(NA_real_, length(x)))
  100 * (x - rng[1]) / diff(rng)
}
for (v in c("Morphology","Biomarker","Performance","Combined")) {
  bio[[paste0("GHI_", v)]] <- rescale_100(bio[[paste0("GHI_", v, "_raw")]])
}

# =============================================================================
# STEP 7 - ROC ANALYSIS AGAINST KNOWN CONTROL / CHALLENGE STATUS
#
# The corpus carries an independent label - whether a sample came from a
# control group or a challenged group - which was never used in building the
# index. That makes it a genuine external check: can the index separate birds
# already known to be healthy from birds already known to be challenged?
#
# AUC is the probability that a randomly chosen control scores healthier than
# a randomly chosen challenged bird. 0.5 is coin-flip; higher is better. The
# reported cutoff is the point maximising sensitivity + specificity (Youden's
# J), which is what turns the index from a number into a decision rule.
#
# Implemented directly rather than via pROC so the script has no extra
# dependency; AUC is computed by the Mann-Whitney identity, which is exact.
# =============================================================================
log_msg("[7/7] ROC analysis vs control/challenge status")

roc_eval <- function(score, is_control) {
  ok <- !is.na(score) & !is.na(is_control)
  s <- score[ok]; g <- is_control[ok]
  n_pos <- sum(g); n_neg <- sum(!g)
  if (n_pos < 5 || n_neg < 5) return(NULL)
  # AUC via rank-sum (equivalent to the trapezoidal ROC area, no ties bias)
  r <- rank(s)
  auc <- (sum(r[g]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  # Youden's J over every observed cut point
  cuts <- sort(unique(s))
  best <- list(J = -Inf)
  for (cp in cuts) {
    pred_ctrl <- s >= cp
    sens <- sum(pred_ctrl & g) / n_pos
    spec <- sum(!pred_ctrl & !g) / n_neg
    J <- sens + spec - 1
    if (J > best$J) best <- list(J = J, cut = cp, sens = sens, spec = spec)
  }
  list(auc = auc, n_control = n_pos, n_challenge = n_neg,
       cut = best$cut, sens = best$sens, spec = best$spec)
}

# Defining the two groups correctly matters more than the ROC maths. The
# corpus labels several distinct conditions, and only some are challenges:
#
#   Control                  - untreated baseline            -> HEALTHY
#   Non_Disruptive_Additive  - probiotic / phytogenic / etc.  -> not a challenge
#   Pathogen_Challenge       - infectious challenge           -> CHALLENGED
#   Antibiotic_Challenge     - antibiotic-induced disruption  -> CHALLENGED
#   Dietary_Feed_Challenge   - dietary insult                 -> CHALLENGED
#
# Beneficial additives are an intervention, not a health state: a bird given a
# probiotic is not "challenged", and should score at least as healthy as an
# untreated control. Counting those rows as challenged inverts the comparison
# and makes a working index look worse than chance. They are therefore
# EXCLUDED from the ROC contrast rather than assigned to either arm - the test
# asks only whether the index separates untreated controls from genuinely
# challenged birds.
# Old vocabulary kept for reference; the ROC section below now scores on
# Intervention_Class, so Bernard_2024's low-energy arms are no longer counted
# as challenged and Ding_2025's depletion cocktail now is.
CHALLENGE_LEVELS_OLD <- c("Pathogen_Challenge", "Antibiotic_Challenge",
                          "Dietary_Feed_Challenge")
CHALLENGE_LEVELS <- c("Dietary challenge", "Infectious challenge",
                      "Toxin challenge", "Antibiotic challenge",
                      "Environmental challenge")
HEALTHY_LEVELS   <- c("Unchallenged control")
CLASS_COL        <- "Intervention_Class"

# Studies whose Control rows are not contemporaneous with their challenged
# rows, and so cannot form a within-study pair. Massacci_2019 challenged every
# bird on d21; its Control rows are pre-challenge d21 faecal rows and its
# challenged rows are d40, so the pair differs by 19 days of growth and shares
# no domain. Excluded from the CONCORDANCE count only -- its d40 rows are
# genuinely challenged birds and remain in the pooled ROC.
NO_CONCURRENT_CONTROL <- c("Massacci_2019")

# -----------------------------------------------------------------------------
# WITHIN-STUDY CONCORDANCE - the primary validation for this corpus
#
# Pooled ROC is not a fair test here, and the reason is structural. Markers are
# z-scored WITHIN study, which centres every study on zero. In a study
# comparing a control against a beneficial additive, the control therefore
# scores BELOW its study mean; in a study comparing a control against a
# pathogen challenge, the control scores ABOVE it. Pooling the two cancels the
# signal, and this corpus contains more additive arms than challenge arms, so
# the cancellation is severe. A pooled AUC near 0.5 under these conditions
# reflects the standardisation scheme, not a failure of the index.
#
# The question the data can actually answer is paired and within-study: inside
# a single study, does the untreated control score healthier than that same
# study's challenged arms? Every term is on the same z-scale, so the
# comparison is valid. Reported as the proportion of eligible studies in which
# the control mean exceeds the challenged mean, which is a sign test - its
# null is 50%, the same null as an AUC of 0.5, so the two are directly
# comparable in interpretation.
# -----------------------------------------------------------------------------
within_study_concordance <- function(score, study, is_ctrl, is_chal) {
  ok <- !is.na(score)
  studies <- unique(study[ok])
  n_eval <- 0; n_correct <- 0; detail <- character(0)
  for (s in studies) {
    sel <- ok & study == s
    c_vals <- score[sel & is_ctrl]
    h_vals <- score[sel & is_chal]
    if (length(c_vals) == 0 || length(h_vals) == 0) next
    n_eval <- n_eval + 1
    correct <- mean(c_vals) > mean(h_vals)
    if (correct) n_correct <- n_correct + 1
    detail <- c(detail, sprintf("    %-24s ctrl=%6.1f  chal=%6.1f  %s",
                s, mean(c_vals), mean(h_vals), if (correct) "OK" else "INVERTED"))
  }
  list(n_eval = n_eval, n_correct = n_correct,
       prop = if (n_eval > 0) n_correct / n_eval else NA_real_,
       detail = detail)
}

conc_results <- list()

roc_results <- list()
if ("Challenge_Type" %in% names(bio)) {
  in_contrast <- bio[[CLASS_COL]] %in% c(HEALTHY_LEVELS, CHALLENGE_LEVELS)
  is_control  <- ifelse(in_contrast, bio[[CLASS_COL]] %in% HEALTHY_LEVELS, NA)
  log_msg("  contrast: %d control vs %d challenged (%d additive/other rows excluded)",
          sum(is_control, na.rm = TRUE), sum(!is_control, na.rm = TRUE),
          sum(!in_contrast))
  pairable <- !bio$Study_ID %in% NO_CONCURRENT_CONTROL
  is_ctrl_v <- bio[[CLASS_COL]] %in% HEALTHY_LEVELS   & pairable
  is_chal_v <- bio[[CLASS_COL]] %in% CHALLENGE_LEVELS & pairable
  for (v in c("Morphology","Biomarker","Performance","Combined")) {
    cc <- within_study_concordance(bio[[paste0("GHI_", v)]], bio$Study_ID,
                                   is_ctrl_v, is_chal_v)
    conc_results[[v]] <- cc
    if (cc$n_eval > 0) {
      log_msg("  GHI_%-12s within-study concordance: %d/%d studies (%.0f%%) score control > challenged",
              v, cc$n_correct, cc$n_eval, 100 * cc$prop)
    } else {
      log_msg("  GHI_%-12s no study has both a control and a challenged arm", v)
    }
    r <- roc_eval(bio[[paste0("GHI_", v)]], is_control)
    roc_results[[v]] <- r
    if (is.null(r)) {
      log_msg("  GHI_%-12s too few control or challenge samples for ROC", v)
    } else {
      log_msg("  GHI_%-12s AUC=%.3f  (n_control=%d, n_challenge=%d)  cutoff=%.1f  sens=%.2f  spec=%.2f",
              v, r$auc, r$n_control, r$n_challenge, r$cut, r$sens, r$spec)
    }
  }
} else {
  log_msg("  Challenge_Type column not found - ROC skipped")
}

# =============================================================================
# EXPORT
# =============================================================================
keep <- c(intersect(c("Study_ID","Sample_ID","Treatment_Group","Segment","Age_d",
                      "Challenge_Type", "Intervention_Class",
                      "Intervention_Subclass",
                      "Challenge_Contrast_Eligible"), names(bio)),
          domain_cols, grep("^GHI_nmark_", names(bio), value = TRUE),
          "GHI_Morphology_raw","GHI_Biomarker_raw","GHI_Performance_raw","GHI_Combined_raw",
          "GHI_Morphology","GHI_Biomarker","GHI_Performance","GHI_Combined","GHI_n_domains")
out <- file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")
fwrite(bio[, keep, drop = FALSE], out, sep = "\t")
log_msg("per-sample scores -> %s", out)

# -----------------------------------------------------------------------------
# PORTABLE INDEX DEFINITION
#
# Everything needed to compute the GHI on a study that was not in this corpus:
# which markers belong to which domain, the direction each is oriented in, and
# the PC1 loading used to combine them where a domain was scored by PCA.
#
# Domains scored by mean-of-available carry no loadings; those are recomputed
# as a simple mean of the oriented markers present. Both cases are marked in
# the method column.
#
# fill_offset and fill_scale apply to PCA domains only. A row missing one or
# more markers of the domain cannot be given a PC1 score, and is instead scored
# as fill_offset + fill_scale x (mean of its oriented markers). The two
# constants place that value on the PC1 scale, and were estimated on the rows
# where both quantities exist.
#
# Note the index is study-relative by construction: markers are z-scored within
# study before orientation, so a new study is scored against its own baseline
# rather than against this corpus.
# -----------------------------------------------------------------------------
defn <- rbindlist(lapply(names(domain_info), function(dom) {
  di <- domain_info[[dom]]
  mk <- sub("^O_", "", di$markers)
  data.table(
    domain    = dom,
    marker    = mk,
    direction = vapply(mk, function(m) DOMAINS[[dom]][[m]], numeric(1)),
    method    = di$method,
    pc1_loading = if (is.null(di$loadings)) NA_real_ else as.numeric(di$loadings),
    pc1_varexp  = if (is.na(di$varexp)) NA_real_ else round(di$varexp, 4),
    n_complete  = di$n,
    # constants for scoring rows that lack one or more markers of the domain
    fill_offset = if (is.null(di$fill_offset)) NA_real_ else di$fill_offset,
    fill_scale  = if (is.null(di$fill_scale))  NA_real_ else di$fill_scale)
}), fill = TRUE)
fdef <- file.path(R_OUT_DIR, "ghi", "ghi_index_definition.tsv")
fwrite(defn, fdef, sep = "\t")
log_msg("INDEX DEFINITION -> %s", fdef)
log_msg("  %d markers across %d domains - needed to score a new study",
        nrow(defn), uniqueN(defn$domain))

lines <- c(
  "=== GHI construction by PCA (registered protocol method) ===", "",
  sprintf("Input rows: %d", nrow(bio)),
  sprintf("PCA used where complete cases >= %d; mean-of-available otherwise.", MIN_N_FOR_PCA),
  "",
  "--- Domain scores ---")
for (dom in names(domain_info)) {
  di <- domain_info[[dom]]
  lines <- c(lines, sprintf("  %-16s %-18s n=%-5d %s", dom, di$method, di$n,
    if (is.na(di$varexp)) "" else sprintf("PC1 explains %.1f%% of variance", 100*di$varexp)))
  if (!is.null(di$loadings)) {
    lines <- c(lines, sprintf("      loadings: %s",
      paste(sprintf("%s=%+.2f", sub("^O_", "", names(di$loadings)), di$loadings), collapse = "  ")))
  }
}
lines <- c(lines, "", "--- Index coverage (0-100 scale, 100 = best observed) ---")
for (v in c("Morphology","Biomarker","Performance","Combined")) {
  x <- bio[[paste0("GHI_", v)]]; x <- x[!is.na(x)]
  lines <- c(lines, sprintf("  GHI_%-12s n=%-5d median=%.1f  IQR %.1f-%.1f", v, length(x),
                            stats::median(x), stats::quantile(x, .25), stats::quantile(x, .75)))
}
lines <- c(lines, "", "--- Within-study concordance (PRIMARY validation) ---",
  "Proportion of studies in which the untreated control scores healthier than",
  "that same study's challenged arms. Null expectation 50%.")
for (v in names(conc_results)) {
  cc <- conc_results[[v]]
  lines <- c(lines, if (cc$n_eval == 0) sprintf("  GHI_%-12s not evaluable", v) else
    c(sprintf("  GHI_%-12s %d/%d studies (%.0f%%)", v, cc$n_correct, cc$n_eval, 100*cc$prop),
      cc$detail))
}
lines <- c(lines, "", "--- Pooled ROC (SECONDARY - see caveat) ---",
  "Interpret with care: within-study z-scoring centres each study on zero, so",
  "a control in an additive study scores below its study mean while a control",
  "in a challenge study scores above it. Pooling cancels much of the signal,",
  "and this corpus has more additive than challenge arms. An AUC near 0.5 here",
  "reflects that standardisation, not necessarily a weak index - the",
  "within-study concordance above is the meaningful test.")
for (v in names(roc_results)) {
  r <- roc_results[[v]]
  lines <- c(lines, if (is.null(r)) sprintf("  GHI_%-12s not evaluable", v) else
    sprintf("  GHI_%-12s AUC=%.3f  cutoff=%.1f  sens=%.2f  spec=%.2f  (n=%d control / %d challenge)",
            v, r$auc, r$cut, r$sens, r$spec, r$n_control, r$n_challenge))
}
lines <- c(lines, "",
  "--- PRE-SPECIFIED PRIMARY TARGET ---",
  "GHI_Combined is the pre-specified primary outcome for the downstream GDI",
  "regression. GHI_Morphology, GHI_Biomarker and GHI_Performance are pre-",
  "specified sensitivity analyses.",
  "",
  "Fixed BEFORE any 16S sequence data was used, on grounds independent of any",
  "validation result:",
  "  (a) closest executable form of the registered protocol, which specifies a",
  "      composite derived from normalised biomarker data - all biomarkers,",
  "      not a subset;",
  "  (b) largest sample of the three versions;",
  "  (c) uses all available biology rather than discarding domains.",
  "",
  "It was NOT chosen because it performed best. Under within-study concordance",
  "GHI_Morphology ranked highest; under random-effects meta-analysis of paired",
  "study differences GHI_Biomarker ranked highest. Selecting whichever version",
  "won a given test would be selection on the outcome. Choosing the most",
  "protocol-faithful version a priori, and reporting the other two as declared",
  "sensitivity analyses, avoids that.",
  "",
  "--- Deviations from the registered protocol, and why ---",
  "PCA is applied within each domain, as registered. It is NOT applied across",
  "domains: zero samples in this corpus carry every biomarker, because which",
  "markers a study measured is determined by what that study was designed to",
  "test. Complete-case PCA across domains is therefore undefined, not merely",
  "underpowered. Domain scores are combined by arithmetic mean instead.",
  "",
  "VH:CD is excluded from the Morphology PCA because it is the arithmetic ratio",
  "of VH to CD; including it alongside its own numerator and denominator would",
  "distort the loadings. It is retained in the data as a descriptive measure.",
  "",
  "0-100 rescaling is anchored to the observed range in THIS corpus, so scores",
  "are comparable within this dataset but not directly transferable to an",
  "external one without re-anchoring.",
  "",
  "The immunoglobulin panel is scored as TWO domains, not one:",
  "  Immunoglobulin_Mucosal   (sIgA, IgA)  - competence, +1 defensible",
  "  Immunoglobulin_Systemic  (IgG, IgM)   - direction unresolved",
  "",
  "The split follows the domain correlation diagnostic (04g), where the",
  "combined domain ran negative against Morphology (r = -0.247, n = 68) and",
  "flat against Barrier and Claudin1 - the signature of markers pulling in",
  "opposite directions inside one domain. Same treatment, same rationale, as",
  "claudin-1 being separated from ZO-1 and occludin.",
  "",
  "The +1 on IgG / IgM orients the systemic domain for PCA. It is NOT a claim",
  "that higher systemic immunoglobulin is healthier. Downstream taxon selection",
  "should be run both ways as a sensitivity analysis.")
outs <- file.path(R_OUT_DIR, "ghi", "ghi_pca_summary.txt")
writeLines(lines, outs)
log_msg("summary -> %s", outs)
log_msg("=== DONE ===")
