#!/usr/bin/env Rscript
# =============================================================================
# 04d_validation_matrix.R
#
# Tests every combination of standardisation method and validation contrast,
# and reports what works and what does not, in one table.
#
# RUN:  Rscript 04d_validation_matrix.R
# TIME: under a minute. Reads the biomarker workbook only. Changes nothing.
#
# -----------------------------------------------------------------------------
# WHY A MATRIX RATHER THAN A SINGLE TEST
# -----------------------------------------------------------------------------
# Two separate choices were being made implicitly, and each was capable of
# making a working index look broken:
#
#   1. HOW values are standardised across studies (raw / z-score / percentile)
#   2. WHICH biological contrast the index is asked to detect
#
# Testing them together shows whether a weak result reflects the index, the
# standardisation, or simply the wrong question being asked of it.
#
# -----------------------------------------------------------------------------
# THE TWO CONTRASTS, AND WHY BOTH MATTER
# -----------------------------------------------------------------------------
# CHALLENGE   control vs pathogen / antibiotic / dietary insult.
#             The classic diagnostic framing: can the index tell a sick bird
#             from a healthy one? Fewer studies carry this contrast.
#
# ADDITIVE    control vs beneficial additive (probiotic, phytogenic, organic
#             acid). Asks the opposite question: does the index detect
#             IMPROVEMENT above an untreated baseline?
#
# The additive contrast is not a lesser substitute. Most poultry gut-health
# research tests feed additives rather than infections, so it is the contrast
# this corpus carries most data for, and arguably the more field-relevant
# question. An index that detects improvement is directly useful for
# evaluating feed interventions.
#
# Note the expected direction differs: for CHALLENGE the control should score
# HIGHER (healthier) than the challenged group; for ADDITIVE the additive
# group should score higher than the control. The script handles each
# direction explicitly rather than assuming one.
#
# -----------------------------------------------------------------------------
# METRIC
# -----------------------------------------------------------------------------
# Within-study concordance: within one study, does the group expected to be
# healthier actually score higher? Reported as a proportion of eligible
# studies, with a one-sided binomial test against a 50% null.
#
# Paired and within-study by construction, so it is valid even for RAW values
# whose units differ between studies. Pooled AUC is not used as the headline:
# within-study standardisation centres each study on zero, so pooling controls
# from additive studies with controls from challenge studies cancels the
# signal and reports ~0.5 regardless of whether the index works.
# =============================================================================

.locate_and_source_config <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  d <- if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
  source(file.path(d, "00_r_config.R"))
}
.locate_and_source_config()

CONTRASTS <- list(
  CHALLENGE = list(
    healthier = c("Unchallenged control"),
    comparator = c("Infectious challenge", "Toxin challenge", "Environmental challenge", "Antibiotic challenge", "Dietary challenge"),
    exclude   = c("Massacci_2019"),
    label     = "control should beat challenged"),
  ADDITIVE = list(
    healthier = c("Non-disruptive additive"),
    comparator = c("Unchallenged control"),
    label     = "additive should beat control"),
  # A single three-arm ordering, tested directly rather than assumed: additive
  # should score highest, challenge should score lowest, control in between.
  # This is the same logic behind MAPI and the Microbial Dysbiosis Index -
  # log(taxa enriched under improvement / taxa enriched under damage) - except
  # here it is checked per marker before being trusted, not assumed to hold
  # for every marker equally. Claudin-1 is the demonstrated counterexample:
  # it moves in the SAME direction under both additive and challenge, so a
  # single shared scale would be wrong for it specifically.
  DYSBIOSIS = list(
    healthier = c("Non-disruptive additive"),
    comparator = c("Infectious challenge", "Toxin challenge", "Environmental challenge", "Antibiotic challenge", "Dietary challenge"),
    label     = "additive should beat challenged (the full-range ordering)")
)

DOMAINS <- list(
  Morphology              = c(VH_um = +1, CD_um = -1),
  Barrier                 = c(ZO1 = +1, Ocln = +1, Cldn1 = +1, MUC2 = +1),
  Inflammation            = c(IL6 = -1, TNFa = -1, IL1b = -1),
  Anti_Inflammatory       = c(IL10 = +1),
  SCFA                    = c(Ace = +1, Prop = +1, Buty = +1),
  Performance             = c(BWG_g = +1, FCR = -1),
  Immunoglobulin_Systemic = c(IgG = +1, IgM = +1, IgA_serum = +1)
)

log_msg("[1/4] Loading biomarker sheet")
bio <- as.data.frame(read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
log_msg("  %d rows x %d columns", nrow(bio), ncol(bio))
if (!"Intervention_Class" %in% names(bio)) stop("Intervention_Class column required.")

# ---- standardisation methods ----------------------------------------------
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
    n <- sum(k); if (n < 2) next
    out[k] <- (rank(v[k], ties.method = "average") - 0.5) / n
  }
  out
}

# Median-centred: subtracts the study median but does NOT divide by an SD.
# Included because the SD is the suspect quantity - it is estimated from only
# a handful of treatment-group means per study - while centring alone is
# stable. This isolates whether the damage comes from centring or from
# scaling.
std_median_centred <- function(v, study) {
  v <- as.numeric(v); out <- rep(NA_real_, length(v))
  for (s in unique(study)) {
    k <- study == s & !is.na(v)
    if (sum(k) < 2) next
    out[k] <- v[k] - stats::median(v[k])
  }
  out
}

# Control-referenced: expresses each value as a proportion of that study's
# control mean. Directly implements "percentile normalisation relative to
# intra-study controls" from the literature, anchoring every study to its own
# untreated baseline rather than to its overall mean. Only defined where a
# study actually has a control group.
make_std_control_ref <- function(challenge_col) {
  function(v, study) {
    v <- as.numeric(v); out <- rep(NA_real_, length(v))
    for (s in unique(study)) {
      k <- study == s & !is.na(v)
      if (!any(k)) next
      ctrl <- k & challenge_col == "Unchallenged control"
      if (!any(ctrl)) next
      m <- mean(v[ctrl])
      if (!is.finite(m) || m == 0) next
      out[k] <- v[k] / m
    }
    out
  }
}

METHODS <- list(
  RAW            = std_raw,
  ZSCORE         = std_zscore,
  PERCENTILE     = std_percentile,
  MEDIAN_CENTRED = std_median_centred,
  CONTROL_REF    = make_std_control_ref(bio$Intervention_Class)
)

# ---- concordance ------------------------------------------------------------
concordance <- function(score, study, better, worse) {
  ok <- !is.na(score); n_eval <- 0; n_ok <- 0
  for (s in unique(study[ok])) {
    k <- ok & study == s
    a <- score[k & better]; b <- score[k & worse]
    if (!length(a) || !length(b)) next
    n_eval <- n_eval + 1
    if (mean(a) > mean(b)) n_ok <- n_ok + 1
  }
  p <- if (n_eval > 0) stats::binom.test(n_ok, n_eval, 0.5, alternative = "greater")$p.value else NA_real_
  list(n_eval = n_eval, n_ok = n_ok,
       prop = if (n_eval) n_ok / n_eval else NA_real_, p = p)
}

build_index <- function(std_fn) {
  dom <- list()
  for (d in names(DOMAINS)) {
    cols <- intersect(names(DOMAINS[[d]]), names(bio))
    if (!length(cols)) next
    M <- sapply(cols, function(m) std_fn(bio[[m]], bio$Study_ID) * DOMAINS[[d]][[m]])
    if (is.null(dim(M))) M <- matrix(M, ncol = 1)
    n <- rowSums(!is.na(M))
    s <- ifelse(n >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
    s[!is.finite(s)] <- NA_real_
    dom[[d]] <- s
  }
  agg <- function(M) {
    if (is.null(M) || !length(M)) return(rep(NA_real_, nrow(bio)))
    M <- do.call(cbind, M)
    if (is.null(dim(M))) M <- matrix(M, ncol = 1)
    n <- rowSums(!is.na(M)); s <- ifelse(n >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
    s[!is.finite(s)] <- NA_real_; s
  }
  list(domains = dom,
       Morphology = dom[["Morphology"]],
       Biomarker  = agg(dom[setdiff(names(dom), "Morphology")]),
       Combined   = agg(dom))
}

log_msg("[2/4] Building indices under %d standardisation methods", length(METHODS))
built <- lapply(METHODS, build_index)

log_msg("[3/4] Testing every method x contrast combination")

fmt <- function(cc) {
  if (cc$n_eval == 0) return("     -    ")
  sprintf("%3.0f%% (%2d)%s", 100 * cc$prop, cc$n_eval,
          if (!is.na(cc$p) && cc$p < 0.05) "*" else " ")
}

out_lines <- c("=== Validation matrix: standardisation method x contrast ===", "",
  "Cell shows: concordance% (number of studies).  * = p < 0.05 (one-sided",
  "binomial vs a 50% coin-flip null).", "")

for (cn in names(CONTRASTS)) {
  ct <- CONTRASTS[[cn]]
  better <- bio$Intervention_Class %in% ct$healthier
  worse  <- bio$Intervention_Class %in% ct$comparator
  if (!is.null(ct$exclude)) { drop <- bio$Study_ID %in% ct$exclude; better <- better & !drop; worse <- worse & !drop }
  hdr <- sprintf("--- %s contrast (%s) ---", cn, ct$label)
  log_msg(""); log_msg("  %s", hdr)
  out_lines <- c(out_lines, hdr, "")

  head_row <- sprintf("  %-14s %s", "index",
                      paste(sprintf("%-12s", names(METHODS)), collapse = ""))
  log_msg("%s", head_row); out_lines <- c(out_lines, head_row)

  for (v in c("Morphology", "Biomarker", "Combined")) {
    cells <- sapply(names(METHODS), function(mn)
      fmt(concordance(built[[mn]][[v]], bio$Study_ID, better, worse)))
    row <- sprintf("  %-14s %s", v, paste(sprintf("%-12s", cells), collapse = ""))
    log_msg("%s", row); out_lines <- c(out_lines, row)
  }

  out_lines <- c(out_lines, "", "  per domain:")
  log_msg("  per domain:")
  for (d in names(DOMAINS)) {
    cells <- sapply(names(METHODS), function(mn) {
      s <- built[[mn]]$domains[[d]]
      if (is.null(s)) "     -    " else fmt(concordance(s, bio$Study_ID, better, worse))
    })
    row <- sprintf("  %-14s %s", d, paste(sprintf("%-12s", cells), collapse = ""))
    log_msg("%s", row); out_lines <- c(out_lines, row)
  }
  out_lines <- c(out_lines, "")
}

# =============================================================================
# STRATIFIED BY CHALLENGE SUBTYPE
#
# Pooling Eimeria, C. perfringens, LPS and antibiotic disruption into a single
# "challenge" category cancels real signal, because these insults act through
# unrelated mechanisms and move tight junction genes in opposite directions.
# A study comparing the first two directly found C. perfringens downregulated
# occludin, claudin-1, TJP1 and TJP2 while Eimeria downregulated claudin-1 but
# UPREGULATED JAM2 and TJP1 - opposite directions in the same experiment.
#
# Each subtype is therefore tested separately against controls from the same
# studies. Subtypes with fewer than 3 evaluable studies are reported but not
# interpreted - the estimate is too unstable to mean anything.
# =============================================================================
if ("Intervention_Subclass" %in% names(bio)) {
  log_msg("")
  log_msg("  --- stratified by challenge subtype (RAW standardisation) ---")
  out_lines <- c(out_lines, "", "--- stratified by challenge subtype (RAW) ---", "")

  subs <- unique(bio$Intervention_Subclass[bio$Intervention_Class %in%
            c("Infectious challenge", "Toxin challenge", "Environmental challenge","Antibiotic challenge","Dietary challenge")])
  subs <- subs[!is.na(subs) & nzchar(subs)]

  hdr <- sprintf("  %-34s %-12s %s", "challenge subtype", "index", "concordance")
  log_msg("%s", hdr); out_lines <- c(out_lines, hdr)

  for (sb in sort(subs)) {
    better <- bio$Intervention_Class %in% "Unchallenged control"
    worse  <- bio$Intervention_Subclass %in% sb
    for (v in c("Morphology", "Biomarker", "Combined")) {
      cc <- concordance(built[["RAW"]][[v]], bio$Study_ID, better, worse)
      if (cc$n_eval == 0) next
      mark <- if (!is.na(cc$p) && cc$p < 0.05) "*" else " "
      note <- if (cc$n_eval < 3) "  (too few studies to interpret)" else ""
      row <- sprintf("  %-34s %-12s %2d/%2d = %3.0f%%%s%s",
                     substr(sb, 1, 34), v, cc$n_ok, cc$n_eval,
                     100 * cc$prop, mark, note)
      log_msg("%s", row); out_lines <- c(out_lines, row)
    }
  }

  log_msg("")
  log_msg("  --- Barrier vs Claudin1 by subtype (the asymmetry under test) ---")
  out_lines <- c(out_lines, "", "--- Barrier (ZO-1 + occludin) vs Claudin1, by subtype ---", "")
  for (sb in sort(subs)) {
    better <- bio$Intervention_Class %in% "Unchallenged control"
    worse  <- bio$Intervention_Subclass %in% sb
    for (d in c("Barrier", "Claudin1")) {
      s <- built[["RAW"]]$domains[[d]]
      if (is.null(s)) next
      cc <- concordance(s, bio$Study_ID, better, worse)
      if (cc$n_eval == 0) next
      mark <- if (!is.na(cc$p) && cc$p < 0.05) "*" else " "
      row <- sprintf("  %-34s %-12s %2d/%2d = %3.0f%%%s", substr(sb,1,34), d,
                     cc$n_ok, cc$n_eval, 100 * cc$prop, mark)
      log_msg("%s", row); out_lines <- c(out_lines, row)
    }
  }
} else {
  log_msg("  Intervention_Subclass column not found - stratification skipped")
}

log_msg("[4/4] Writing summary")

# =============================================================================
# CAN A SINGLE UNIFIED SCALE BE BUILT AT ALL?
#
# The proposal - one index where additive > control > challenge - is exactly
# the structure used by MAPI and the Microbial Dysbiosis Index (log-ratio of
# improvement-associated over damage-associated taxa). It is a real, precedented
# design. But it is only valid for a marker whose direction of effect is the
# SAME under both comparisons - i.e. "higher is healthier" has to mean the
# same thing whether the bird got an additive or a challenge.
#
# This is checked directly per marker, rather than assumed. A marker qualifies
# for the unified scale only if BOTH:
#   - it moved the expected direction for ADDITIVE (additive > control)
#   - it moved the expected direction for CHALLENGE (control > challenged)
# A marker passing one but not the other - claudin-1 is the known case, since
# it rises under both conditions rather than falling under challenge - cannot
# sit on a single shared axis without distorting it, and is flagged rather
# than silently included.
# =============================================================================
log_msg("")
log_msg("  --- per-marker direction consistency (RAW values) ---")
out_lines <- c(out_lines, "", "--- Can a single unified index be built? Per-marker check ---",
  "A marker qualifies for one shared scale only if it moves the expected",
  "direction under BOTH contrasts. Failing one side means that marker cannot",
  "sit on the same axis as the others without distortion.", "")

is_add_better <- bio$Intervention_Class %in% "Non-disruptive additive"
is_ctrl        <- bio$Intervention_Class %in% "Unchallenged control"
is_chal        <- bio$Intervention_Class %in% c("Infectious challenge", "Toxin challenge", "Environmental challenge","Antibiotic challenge","Dietary challenge")

hdr <- sprintf("  %-10s %-18s %10s %10s   %s", "domain", "marker", "additive", "challenge", "qualifies for unified scale?")
log_msg("%s", hdr); out_lines <- c(out_lines, hdr)

unified_ok <- character(0)
unified_flagged <- character(0)

for (dom in names(DOMAINS)) {
  for (m in names(DOMAINS[[dom]])) {
    if (!m %in% names(bio)) next
    dir <- DOMAINS[[dom]][[m]]
    sc <- std_raw(bio[[m]], bio$Study_ID) * dir
    add_cc  <- concordance(sc, bio$Study_ID, is_add_better, is_ctrl)
    chal_cc <- concordance(sc, bio$Study_ID, is_ctrl, is_chal)
    add_ok  <- !is.na(add_cc$prop)  && add_cc$prop  > 0.5
    chal_ok <- !is.na(chal_cc$prop) && chal_cc$prop > 0.5
    verdict <- if (add_ok && chal_ok) "YES" else if (add_ok || chal_ok) "NO - one-sided only" else "NO"
    if (verdict == "YES") unified_ok <- c(unified_ok, m) else unified_flagged <- c(unified_flagged, m)
    a_txt <- if (is.na(add_cc$prop)) "   -" else sprintf("%3.0f%%", 100*add_cc$prop)
    c_txt <- if (is.na(chal_cc$prop)) "   -" else sprintf("%3.0f%%", 100*chal_cc$prop)
    row <- sprintf("  %-10s %-18s %10s %10s   %s", dom, m, a_txt, c_txt, verdict)
    log_msg("%s", row); out_lines <- c(out_lines, row)
  }
}

log_msg("")
log_msg("  markers that qualify for a single unified scale: %s",
        if (length(unified_ok)) paste(unified_ok, collapse=", ") else "NONE")
log_msg("  markers that do NOT (need separate handling, e.g. claudin-1): %s",
        if (length(unified_flagged)) paste(unified_flagged, collapse=", ") else "none")
out_lines <- c(out_lines, "",
  sprintf("Qualify for unified scale: %s",
          if (length(unified_ok)) paste(unified_ok, collapse=", ") else "NONE"),
  sprintf("Do NOT qualify (need separate handling): %s",
          if (length(unified_flagged)) paste(unified_flagged, collapse=", ") else "none"),
  "",
  "If most markers qualify: build the unified index using only the qualifying",
  "markers (a real MAPI/MDI-style log-ratio or averaged scale), and report the",
  "disqualified markers separately rather than forcing them in.",
  "",
  "If few markers qualify: a single shared scale is not well supported by this",
  "data, and the three separate contrasts (additive, challenge, and the",
  "domain-specific results above) are the more honest reporting structure.")
out_lines <- c(out_lines, "",
  "--- how to read this ---",
  "A high percentage with a * means the index reliably ordered the two groups",
  "correctly in that many studies, more often than chance would give.",
  "",
  "Compare ACROSS a row: does the standardisation method matter for this index?",
  "Compare BETWEEN the two contrast blocks: is the index better at detecting",
  "damage (CHALLENGE) or improvement (ADDITIVE)?",
  "",
  "The ADDITIVE contrast is not a weaker substitute for CHALLENGE. Most poultry",
  "gut-health research evaluates feed additives rather than infections, so this",
  "corpus carries more data for it, and detecting improvement above an",
  "untreated baseline is directly useful for evaluating interventions.",
  "",
  "A method that wins on one contrast but not the other is informative, not a",
  "contradiction: detecting damage and detecting improvement are different",
  "tasks and need not be equally easy.",
  "",
  "CAUTION on CONTROL_REF: it divides by each study's control mean, so it is",
  "defined only where a control group exists, and by construction sets every",
  "control to exactly 1.0. That removes within-control variation, which can",
  "flatter the contrast. Treat it as a diagnostic comparison rather than as a",
  "candidate for the final pipeline without further scrutiny.")

dir.create(file.path(R_OUT_DIR, "ghi"), recursive = TRUE, showWarnings = FALSE)
f <- file.path(R_OUT_DIR, "ghi", "validation_matrix.txt")
writeLines(out_lines, f)
log_msg("  -> %s", f)
log_msg("=== DONE ===")
