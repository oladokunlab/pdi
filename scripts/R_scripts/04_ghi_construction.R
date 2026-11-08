#!/usr/bin/env Rscript
# =============================================================================
# 04_ghi_construction.R
#
# Builds the Gut Health Index (GHI) from measured biomarkers, using a
# hierarchical (domain-based) composite structure.
#
# HOW TO RUN: from a terminal (Rscript 04_ghi_construction.R) or inside
# RStudio (open the file, click Source). No path editing required; see
# 00_r_config.R for the expected repository layout.
#
# -----------------------------------------------------------------------------
# METHOD AND ITS PRECEDENT
# -----------------------------------------------------------------------------
# Individual biomarkers are grouped into biological domains; a sub-score is
# computed within each domain from whichever markers that sample has; and the
# domain sub-scores are then combined into one overall index. This is the
# standard construction used for composite indices where the underlying
# indicators are not all measured on every unit - the Human Development Index
# (UNDP) is built this way (health, education, and living-standards
# sub-indices computed separately, then aggregated), precisely because
# component data come from different sources with different coverage.
#
# The alternative - requiring every biomarker to be present in the same
# sample - is complete-case analysis. It is rejected here for a specific
# reason: in this corpus, whether a biomarker was measured depends on what
# kind of study it was (challenge studies measure cytokines; nutrition
# studies measure morphology). Missingness is therefore informative, not
# random, which is the textbook condition under which complete-case analysis
# introduces selection bias. Requiring joint completeness across the full
# panel also yields an empty dataset in this corpus.
#
# Data completeness is reported per sample (GHI_n_domains, GHI_n_markers)
# rather than hidden, following the convention that a composite score should
# travel with a record of how much data it was actually built from.
#
# -----------------------------------------------------------------------------
# WHAT THIS SCRIPT PRODUCES
# -----------------------------------------------------------------------------
#   GHI_Morphology  - morphology domain only (villus height, crypt depth)
#   GHI_Biomarker   - mean of the non-morphology domains present
#   GHI_Combined    - mean of ALL domains present, morphology included
# plus per-sample domain sub-scores and completeness counts.
#
# The three outputs correspond to the three index versions to be validated
# separately; the domain structure is what makes the second and third
# computable at all, since neither is obtainable under joint completeness.
# =============================================================================

# ---- Locate this script's folder and load shared configuration -------------
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
# CONFIGURATION - REVIEW THIS SECTION BEFORE RUNNING
# =============================================================================
#
# DIRECTION OF EFFECT (+1 / -1) determines whether a HIGH value of a marker
# is treated as BETTER (+1) or WORSE (-1) gut health. Every marker is
# multiplied by its direction so that, after orientation, a higher score
# always means better gut health. Getting a direction wrong does not cause
# an error - it silently inverts that marker's contribution - so these are
# listed explicitly here rather than assumed anywhere in the code.
#
# Directions marked REVIEW below are genuinely debatable and should be
# confirmed before the index is used for any published result.

DOMAINS <- list(
  Morphology = c(
    VH_um  = +1,   # villus height: taller villi = greater absorptive surface
    CD_um  = -1,   # crypt depth: deeper crypts = higher epithelial turnover,
                   #   conventionally interpreted as a damage/stress response
    VH_CD  = +1    # villus height : crypt depth ratio
  ),
  Barrier = c(
    ZO1    = +1,   # zonula occludens-1: tight junction integrity
    Ocln   = +1,   # occludin: tight junction integrity
    Cldn1  = +1,   # claudin-1: tight junction integrity
    MUC2   = +1,   # mucin-2: protective mucus layer
    FITC   = -1    # FITC-dextran translocation: higher = leakier gut
  ),
  # Inflammation is split into two domains rather than one, because this
  # corpus reports cytokines in two biologically distinct compartments:
  # circulating serum concentration versus local tissue/mucosal expression
  # (protein-normalized or relative mRNA). These are not interchangeable
  # scales of the same signal. The compartment-splitting step above
  # (COMPARTMENT SPLITTING, run before this config is evaluated) creates the
  # _FLUID/_TISSUE columns referenced here automatically from the raw data;
  # a marker that turned out to be single-compartment in this corpus
  # (confirmed by audit: IgA, IgY, IL10) is referenced by its plain name.
  Inflammation_Systemic = c(
    IL1b_FLUID = -1,   # pro-inflammatory, serum
    IL6_FLUID  = -1,   # pro-inflammatory, serum
    TNFa_FLUID = -1,   # pro-inflammatory, serum
    IFNg_FLUID = -1,   # pro-inflammatory, serum
    IL10       = +1    # anti-inflammatory / regulatory, serum (confirmed
                        #   single-compartment in this corpus)
  ),
  Inflammation_Mucosal = c(
    IL1b_TISSUE = -1,  # pro-inflammatory, tissue/mucosal
    IL6_TISSUE  = -1,  # pro-inflammatory, tissue/mucosal
    TNFa_TISSUE = -1,  # pro-inflammatory, tissue/mucosal
    IFNg_TISSUE = -1   # pro-inflammatory, tissue/mucosal
  ),
  # Immunoglobulin is split the same way and for the same reason: systemic
  # (serum) antibody status reflects a different biological signal than
  # local mucosal antibody status, and averaging them would make the domain
  # mean different things depending on which compartment a given sample
  # happened to have measured.
  Immunoglobulin_Systemic = c(
    sIgA_FLUID = +1,   # secretory IgA, serum
    IgA        = +1,   # confirmed single-compartment (serum) in this corpus
    IgG_FLUID  = +1,   # REVIEW: higher may indicate immune competence (+1)
                        #   OR an active systemic response to challenge (-1).
                        #   Direction is not settled in the literature and
                        #   depends on context.
    IgM_FLUID  = +1,   # REVIEW: same ambiguity as IgG
    IgY        = +1    # REVIEW: same ambiguity as IgG; confirmed
                        #   single-compartment (serum), but sparse in this
                        #   corpus, so its influence is minimal either way
  ),
  Immunoglobulin_Mucosal = c(
    sIgA_TISSUE = +1,  # secretory IgA: mucosal immune competence
    IgG_TISSUE  = +1,  # REVIEW: same directional ambiguity as serum IgG,
                        #   evaluated locally rather than systemically
    IgM_TISSUE  = +1   # REVIEW: same ambiguity as mucosal IgG
  ),
  SCFA = c(
    Ace    = +1,   # acetate: microbial fermentation product
    Prop   = +1,   # propionate
    Buty   = +1    # butyrate: primary colonocyte energy source
  )
)

# Domains counted as "biomarker" (i.e. everything except morphology) for the
# GHI_Biomarker version.
MORPHOLOGY_DOMAIN <- "Morphology"

# ---- Z-SCORE SOURCE --------------------------------------------------------
# "existing_first" : use the workbook's pre-computed *_zscore column where one
#                    exists, and compute within-study z-scores only for markers
#                    lacking one. CONSERVATIVE - the workbook's cytokine
#                    z-scores were computed only on values successfully
#                    converted to a common unit (pg/mL), so values recorded in
#                    other units are excluded.
# "recompute_all"  : ignore the pre-computed columns and z-score every marker
#                    within study from the raw values. RECOVERS MORE DATA -
#                    within-study standardisation does not require cross-study
#                    unit conversion, provided each study uses one consistent
#                    unit internally. Verify that assumption with the unit
#                    audit written by this script before relying on it.
#
# This choice materially changes sample size for the cytokine markers (roughly
# a two-fold difference in this corpus) and is therefore exposed here rather
# than hard-coded.
ZSCORE_SOURCE <- "existing_first"

# A domain sub-score is computed for a sample if at least this many of the
# domain's markers are present. Set to 1 so that a domain contributes
# whenever any of its markers was measured; raising it makes each domain
# sub-score more robust but reduces sample coverage.
MIN_MARKERS_PER_DOMAIN <- 1

# An overall GHI is reported for a sample only if at least this many domains
# have a sub-score. Set to 1 to maximise coverage; the per-sample
# GHI_n_domains column records the actual number so that low-coverage
# samples can be filtered or down-weighted downstream.
MIN_DOMAINS_FOR_GHI <- 1

# =============================================================================
# STEP 1 - LOAD BIOMARKER DATA
# =============================================================================
log_msg("[1/6] Loading biomarker sheet '%s' from %s", BIOMARKER_SHEET, METADATA_XLSX)
if (!file.exists(METADATA_XLSX)) {
  stop(sprintf("Biomarker workbook not found at %s", METADATA_XLSX))
}
bio <- as.data.frame(read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
log_msg("  %d rows x %d columns", nrow(bio), ncol(bio))

for (req in c("Study_ID", "Sample_ID")) {
  if (!req %in% names(bio)) stop(sprintf("Required column '%s' not found in %s", req, BIOMARKER_SHEET))
}

# =============================================================================
# COMPARTMENT SPLITTING (fluid/serum vs tissue/mucosal)
#
# Several markers in this workbook are measured in more than one biological
# compartment across the corpus - e.g. a cytokine reported as a serum
# concentration (pg/mL) in some studies and as a tissue-homogenate,
# protein-normalized value (ng/mg protein) or relative gene expression
# (fold-change) in others. These are not different units of the same
# measurement; they are different biological signals (circulating protein
# vs local tissue transcription) that happen to share a marker name.
#
# Averaging them under one domain score would make that domain mean
# different things for different samples, which defeats the purpose of a
# domain-based composite index (see script header). Every marker with a
# *_Unit column is therefore classified automatically, per value, as FLUID
# (concentration per volume - serum/plasma) or TISSUE (protein-normalized or
# relative expression - mucosal/tissue), and split into two columns before
# any z-scoring happens. This runs on the whole corpus automatically, so a
# newly added study with the same ambiguity is handled without manual
# spreadsheet edits.
#
# A value whose unit does not match either pattern is left untouched under
# its original column name and reported below - if new unrecognized units
# appear (e.g. a new study with a unit variant not yet seen), they will show
# up here rather than being silently miscategorised.
# =============================================================================
FLUID_UNIT_SUFFIXES <- c("/ml", "/l")
TISSUE_UNIT_MARKERS  <- c("mgprot", "mg protein", "relativemrna", "foldchange", "/mg")

classify_unit <- function(u) {
  if (is.na(u) || !nzchar(trimws(u))) return(NA_character_)
  n <- tolower(gsub("[[:space:]]", "", u))
  n <- gsub("\u00b5|\u03bc", "u", n)   # normalise micro-sign to "u"
  if (any(sapply(TISSUE_UNIT_MARKERS, function(p) grepl(gsub(" ", "", p), n, fixed = TRUE))) ||
      grepl("fold", n)) {
    return("TISSUE")
  }
  if (any(sapply(FLUID_UNIT_SUFFIXES, function(s) endsWith(n, s)))) {
    return("FLUID")
  }
  NA_character_
}

log_msg("[1b/6] Auditing marker units for fluid/tissue compartment mixing")

unit_cols <- grep("_Unit$", names(bio), value = TRUE)
compartment_report <- list()
split_any <- FALSE

for (uc in unit_cols) {
  m <- sub("_Unit$", "", uc)
  if (!(m %in% names(bio))) next

  cls <- vapply(bio[[uc]], classify_unit, character(1))
  has_val <- !is.na(bio[[m]])

  fluid_studies  <- unique(bio$Study_ID[has_val & !is.na(cls) & cls == "FLUID"])
  tissue_studies <- unique(bio$Study_ID[has_val & !is.na(cls) & cls == "TISSUE"])
  unknown_rows   <- has_val & is.na(cls) & !is.na(bio[[uc]]) & nzchar(trimws(bio[[uc]]))

  if (length(fluid_studies) > 0 && length(tissue_studies) > 0) {
    split_any <- TRUE
    fcol <- paste0(m, "_FLUID"); tcol <- paste0(m, "_TISSUE")
    bio[[fcol]] <- ifelse(has_val & !is.na(cls) & cls == "FLUID", bio[[m]], NA_real_)
    bio[[tcol]] <- ifelse(has_val & !is.na(cls) & cls == "TISSUE", bio[[m]], NA_real_)
    compartment_report[[m]] <- list(
      fluid_studies = fluid_studies, tissue_studies = tissue_studies,
      fluid_col = fcol, tissue_col = tcol
    )
    log_msg("  %-8s SPLIT: -> %s (n studies=%d), %s (n studies=%d)",
            m, fcol, length(fluid_studies), tcol, length(tissue_studies))
  }

  if (any(unknown_rows)) {
    bad_units <- unique(bio[[uc]][unknown_rows])
    log_msg("  *** %-8s has unrecognised unit string(s), left unsplit: %s",
            m, paste(bad_units, collapse = "; "))
  }
}
if (!split_any) log_msg("  no fluid/tissue compartment mixing detected")

# =============================================================================
# STEP 2 - WITHIN-STUDY Z-SCORING
#
# Markers are standardised WITHIN each study, not across the whole corpus.
# This is essential in a meta-analysis: absolute biomarker values are not
# comparable across studies (different assay kits, units, reference ranges,
# bird ages), so a raw value that is high in one study may be ordinary in
# another. Within-study z-scoring expresses each sample relative to its own
# study's samples, which is what makes values poolable at all.
#
# Some markers already have a *_zscore column in the workbook; those are
# verified to be within-study and reused. Markers without one are z-scored
# here using the same within-study convention.
#
# A marker cannot be z-scored within a study that has fewer than 2 samples
# with that marker, or where all values are identical (zero variance). Those
# cases yield NA rather than a fabricated value, and are counted in the log.
# =============================================================================
log_msg("[2/6] Computing within-study z-scores")

all_markers <- unlist(lapply(DOMAINS, names), use.names = FALSE)
missing_from_sheet <- setdiff(all_markers, names(bio))
if (length(missing_from_sheet) > 0) {
  log_msg("  NOTE: %d configured marker(s) absent from the sheet and skipped: %s",
          length(missing_from_sheet), paste(missing_from_sheet, collapse = ", "))
}
usable_markers <- intersect(all_markers, names(bio))

# Unit audit: within-study z-scoring is only valid if each study uses one
# consistent unit for a given marker internally. Any study mixing units for a
# marker is reported here, because its z-scores would compare values that are
# not on the same scale.
unit_problems <- character(0)
for (m in usable_markers) {
  ucol <- paste0(m, "_Unit")
  if (!ucol %in% names(bio)) next
  for (s in unique(bio$Study_ID)) {
    sel <- bio$Study_ID == s & !is.na(bio[[m]])
    if (sum(sel) == 0) next
    us <- unique(stats::na.omit(bio[[ucol]][sel]))
    us <- us[nzchar(trimws(us))]
    if (length(us) > 1) {
      unit_problems <- c(unit_problems,
        sprintf("%s / %s: %s", s, m, paste(us, collapse = " | ")))
    }
  }
}
if (length(unit_problems) > 0) {
  log_msg("  *** UNIT WARNING: %d study/marker combination(s) mix units internally.", length(unit_problems))
  log_msg("  *** Within-study z-scores for these compare values on different scales:")
  for (u in unit_problems) log_msg("  ***   %s", u)
} else {
  log_msg("  unit audit: no study mixes units within a marker")
}

zscore_within_study <- function(values, study) {
  out <- rep(NA_real_, length(values))
  for (s in unique(study)) {
    sel <- study == s & !is.na(values)
    v <- values[sel]
    if (length(v) < 2) next          # cannot standardise a single observation
    sdv <- stats::sd(v)
    if (!is.finite(sdv) || sdv == 0) next   # zero variance: no meaningful z
    out[sel] <- (v - mean(v)) / sdv
  }
  out
}

z_used_existing <- character(0)
z_computed     <- character(0)
for (m in usable_markers) {
  existing <- paste0(m, "_zscore")
  target   <- paste0("Z_", m)
  use_existing <- (ZSCORE_SOURCE == "existing_first") &&
                  existing %in% names(bio) &&
                  sum(!is.na(bio[[existing]])) > 0
  if (use_existing) {
    bio[[target]] <- suppressWarnings(as.numeric(bio[[existing]]))
    z_used_existing <- c(z_used_existing, m)
  } else {
    raw <- suppressWarnings(as.numeric(bio[[m]]))
    bio[[target]] <- zscore_within_study(raw, bio$Study_ID)
    z_computed <- c(z_computed, m)
  }
}
log_msg("  reused existing z-scores for %d marker(s): %s",
        length(z_used_existing), paste(z_used_existing, collapse = ", "))
log_msg("  computed new within-study z-scores for %d marker(s): %s",
        length(z_computed), paste(z_computed, collapse = ", "))

# Report markers where z-scoring lost data relative to raw availability,
# which happens when a study contributed too few samples to standardise.
for (m in z_computed) {
  n_raw <- sum(!is.na(suppressWarnings(as.numeric(bio[[m]]))))
  n_z   <- sum(!is.na(bio[[paste0("Z_", m)]]))
  if (n_raw > n_z) {
    log_msg("    %-6s %d raw values -> %d z-scored (%d in studies with <2 values or zero variance)",
            m, n_raw, n_z, n_raw - n_z)
  }
}

# =============================================================================
# STEP 3 - ORIENT MARKERS SO THAT HIGHER ALWAYS MEANS BETTER
# Each z-score is multiplied by its configured direction (+1 / -1).
# =============================================================================
log_msg("[3/6] Orienting markers by direction of effect")
n_flipped <- 0
for (dom in names(DOMAINS)) {
  for (m in names(DOMAINS[[dom]])) {
    if (!m %in% usable_markers) next
    dir <- DOMAINS[[dom]][[m]]
    bio[[paste0("O_", m)]] <- bio[[paste0("Z_", m)]] * dir
    if (dir < 0) n_flipped <- n_flipped + 1
  }
}
log_msg("  %d marker(s) inverted (higher raw value = worse gut health)", n_flipped)

# =============================================================================
# STEP 4 - DOMAIN SUB-SCORES
# Each domain sub-score is the mean of that sample's available oriented
# z-scores within the domain. Using the mean of whatever is present (rather
# than requiring all markers) is what allows partially-measured samples to
# contribute, and keeps the sub-score on the same scale regardless of how
# many markers underlie it.
# =============================================================================
log_msg("[4/6] Computing domain sub-scores (>= %d marker(s) required per domain)",
        MIN_MARKERS_PER_DOMAIN)

domain_score_cols <- character(0)
for (dom in names(DOMAINS)) {
  cols <- paste0("O_", intersect(names(DOMAINS[[dom]]), usable_markers))
  if (length(cols) == 0) {
    log_msg("  %-16s SKIPPED - no markers available", dom)
    next
  }
  M <- as.matrix(bio[, cols, drop = FALSE])
  n_present <- rowSums(!is.na(M))
  score <- ifelse(n_present >= MIN_MARKERS_PER_DOMAIN, rowMeans(M, na.rm = TRUE), NA_real_)
  score[!is.finite(score)] <- NA_real_
  scol <- paste0("GHI_dom_", dom)
  ncol_ <- paste0("GHI_nmark_", dom)
  bio[[scol]]  <- score
  bio[[ncol_]] <- n_present
  domain_score_cols <- c(domain_score_cols, scol)
  log_msg("  %-16s %4d samples scored (from %d marker column(s))",
          dom, sum(!is.na(score)), length(cols))
}

# =============================================================================
# STEP 5 - COMBINE DOMAINS INTO THE THREE INDEX VERSIONS
#
# Domains are combined by arithmetic mean of the available domain sub-scores.
# Note this differs deliberately from the Human Development Index, which uses
# a geometric mean: a geometric mean requires strictly positive values, and
# these sub-scores are z-score based and therefore centred on zero and often
# negative. An arithmetic mean is the appropriate analogue on this scale.
# =============================================================================
log_msg("[5/6] Combining domains (>= %d domain(s) required for an overall score)",
        MIN_DOMAINS_FOR_GHI)

combine_domains <- function(dt, cols) {
  if (length(cols) == 0) return(list(score = rep(NA_real_, nrow(dt)), n = rep(0L, nrow(dt))))
  M <- as.matrix(dt[, cols, drop = FALSE])
  n <- rowSums(!is.na(M))
  s <- ifelse(n >= MIN_DOMAINS_FOR_GHI, rowMeans(M, na.rm = TRUE), NA_real_)
  s[!is.finite(s)] <- NA_real_
  list(score = s, n = n)
}

morph_col      <- paste0("GHI_dom_", MORPHOLOGY_DOMAIN)
morph_cols     <- intersect(morph_col, domain_score_cols)
biomarker_cols <- setdiff(domain_score_cols, morph_col)

res_morph <- combine_domains(bio, morph_cols)
res_bio   <- combine_domains(bio, biomarker_cols)
res_comb  <- combine_domains(bio, domain_score_cols)

bio$GHI_Morphology <- res_morph$score
bio$GHI_Biomarker  <- res_bio$score
bio$GHI_Combined   <- res_comb$score
bio$GHI_n_domains  <- res_comb$n

marker_count_cols <- grep("^GHI_nmark_", names(bio), value = TRUE)
bio$GHI_n_markers <- rowSums(as.matrix(bio[, marker_count_cols, drop = FALSE]), na.rm = TRUE)

log_msg("  GHI_Morphology : %4d samples", sum(!is.na(bio$GHI_Morphology)))
log_msg("  GHI_Biomarker  : %4d samples", sum(!is.na(bio$GHI_Biomarker)))
log_msg("  GHI_Combined   : %4d samples", sum(!is.na(bio$GHI_Combined)))

# =============================================================================
# STEP 6 - EXPORT
# =============================================================================
log_msg("[6/6] Writing output")

keep_cols <- c(
  intersect(c("Study_ID", "Sample_ID", "Treatment_Group", "Segment", "Age_d"), names(bio)),
  domain_score_cols, marker_count_cols,
  "GHI_Morphology", "GHI_Biomarker", "GHI_Combined",
  "GHI_n_domains", "GHI_n_markers"
)
out_scores <- file.path(R_OUT_DIR, "ghi", "ghi_scores.tsv")
fwrite(bio[, keep_cols, drop = FALSE], out_scores, sep = "\t")
log_msg("  per-sample GHI scores -> %s", out_scores)

# Coverage breakdown, so the completeness behind each version is on record.
cov_lines <- c(
  "=== GHI construction summary ===",
  "",
  sprintf("Input rows: %d", nrow(bio)),
  sprintf("Minimum markers per domain: %d", MIN_MARKERS_PER_DOMAIN),
  sprintf("Minimum domains for an overall score: %d", MIN_DOMAINS_FOR_GHI),
  "",
  "Domain coverage (samples with a sub-score):",
  sapply(domain_score_cols, function(c_) sprintf("  %-24s %d", sub("GHI_dom_", "", c_), sum(!is.na(bio[[c_]])))),
  "",
  "Index coverage:",
  sprintf("  GHI_Morphology  %d", sum(!is.na(bio$GHI_Morphology))),
  sprintf("  GHI_Biomarker   %d", sum(!is.na(bio$GHI_Biomarker))),
  sprintf("  GHI_Combined    %d", sum(!is.na(bio$GHI_Combined))),
  "",
  "Samples by number of domains contributing to GHI_Combined:",
  sapply(sort(unique(bio$GHI_n_domains)), function(k)
    sprintf("  %d domain(s): %d samples", k, sum(bio$GHI_n_domains == k))),
  "",
  "Markers with a REVIEW note on direction of effect (see script header):",
  "  IgG, IgM, IgY - currently +1 (higher = better). Whether elevated systemic",
  "  immunoglobulin indicates immune competence or an active response to",
  "  challenge is context-dependent and should be confirmed before publication.",
  "",
  "Note: domains are combined by arithmetic mean, not the geometric mean used",
  "by the Human Development Index, because these sub-scores are z-score based",
  "and can be negative."
)
out_summary <- file.path(R_OUT_DIR, "ghi", "ghi_construction_summary.txt")
writeLines(unlist(cov_lines), out_summary)
log_msg("  summary -> %s", out_summary)

log_msg("=== DONE ===")
