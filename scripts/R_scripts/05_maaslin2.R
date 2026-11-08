#!/usr/bin/env Rscript
# =============================================================================
# 05_maaslin2.R
#
# Stage 1 of GDI construction: identify which bacterial genera are robustly
# associated with the Gut Health Index, before any predictive model is fitted.
#
# RUN:  Rscript 05_maaslin2.R
# TIME: a few minutes.
# NEEDS: 03_batch_correction.R and 04b_ghi_pca.R to have run first.
#
# -----------------------------------------------------------------------------
# WHY THIS RUNS BEFORE LASSO
# -----------------------------------------------------------------------------
# The registered protocol specifies LASSO regression "of the GHI on the taxa
# identified as robustly associated with gut-health status" - association
# testing first, penalised regression second. That order matters for two
# reasons beyond protocol fidelity:
#
#   1. LASSO has no false-discovery-rate control. Left to choose freely among
#      all genera it will select background features purely to reduce residual
#      variance. MaAsLin2 applies Benjamini-Hochberg FDR correction first, so
#      only genera with a defensible statistical association reach the model.
#
#   2. Gut genera are heavily collinear. Given a correlated group, LASSO picks
#      one essentially arbitrarily and discards the rest, which produces
#      signatures that do not reproduce across studies. Pre-filtering on
#      per-genus association reduces how much that arbitrariness matters.
#
# -----------------------------------------------------------------------------
# THE ID PROBLEM, AND HOW IT IS HANDLED
# -----------------------------------------------------------------------------
# The genus table is indexed by sequencing run (e.g. SRR10567010). The GHI is
# indexed by biomarker row (e.g. Study_Treatment_Segment_Age). These are
# different identifier spaces and cannot be joined directly.
#
# The link is Study_ID + Treatment_Group + Segment: a sequencing run and a
# biomarker measurement belong together when they come from the same study,
# the same treatment arm, and the same gut segment. That mapping is read from
# the Crosswalk_16S sheet, which carries Study_ID/Treatment_Group/Segment for
# every run.
#
# This join is where silent data loss would occur, so match rates are reported
# explicitly at every step and unmatched runs are listed rather than dropped
# quietly.
#
# A consequence worth stating in the manuscript: several sequencing runs map
# to a single biomarker row, because papers report one group mean for a group
# that was sequenced in replicate. Every run in a treatment arm therefore
# receives that arm's GHI value. This is measurement error in the outcome,
# not in the predictors, and it attenuates the achievable correlation - the
# ceiling on GDI performance is set partly by this, independent of how good
# the microbiome signal is.
# =============================================================================

.locate_and_source_config <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  d <- if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
  source(file.path(d, "00_r_config.R"))
}
.locate_and_source_config()

for (p in c("Maaslin2")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("Package '%s' not installed. BiocManager::install('%s')", p, p))
}

# =============================================================================
# CONFIGURATION
# =============================================================================

# Pre-specified primary outcome, fixed before any sequence data was examined.
# See ghi_pca_summary.txt for the justification. Changing this after seeing
# results below would be selection on the outcome.
GHI_TARGET <- "GHI_Combined"

# How strictly a sequencing run is matched to a biomarker measurement.
#
#   "segment" - run and biomarker must share Study + Treatment arm + Segment.
#               Biologically exact: the microbiome and the phenotype come from
#               the same gut region. But only 45% of runs (21 of 50 studies)
#               satisfy it, because many studies sequenced one segment and
#               measured histology in another.
#
#   "arm"     - run and biomarker must share Study + Treatment arm only.
#               Recovers 86% of runs across 47 studies. The assumption is that
#               a treatment arm's gut-health status is a property of the bird,
#               not of the segment sampled - defensible for systemic markers
#               (serum immunoglobulins, cytokines) and weaker for segment-local
#               ones (villus height in a specific region).
#
# "arm" is the default because leave-one-study-out validation needs studies,
# and 21 is too few to be stable. The stricter join is run as a sensitivity
# analysis. Whichever is used, state it in the manuscript - it materially
# changes the sample.
JOIN_LEVEL <- "arm"     # "arm" or "segment"

# -----------------------------------------------------------------------------
# 16S REGION - ENFORCED UPSTREAM, VERIFIED HERE
#
# The restriction is applied in 01_genus_agglomeration.R, before the
# prevalence filter, the CLR/zero-replacement step and MMUPHin batch
# correction, so that none of those are fitted on runs from other amplicon
# regions. See the block at the top of 01 for the reasoning.
#
# The check below is a guard, not the mechanism: if a non-V3-V4 run reaches
# this point it means 01 was run with a different setting than 05, and the
# two would silently disagree about what the analysis set is. Keep both in
# step, or set both to NA together for the pooled-region sensitivity run.
# -----------------------------------------------------------------------------
REGION_FILTER <- "V3-V4"

# Genera reaching this FDR threshold pass to the LASSO stage.
FDR_THRESHOLD <- 0.25   # MaAsLin2's own default; deliberately permissive at
                        # a screening step, since LASSO applies further
                        # selection afterwards. A stricter 0.05 filter is run
                        # alongside and reported for comparison.

# =============================================================================
# STEP 1 - LOAD THE BATCH-CORRECTED GENUS MATRIX
# =============================================================================
log_msg("[1/6] Loading batch-corrected genus abundances")
clr_path <- file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv")
if (!file.exists(clr_path)) {
  stop(sprintf("Not found: %s\nRun 03_batch_correction.R first.", clr_path))
}
clr <- fread(clr_path, sep = "\t", header = TRUE)
genus_names <- clr[[1]]
clr_mat <- as.matrix(clr[, -1, with = FALSE])
rownames(clr_mat) <- genus_names
log_msg("  %d genera x %d sequencing runs", nrow(clr_mat), ncol(clr_mat))

# =============================================================================
# STEP 2 - MAP SEQUENCING RUNS TO STUDY / TREATMENT / SEGMENT
# =============================================================================
log_msg("[2/6] Mapping runs to study, treatment arm and segment")
cw <- as.data.frame(read_excel(METADATA_XLSX, sheet = "Crosswalk_16S"),
                    stringsAsFactors = FALSE)
need <- c("Run", "Study_ID", "Treatment_Group", "Segment")
missing <- setdiff(need, names(cw))
if (length(missing)) stop(sprintf("Crosswalk_16S missing column(s): %s",
                                  paste(missing, collapse = ", ")))
cw <- cw[, need]
cw$Run <- trimws(as.character(cw$Run))

runs <- colnames(clr_mat)
cw_matched <- cw[match(runs, cw$Run), ]
n_unmapped <- sum(is.na(cw_matched$Study_ID))
if (n_unmapped > 0) {
  bad <- runs[is.na(cw_matched$Study_ID)]
  log_msg("  WARNING: %d run(s) absent from the crosswalk and dropped: %s%s",
          n_unmapped, paste(head(bad, 5), collapse = ", "),
          if (length(bad) > 5) sprintf(" ... (+%d more)", length(bad) - 5) else "")
}
keep <- !is.na(cw_matched$Study_ID)
clr_mat <- clr_mat[, keep, drop = FALSE]
cw_matched <- cw_matched[keep, , drop = FALSE]
log_msg("  %d runs mapped across %d studies",
        ncol(clr_mat), length(unique(cw_matched$Study_ID)))

# ---- 16S region: verify 01 already restricted (see note above) --------------
if (!is.na(REGION_FILTER)) {
  reg <- study_regions()
  sid <- trimws(cw_matched$Study_ID)
  run_region <- reg[sid]
  unknown <- unique(sid[is.na(run_region)])
  if (length(unknown))
    stop("No region in studies.tsv for: ", paste(unknown, collapse = ", "))

  off <- run_region != REGION_FILTER
  if (any(off)) {
    log_msg("  *** %d run(s) from %d non-%s stud%s reached this stage:",
            sum(off), length(unique(sid[off])), REGION_FILTER,
            if (length(unique(sid[off])) == 1) "y" else "ies")
    for (s in sort(unique(sid[off]))) log_msg("      %-22s %s", s, reg[[s]])
    log_msg("  *** 01_genus_agglomeration.R was run with a DIFFERENT region")
    log_msg("  *** setting than this script. Re-run 01 (and 03) so the")
    log_msg("  *** prevalence filter and batch correction match the analysis set.")
    log_msg("  *** Dropping them here would leave the upstream steps fitted on")
    log_msg("  *** data this analysis excludes - a partial restriction.")
    stop("Region setting inconsistent between 01 and 05.")
  }
  log_msg("  region check: all %d runs are %s, consistent with 01",
          ncol(clr_mat), REGION_FILTER)
} else {
  log_msg("  REGION_FILTER is NA - pooled-region sensitivity analysis")
}

# =============================================================================
# STEP 3 - ATTACH THE GHI VIA STUDY + TREATMENT + SEGMENT
# =============================================================================
log_msg("[3/6] Joining GHI at %s level", JOIN_LEVEL)
ghi_path <- file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")
if (!file.exists(ghi_path)) {
  stop(sprintf("Not found: %s\nRun 04b_ghi_pca.R first.", ghi_path))
}
ghi <- fread(ghi_path, sep = "\t", header = TRUE)
if (!GHI_TARGET %in% names(ghi)) {
  stop(sprintf("Target '%s' not present in %s", GHI_TARGET, ghi_path))
}

mk_key <- if (JOIN_LEVEL == "segment") {
  function(s, t, g) paste(trimws(s), trimws(t), trimws(g), sep = "||")
} else {
  function(s, t, g) paste(trimws(s), trimws(t), sep = "||")
}
log_msg("  join level: %s", JOIN_LEVEL)

# -----------------------------------------------------------------------------
# AGE AS AN ADDITIONAL FIXED EFFECT
#
# The registered pipeline test showed age explains far more of the community
# variance than GHI does (PERMANOVA R2 0.12 vs 0.003). Screening genera
# against GHI alone, with age uncontrolled, risks two things: a genus that
# happens to track age gets falsely credited to GHI, and a genus that genuinely
# tracks GHI but is anticorrelated with age in this sample can be masked by it.
# Entering Age_d as a second fixed effect asks MaAsLin2 the sharper question -
# does this genus associate with GHI AFTER accounting for age - rather than
# the marginal, confounded one it was asking before.
#
# AGE_COVARIATE=no restores the original single-fixed-effect model, so the
# two can be compared directly as a sensitivity analysis.
# -----------------------------------------------------------------------------
AGE_COVARIATE <- { e <- tolower(Sys.getenv("AGE_COVARIATE", "")); !(e %in% c("no","false","0")) }

# CLASS_COVARIATE and SEGMENT_COVARIATE: two more fixed effects the guide
# specifies alongside age (Challenge_model, Intestinal_Segment). Segment is
# already carried through this join as an identity column; CLASS_COVARIATE
# controls whether Challenge_Type is entered as a covariate ALONGSIDE GHI, as
# distinct from the stratification done elsewhere (06e/06f/12/16), which
# fits separate models per class. Entering it as a covariate here asks
# whether a genus tracks GHI net of treatment-type differences, using all
# arms in one fit rather than splitting the corpus.
CLASS_COVARIATE   <- { e <- tolower(Sys.getenv("CLASS_COVARIATE", ""));   !(e %in% c("no","false","0")) }
SEGMENT_COVARIATE <- { e <- tolower(Sys.getenv("SEGMENT_COVARIATE", "")); !(e %in% c("no","false","0")) }

ghi_sub <- ghi[!is.na(ghi[[GHI_TARGET]]),
               c("Study_ID", "Treatment_Group", "Segment", GHI_TARGET,
                 intersect("Age_d", names(ghi)),
                 intersect("Challenge_Type", names(ghi))), with = FALSE]
ghi_sub$.key <- mk_key(ghi_sub$Study_ID, ghi_sub$Treatment_Group, ghi_sub$Segment)
# At arm level several segments collapse to one key, so duplicates are the
# norm rather than an anomaly; they are averaged below and the count reported.

# One biomarker row per key is expected. Where a study reports the same arm
# and segment more than once (e.g. two ages), the values are averaged and the
# collapse is reported rather than silently taking the first.
has_age <- "Age_d" %in% names(ghi_sub)
has_cls <- "Challenge_Type" %in% names(ghi_sub)
# most-common-value collapse for the class label, matching how duplicate keys
# are handled elsewhere in this project (mode rather than mean, since this is
# categorical)
mode1 <- function(v) { v <- v[!is.na(v) & nzchar(v)]
  if (!length(v)) return(NA_character_); names(sort(table(v), decreasing = TRUE))[1] }
dupes <- sum(duplicated(ghi_sub$.key))
if (dupes > 0) {
  log_msg("  NOTE: %d duplicate study/arm/segment key(s) in the GHI table;", dupes)
  log_msg("        averaging their %s values.", GHI_TARGET)
  agg <- list(val = quote(mean(get(GHI_TARGET), na.rm = TRUE)))
  if (has_age) agg$age <- quote(mean(suppressWarnings(as.numeric(Age_d)), na.rm = TRUE))
  if (has_cls) agg$cls <- quote(mode1(as.character(Challenge_Type)))
  ghi_sub <- ghi_sub[, eval(as.call(c(quote(list), agg))), by = .key]
} else {
  agg <- list(.key = quote(.key), val = quote(get(GHI_TARGET)))
  if (has_age) agg$age <- quote(suppressWarnings(as.numeric(Age_d)))
  if (has_cls) agg$cls <- quote(as.character(Challenge_Type))
  ghi_sub <- ghi_sub[, eval(as.call(c(quote(list), agg)))]
}

run_key <- mk_key(cw_matched$Study_ID, cw_matched$Treatment_Group, cw_matched$Segment)
ghi_val <- ghi_sub$val[match(run_key, ghi_sub$.key)]
age_val <- if (has_age) ghi_sub$age[match(run_key, ghi_sub$.key)] else rep(NA_real_, length(run_key))
cls_val <- if (has_cls) ghi_sub$cls[match(run_key, ghi_sub$.key)] else rep(NA_character_, length(run_key))
if (AGE_COVARIATE && !has_age)
  log_msg("  NOTE: Age_d not present in the GHI table - AGE_COVARIATE requested but skipped")
if (CLASS_COVARIATE && !has_cls)
  log_msg("  NOTE: Challenge_Type not present - CLASS_COVARIATE requested but skipped")

# -----------------------------------------------------------------------------
# COARSE-LABEL FALLBACK
#
# Some studies deposited sequences under a coarser treatment label than the
# one used to report phenotype. Jansseune_2024 is the case here: all six arms
# (Ctrl/Pro/Post x Standard/Challenge diet) were sequenced and all six have
# biomarker data, but the deposit records only the diet - the run-level
# Ctrl/Pro/Post assignment is not in the public record, and the sample aliases
# are opaque internal lab IDs.
#
# Rather than drop 249 runs, or invent an arm assignment that would attach the
# control group's phenotype to probiotic and postbiotic birds, any run whose
# own label finds no match is given the MEAN GHI of every arm in that study
# whose label CONTAINS the run's label as a substring. For a run labelled
# "Standard diet", that averages Ctrl SD, Pro SD and Post SD.
#
# What this asserts is only what is actually known: this bird was on the
# standard diet, and standard-diet birds collectively scored X. It does not
# claim to know which treatment the bird received.
#
# The cost is precision: a diet-level mean is a noisier outcome than an
# arm-level value, which attenuates any association involving these runs.
# Runs assigned this way are flagged in the output so the analysis can be
# repeated without them as a sensitivity check.
# -----------------------------------------------------------------------------
coarse_assigned <- rep(FALSE, length(ghi_val))
if (any(is.na(ghi_val))) {
  ghi_arm <- ghi[!is.na(ghi[[GHI_TARGET]]),
                 c("Study_ID", "Treatment_Group", GHI_TARGET), with = FALSE]
  for (i in which(is.na(ghi_val))) {
    s  <- trimws(cw_matched$Study_ID[i])
    lab <- trimws(cw_matched$Treatment_Group[i])
    if (is.na(lab) || !nzchar(lab)) next
    cand <- ghi_arm[trimws(ghi_arm$Study_ID) == s, ]
    if (!nrow(cand)) next
    # Match the diet token in either direction: the run label may be a
    # substring of the arm label ("Standard diet" vs "Ctrl SD") or vice versa.
    toks <- unlist(strsplit(tolower(lab), "[^a-z0-9]+"))
    toks <- toks[nchar(toks) > 2 & !toks %in% c("diet")]
    if (!length(toks)) next
    arm_l <- tolower(cand$Treatment_Group)
    hit <- rep(FALSE, nrow(cand))
    for (tk in toks) hit <- hit | grepl(tk, arm_l, fixed = TRUE)
    # Fall back to an abbreviation match (standard -> sd, challenge -> cd)
    if (!any(hit)) {
      abbr <- c(standard = "sd", challenge = "cd", control = "ctrl")
      for (tk in toks) if (tk %in% names(abbr))
        hit <- hit | grepl(paste0("\\b", abbr[[tk]], "\\b"), arm_l)
    }
    if (any(hit)) {
      ghi_val[i] <- mean(cand[[GHI_TARGET]][hit], na.rm = TRUE)
      coarse_assigned[i] <- TRUE
    }
  }
  n_coarse <- sum(coarse_assigned)
  if (n_coarse > 0) {
    log_msg("  coarse-label fallback: %d run(s) assigned a study/diet-level mean GHI",
            n_coarse)
    for (s in unique(cw_matched$Study_ID[coarse_assigned])) {
      log_msg("    %-24s %d runs", s, sum(coarse_assigned & cw_matched$Study_ID == s))
    }
    log_msg("    (flagged in the output; repeat without them as a sensitivity check)")
  }
}

n_with <- sum(!is.na(ghi_val))
log_msg("  %d of %d runs received a %s value (%.0f%%)",
        n_with, length(ghi_val), GHI_TARGET, 100 * n_with / length(ghi_val))
if (n_with == 0) {
  log_msg("")
  log_msg("  No runs matched. The two tables likely disagree on how treatment")
  log_msg("  arms or segments are named. Examples of each side:")
  log_msg("    crosswalk keys: %s", paste(head(unique(run_key), 3), collapse = " | "))
  log_msg("    GHI keys:       %s", paste(head(unique(ghi_sub$.key), 3), collapse = " | "))
  stop("Join produced no matches - reconcile the labels before proceeding.")
}

# Report which studies contribute, since a genus association can only be
# estimated from studies that have both sequence data and a GHI value.
contributing <- unique(cw_matched$Study_ID[!is.na(ghi_val)])
log_msg("  %d studies contribute both sequence data and a GHI value",
        length(contributing))
if (length(contributing) < 10) {
  log_msg("  *** WARNING: fewer than 10 contributing studies. Leave-one-study-out")
  log_msg("  *** validation downstream will be unstable at this sample size.")
}

# When age is used as a covariate, a run also needs a non-missing age to
# enter the model - otherwise MaAsLin2 would drop it internally anyway, but
# silently, and the reported n would not match what was intended.
use_age <- AGE_COVARIATE && any(!is.na(age_val))
# A run needs a non-missing SEGMENT to enter when SEGMENT_COVARIATE is on -
# Segment is already required upstream for the join itself (see mk_key), so
# in practice this adds no further exclusions; kept explicit for clarity.
use_seg <- SEGMENT_COVARIATE && any(!is.na(cw_matched$Segment) & nzchar(cw_matched$Segment))
use_cls <- CLASS_COVARIATE && any(!is.na(cls_val) & nzchar(cls_val))

sel <- !is.na(ghi_val) &
       (!use_age | !is.na(age_val)) &
       (!use_cls | (!is.na(cls_val) & nzchar(cls_val))) &
       (!use_seg | (!is.na(cw_matched$Segment) & nzchar(cw_matched$Segment)))

if (use_age) log_msg("  age covariate requested: %d additional run(s) excluded for missing age",
                     sum(!is.na(ghi_val)) - sum(!is.na(ghi_val) & (!is.na(age_val))))
if (use_cls) log_msg("  class covariate requested: %d run(s) have a usable Challenge_Type",
                     sum(!is.na(cls_val) & nzchar(cls_val)))
if (use_seg) log_msg("  segment covariate requested: using Segment as an additional fixed effect")

clr_mat <- clr_mat[, sel, drop = FALSE]
meta <- data.frame(
  Study_ID  = cw_matched$Study_ID[sel],
  Segment   = cw_matched$Segment[sel],
  GHI       = ghi_val[sel],
  CoarseGHI = coarse_assigned[sel],
  row.names = colnames(clr_mat),
  stringsAsFactors = FALSE
)
if (use_age) meta$Age_d <- age_val[sel]
if (use_cls) meta$Challenge_Type <- factor(cls_val[sel])
if (use_seg) meta$Segment <- factor(meta$Segment)
log_msg("  final matrix: %d genera x %d runs", nrow(clr_mat), ncol(clr_mat))

# =============================================================================
# STEP 4 - MaAsLin2
#
# Study identity enters as a random effect, per the registered protocol, so
# that repeated contribution of a study does not inflate significance and
# residual between-study heterogeneity is absorbed.
#
# Normalisation and transformation are both set to NONE because the input is
# already CLR-transformed (03_batch_correction.R). MaAsLin2's defaults (TSS
# then LOG) assume raw or proportional abundances; applying them to CLR values
# would transform already-transformed data.
# =============================================================================
fixed_effects_used <- c("GHI",
                        if (use_age) "Age_d",
                        if (use_cls) "Challenge_Type",
                        if (use_seg) "Segment")
log_msg("[4/6] Running MaAsLin2 (fixed effect(s): %s; random effect: study)",
        paste(fixed_effects_used, collapse = " + "))
out_tag <- paste0(if (use_age) "_age" else "", if (use_cls) "_class" else "",
                  if (use_seg) "_seg" else "")
out_dir <- file.path(R_OUT_DIR, "gdi", paste0("maaslin2", out_tag))
if (out_tag == "") out_dir <- file.path(R_OUT_DIR, "gdi", "maaslin2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit <- Maaslin2::Maaslin2(
  input_data      = as.data.frame(t(clr_mat)),
  input_metadata  = meta,
  output          = out_dir,
  fixed_effects   = fixed_effects_used,
  random_effects  = c("Study_ID"),
  normalization   = "NONE",
  transform       = "NONE",
  min_prevalence  = 0,      # prevalence filtering already applied upstream
  min_abundance   = 0,
  max_significance = FDR_THRESHOLD,
  plot_heatmap    = FALSE,
  plot_scatter    = FALSE,
  standardize     = FALSE
)

# =============================================================================
# STEP 5 - EXTRACT THE GENERA THAT PASS
# =============================================================================
log_msg("[5/6] Extracting significant genera")
res_file <- file.path(out_dir, "all_results.tsv")
if (!file.exists(res_file)) stop("MaAsLin2 produced no results file.")
res <- fread(res_file, sep = "\t", header = TRUE)
res <- res[res$metadata == "GHI", ]
setorder(res, qval)

pass_loose  <- res[res$qval <= FDR_THRESHOLD, ]
pass_strict <- res[res$qval <= 0.05, ]

log_msg("  %d genera tested", nrow(res))
log_msg("  %d pass FDR <= %.2f  (these proceed to LASSO)", nrow(pass_loose), FDR_THRESHOLD)
log_msg("  %d pass FDR <= 0.05  (stricter comparison)", nrow(pass_strict))

if (nrow(pass_loose) == 0) {
  log_msg("")
  log_msg("  *** No genus is associated with %s at FDR <= %.2f.", GHI_TARGET, FDR_THRESHOLD)
  log_msg("  *** This is a substantive result, not a technical failure: it means")
  log_msg("  *** no single genus tracks the index once study is accounted for.")
  log_msg("  *** LASSO on an unfiltered genus set remains possible as a declared")
  log_msg("  *** deviation, but a null screen is itself worth reporting.")
} else {
  log_msg("")
  log_msg("  top genera by FDR:")
  n_show <- min(15, nrow(pass_loose))
  for (i in seq_len(n_show)) {
    log_msg("    %-42s coef=%+.4f  q=%.4g", pass_loose$feature[i],
            pass_loose$coef[i], pass_loose$qval[i])
  }
}

# =============================================================================
# STEP 6 - WRITE OUTPUTS FOR THE LASSO STAGE
# =============================================================================
log_msg("[6/6] Writing outputs")
sel_file <- file.path(R_OUT_DIR, "gdi", "maaslin2_selected_genera.txt")
writeLines(pass_loose$feature, sel_file)
fwrite(res, file.path(R_OUT_DIR, "gdi", "maaslin2_all_results.tsv"), sep = "\t")

# The joined matrix is saved so the LASSO stage uses exactly this data,
# rather than repeating the join and risking a different result.
saveRDS(list(clr = clr_mat, meta = meta, target = GHI_TARGET,
             n_coarse = sum(coarse_assigned[sel]),
             selected = pass_loose$feature),
        file.path(R_OUT_DIR, "gdi", "gdi_input.rds"))

summ <- c("=== MaAsLin2: genus associations with the GHI ===", "",
  sprintf("Target (pre-specified): %s", GHI_TARGET),
  sprintf("Runs analysed: %d", ncol(clr_mat)),
  sprintf("Genera tested: %d", nrow(res)),
  # `contributing` is computed before the analysed subset is taken, so it
  # counts the studies of the full CLR matrix, not of the runs actually
  # modelled. Both are reported, each labelled with what it counts.
  sprintf("Studies represented in the analysed runs: %d",
          length(unique(cw$Study_ID[match(colnames(clr_mat), cw$Run)]))),
  sprintf("Studies with sequence data and a GHI value, before filtering: %d",
          length(contributing)),
  "",
  sprintf("Pass FDR <= %.2f : %d  (proceed to LASSO)", FDR_THRESHOLD, nrow(pass_loose)),
  sprintf("Pass FDR <= 0.05 : %d", nrow(pass_strict)),
  "",
  "Model: per-genus linear mixed model, GHI as fixed effect, study identity",
  "as random effect. Input already CLR-transformed, so MaAsLin2 normalisation",
  "and transformation are both set to NONE.",
  "",
  "--- limitations ---",
  "1. Sequencing runs are joined to the GHI on Study_ID + Treatment_Group +",
  "   Segment. Several runs share one biomarker row, because papers report a",
  "   single group mean per arm. Every run in an arm therefore carries that",
  "   arm's GHI. This is outcome measurement error and it caps the achievable",
  "   correlation independently of the microbiome signal.",
  "",
  "2. Batch correction (MMUPHin) has already been applied, and study identity",
  "   also enters here as a random effect. This is deliberately conservative -",
  "   both are specified in the registered protocol - but it may remove some",
  "   genuine between-study biological variation along with technical variation.",
  "",
  "3. FDR 0.25 is MaAsLin2's default and is permissive for a screening step.",
  "   The 0.05 count is reported alongside so the effect of that choice is",
  "   visible rather than buried.")
writeLines(summ, file.path(R_OUT_DIR, "gdi", "maaslin2_summary.txt"))

log_msg("  selected genera -> %s", sel_file)
log_msg("  full results    -> %s", file.path(R_OUT_DIR, "gdi", "maaslin2_all_results.tsv"))
log_msg("  input for LASSO -> %s", file.path(R_OUT_DIR, "gdi", "gdi_input.rds"))
log_msg("=== DONE ===")
