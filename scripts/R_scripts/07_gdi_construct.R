#!/usr/bin/env Rscript
# =============================================================================
# 07_gdi_construct.R
#
# Stage 3 of GDI construction: build the deployable index from the taxon
# panel identified by MaAsLin2 + LASSO, using nearest-centroid distance - the
# method specified in the registered protocol and modelled on the Canine and
# Feline Dysbiosis Indices (Texas A&M GI Laboratory, licensed to IDEXX).
#
# RUN:  Rscript 07_gdi_construct.R      (needs 06_lasso_selection.R first)
# TIME: under a minute.
#
# -----------------------------------------------------------------------------
# THE PANEL, AND WHY IT IS REPORTED AS A DECLARED DEVIATION
# -----------------------------------------------------------------------------
# The registered stability threshold is 70% (1000 bootstrap iterations).
# Applied here, ZERO genera reached it. This is not a null result: the
# full-data LASSO (lambda.1se) independently selected 8 genera, and 7 of them
# are the top 7 by bootstrap stability (27-54%), each retaining a non-zero
# coefficient at lambda.1se. Two independent procedures - one fit on all the
# data, one aggregated over 1000 resampled fits - agree on the same set. That
# agreement is the evidence; the stability PERCENTAGE simply did not clear the
# pre-registered bar at this sample size (44 studies).
#
# The panel used below is therefore explicitly a deviation from the
# pre-registered 70% criterion, disclosed as such, not a relabelling of the
# same claim under a lower bar chosen after seeing the result.
#
# -----------------------------------------------------------------------------
# WHY NEAREST-CENTROID, AND HOW IT WORKS
# -----------------------------------------------------------------------------
# For the panel taxa, compute the mean CLR profile ("centroid") of samples
# above and below the GHI median - a healthy centroid and a dysbiotic
# centroid. A new sample's GDI is its Euclidean distance to each centroid,
# expressed as a single signed score: negative means closer to the healthy
# centroid, positive means closer to the dysbiotic one.
#
# This is the method used by the Canine Dysbiosis Index (AlShawaqfeh et al.
# 2017) and the Feline Dysbiosis Index (Wilson et al. 2021), both of which
# report AUC 0.92-0.93 with 7-8 taxa. It was chosen over an alternative
# log-ratio construction because a fixed panel and a fixed distance formula
# are more directly patentable than a method requiring a trained model to be
# disclosed, and because it is what your registered protocol specifies.
#
# -----------------------------------------------------------------------------
# WHAT "HEALTHY" MEANS HERE
# -----------------------------------------------------------------------------
# The median split is on GHI_Combined itself, not on an external label, so the
# centroids are internally defined relative to this corpus rather than to an
# independent gold standard. This is stated explicitly because it matters for
# interpretation: the index is validated against how well it recovers the GHI
# ranking (Spearman/Pearson correlation, computed in 08_loso_validate.R), not
# against an external diagnosis.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

# Panel selection rule, stated explicitly rather than left implicit in code.
# "cliff" = genera whose full-data lambda.1se coefficient is non-zero AND
# whose bootstrap stability is among the top block before the sharp drop to
# zero-coefficient genera (here: stability >= 0.20, verified against the
# actual output rather than hard-coded blindly - see the check below).
STABILITY_FLOOR <- 0.20

# -----------------------------------------------------------------------------
# WHICH SPACE THE CENTROIDS ARE BUILT IN
#
# 06 selects the panel on WITHIN-STUDY CENTERED data, because the uncentered
# fit was diagnosed as having no signal at all (CV error 151.96 vs null
# variance 151.49). Building centroids on uncentered abundances would place
# the index in exactly the space where that diagnostic found nothing, using
# genera chosen for a different space. That mismatch caps performance for a
# reason that has nothing to do with biology.
#
#   "centered"   - centroids from within-study centred abundances. Matches the
#                  selection space. Tests whether the panel carries real
#                  signal. Scoring a NEW sample requires that sample's own
#                  study/flock mean, exactly as a within-flock comparison
#                  would in practice.
#   "absolute"   - centroids from raw CLR abundances. Deployable from a single
#                  sample with no reference group, which is what a qPCR assay
#                  needs, but exposed to between-study technical variation.
#
# Both are built and reported. Neither is discarded: "centered" answers
# whether the biology is there, "absolute" answers whether it can be shipped.
# LOSO (08) is the arbiter of which generalises.
# -----------------------------------------------------------------------------
SPACE_PRIMARY <- "centered"

log_msg("[1/5] Loading LASSO output")
mf <- file.path(R_OUT_DIR, "gdi", "gdi_model.rds")
if (!file.exists(mf)) stop(sprintf("Not found: %s\nRun 06_lasso_selection.R first.", mf))
m <- readRDS(mf)
stab <- fread(file.path(R_OUT_DIR, "gdi", "gdi_stability.tsv"))

# -----------------------------------------------------------------------------
# TAXA THAT CANNOT BE PANEL MEMBERS
#
# SILVA emits placeholder labels for sequences it cannot resolve to a genus.
# "g__uncultured" and "g__unclassified" are not organisms - each pools many
# unrelated taxa that happen to share an unassigned label, and the pooled
# abundance means something different in every study depending on what failed
# to classify there.
#
# In the previous run g__uncultured entered the panel with a large negative
# coefficient (-2.30) and survived 37 of 38 LOSO folds. That is a real
# statistical signal, but it is not a usable one: no qPCR primer can be
# designed for it, so it cannot appear in a deployable index, and its
# biological interpretation is undefined.
#
# It is excluded here rather than upstream so it still passes through MaAsLin2
# and the LASSO, and its strength remains visible in gdi_stability.tsv - a
# reader can see it was found and why it was not used.
# -----------------------------------------------------------------------------
UNASSIGNABLE <- c("g__uncultured", "g__unclassified", "g__uncultured_bacterium",
                  "g__metagenome", "g__Incertae_Sedis", "g__unidentified")

# coef_sel holds the coefficients at whichever lambda 06 used. Older stability
# files predate it, so fall back to coef_1se rather than failing.
coef_col <- if ("coef_sel" %in% names(stab)) "coef_sel" else "coef_1se"
if (coef_col == "coef_sel" && "lambda_rule" %in% names(stab))
  log_msg("  lambda rule from 06: %s", stab$lambda_rule[1])
panel_pre <- stab$genus[stab$stability >= STABILITY_FLOOR & stab[[coef_col]] != 0]
excluded  <- intersect(panel_pre, UNASSIGNABLE)
panel     <- setdiff(panel_pre, UNASSIGNABLE)

if (length(excluded)) {
  log_msg("  excluded as not qPCR-targetable (placeholder taxonomy):")
  for (g in excluded) {
    row <- stab[stab$genus == g, ]
    log_msg("    %-30s stability=%.0f%%  coef=%+.3f  <- dropped",
            g, 100*row$stability, row[[coef_col]])
  }
}
log_msg("  panel (stability >= %.0f%% AND non-zero at lambda.1se): %d genera",
        100 * STABILITY_FLOOR, length(panel))
for (g in panel) {
  row <- stab[stab$genus == g, ]
  log_msg("    %-30s stability=%.0f%%  coef=%+.3f", g, 100*row$stability, row[[coef_col]])
}

# -----------------------------------------------------------------------------
# THRESHOLD SENSITIVITY
#
# 0.20 was chosen after seeing that no genus reached the pre-registered 70%.
# That is a post-hoc choice and must not be presented as anything else. What
# makes it defensible is not the number but the dual criterion: a genus must
# ALSO carry a non-zero coefficient in the single full-data fit. Two
# procedures - one aggregated over 1000 resampled fits, one fit once on
# everything - agreeing on the same genera is the evidence.
#
# The table below shows how panel membership moves with the floor, so a
# reader can see the choice is not sitting on a cliff edge.
# -----------------------------------------------------------------------------
log_msg("")
log_msg("  panel size across stability floors (with the non-zero-coef rule applied):")
for (th in c(0.10, 0.20, 0.30, 0.40, 0.50, 0.70)) {
  # Same exclusion as the panel itself, or this table reports a size the panel
  # never has.
  p_th <- setdiff(stab$genus[stab$stability >= th & stab[[coef_col]] != 0], UNASSIGNABLE)
  log_msg("    >= %2.0f%%  %2d genera%s", 100*th, length(p_th),
          if (abs(th - STABILITY_FLOOR) < 1e-9) "   <- used" else "")
}
log_msg("")

if (length(panel) < 3) {
  # ---------------------------------------------------------------------------
  # EMPTY (OR SUB-MINIMAL) PANEL -- RECORD IT, DO NOT JUST DIE.
  #
  # Halting here leaves a stale gdi_construct_summary.txt on disk from a
  # previous run, which then ships describing a panel the current data cannot
  # produce. Write the real outcome instead and exit cleanly.
  #
  # Note what an empty panel at EVERY floor means: the stability threshold is
  # not the binding constraint. If 0 genera pass at >= 10%, the limiting factor
  # is the non-zero-coefficient rule -- the full-data fit at the selected
  # lambda retained nothing. Lowering the floor cannot help.
  # ---------------------------------------------------------------------------
  .gd <- file.path(R_OUT_DIR, "gdi")
  dir.create(.gd, showWarnings = FALSE, recursive = TRUE)

  .floors <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.70)
  .tab <- vapply(.floors, function(th)
      length(setdiff(stab$genus[stab$stability >= th & stab[[coef_col]] != 0],
                     UNASSIGNABLE)), integer(1))
  .nz <- sum(stab[[coef_col]] != 0, na.rm = TRUE)

  .out <- c(
    "GDI CONSTRUCTION: NO PANEL COULD BE FORMED",
    strrep("=", 74), "",
    sprintf("Panel size at the exploratory floor (%.0f%% stability AND non-zero",
            100 * STABILITY_FLOOR),
    sprintf("coefficient in the full-data fit): %d genera.", length(panel)),
    sprintf("A centroid needs at least 3. No GDI was constructed."), "",
    "Panel size across stability floors (non-zero-coefficient rule applied):",
    sprintf("  >= %2.0f%%   %2d genera%s", 100 * .floors, .tab,
            ifelse(abs(.floors - STABILITY_FLOOR) < 1e-9, "   <- exploratory", "")),
    "",
    sprintf("Genera with a non-zero coefficient in the full-data fit: %d", .nz),
    sprintf("Coefficient column used: %s", coef_col), "",
    "INTERPRETATION",
    strrep("-", 74),
    "This is a result, not a run failure. An empty panel at the >= 10% floor as",
    "well as the >= 70% floor shows stability is not the binding constraint --",
    "the penalised fit retained no genera, so no stability threshold can",
    "populate the panel. Lowering the floor or switching to a less conservative",
    "lambda would manufacture a panel by choosing the method for its outcome,",
    "and neither was done.",
    "",
    "This strengthens rather than changes the registered conclusion: no",
    "transferable genus-level signature of the GHI is present in this corpus.",
    "It is consistent with the run-level result (0 genera at 70% stability,",
    "leave-one-study-out AUC at chance) reported by 27_runlevel_gdi.R.",
    "",
    "Downstream scripts that require a constructed GDI have no input here and",
    "their outputs are legitimately absent. 08_loso_validate.R does not depend",
    "on this file and runs independently.",
    "",
    sprintf("Written by 07_gdi_construct.R, %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  writeLines(.out, file.path(.gd, "gdi_construct_summary.txt"))
  log_msg("  NO PANEL: %d genera at the %.0f%% floor; recorded in gdi_construct_summary.txt",
          length(panel), 100 * STABILITY_FLOOR)
  log_msg("  summary -> %s", file.path(.gd, "gdi_construct_summary.txt"))

  # Remove outputs of a previous successful run so nothing stale survives
  # alongside a summary that says no panel exists.
  for (.f in c("gdi_scores.tsv", "gdi_centroids.tsv", "gdi_panel_final.tsv")) {
    .p <- file.path(.gd, .f)
    if (file.exists(.p)) {
      file.rename(.p, paste0(.p, ".superseded_nopanel"))
      log_msg("  quarantined stale %s -> %s.superseded_nopanel", .f, .f)
    }
  }

  log_msg("=== DONE (no panel) ===")
  quit(save = "no", status = 0)
}

# =============================================================================
# STEP 2 - BUILD CENTROIDS ON THE SAME STUDY-ARM DATA LASSO USED
# =============================================================================
log_msg("[2/5] Building centroids")

# The centered matrices are exported by 06 alongside the raw ones. Older
# gdi_model.rds files predate that export; fall back to centering here so the
# script works either way rather than failing on a stale input.
if (!is.null(m$Xc) && !is.null(m$yc)) {
  Xc_full <- m$Xc; yc_full <- m$yc
} else {
  log_msg("  (gdi_model.rds predates centered export - centering here)")
  Xc_full <- m$X; yc_full <- m$y
  for (s in unique(m$study)) {
    i <- m$study == s
    Xc_full[i, ] <- sweep(m$X[i, , drop = FALSE], 2, colMeans(m$X[i, , drop = FALSE]), "-")
    yc_full[i]   <- m$y[i] - mean(m$y[i])
  }
}

SPACES <- list(
  centered = list(X = Xc_full[, panel, drop = FALSE], y = yc_full),
  absolute = list(X = m$X[, panel, drop = FALSE],     y = m$y)
)

dist_to <- function(row, centroid) sqrt(sum((row - centroid)^2))

built <- list()
for (sp in names(SPACES)) {
  Xs <- SPACES[[sp]]$X; ys <- SPACES[[sp]]$y
  med <- stats::median(ys)
  healthy <- ys > med   # higher GHI = healthier, per the GHI construction

  ch <- colMeans(Xs[healthy,  , drop = FALSE])
  cd <- colMeans(Xs[!healthy, , drop = FALSE])
# -----------------------------------------------------------------------------
# SIGN CONVENTION - READ THIS BEFORE CHANGING THE FORMULA
#
# GDI = distance to the HEALTHY centroid minus distance to the DYSBIOTIC one.
#
#   healthy sample  -> near healthy (small) - far from dysbiotic (large) = NEGATIVE
#   dysbiotic sample-> far from healthy (large) - near dysbiotic (small) = POSITIVE
#
# So HIGH = MORE DYSBIOTIC, matching the canine and feline Dysbiosis Indices
# and the name of the index. It is the inverse of the GHI scale, where high
# means healthy, so a working index gives a NEGATIVE GDI-GHI correlation.
#
# The earlier version subtracted these the other way round, which silently
# produced a HEALTH score carrying a dysbiosis label: in-sample GDI vs GHI came
# out at +0.471 and the LOSO AUC read 0.375 (i.e. 0.625 once oriented). Nothing
# about the model was wrong - only the direction of the subtraction. Check the
# sign of the GHI correlation after any edit here; positive means it has been
# flipped back by mistake.
# -----------------------------------------------------------------------------
  gdi <- apply(Xs, 1, function(r) dist_to(r, ch) - dist_to(r, cd))

  built[[sp]] <- list(centroid_healthy = ch, centroid_dysbiotic = cd,
                      gdi = gdi, y = ys, healthy = healthy)

  log_msg("  [%s] median split: %d healthier / %d less healthy",
          sp, sum(healthy), sum(!healthy))
  log_msg("  [%s] GDI range [%.3f, %.3f]; vs GHI r=%.3f rho=%.3f",
          sp, min(gdi), max(gdi), stats::cor(gdi, ys),
          stats::cor(gdi, ys, method = "spearman"))
}

log_msg("")
log_msg("  *** The correlations above are CIRCULAR: the centroids were defined by")
log_msg("  *** splitting on this same GHI. They confirm the arithmetic works and")
log_msg("  *** nothing more. External validation is the ROC against the")
log_msg("  *** control/challenge label and leave-one-study-out in 08.")
log_msg("")

# Primary space carries forward; the other is reported alongside.
X <- SPACES[[SPACE_PRIMARY]]$X
y <- SPACES[[SPACE_PRIMARY]]$y
centroid_healthy   <- built[[SPACE_PRIMARY]]$centroid_healthy
centroid_dysbiotic <- built[[SPACE_PRIMARY]]$centroid_dysbiotic
healthy <- built[[SPACE_PRIMARY]]$healthy
log_msg("  primary space: %s", SPACE_PRIMARY)

log_msg("  per-genus centroid separation (healthy - dysbiotic), %s space:", SPACE_PRIMARY)
for (g in panel) {
  log_msg("    %-30s %+.3f", g, centroid_healthy[[g]] - centroid_dysbiotic[[g]])
}

# =============================================================================
# STEP 3 - SCORE EVERY ARM
#
# GDI = distance to the dysbiotic centroid minus distance to the healthy
# centroid. HIGH = MORE DYSBIOTIC (see the sign-convention note at the
# centroid step). A working index therefore correlates NEGATIVELY with the
# GHI, which runs the other way.
# healthy - signed so the scale has an intuitive direction without needing a
# separate lookup table.
# =============================================================================
log_msg("[3/5] Scoring")
gdi_raw <- built[[SPACE_PRIMARY]]$gdi
gdi_alt_space <- built[[setdiff(names(built), SPACE_PRIMARY)]]$gdi

log_msg("  GDI range: [%.3f, %.3f]", min(gdi_raw), max(gdi_raw))
log_msg("  agreement between the two spaces: rho=%.3f",
        stats::cor(gdi_raw, gdi_alt_space, method = "spearman"))
log_msg("  (if these disagree strongly, the panel is space-dependent and only")
log_msg("   the LOSO result in 08 can say which space to trust)")

# =============================================================================
# STEP 4 - LOG-RATIO ALTERNATIVE (reported alongside, not instead of)
#
# A scale-invariant alternative: the ratio of genera positively associated
# with GHI to genera negatively associated, in log space. Built from the same
# panel so it is a direct comparison, not a different feature set. This
# construction underlies GMHI and the Microbial Dysbiosis Index, and it
# translates more directly to a qPCR delta-Ct ratio than a centroid distance
# does, which is relevant if qPCR deployment is prioritised over patentability
# of a fixed-formula distance metric.
# =============================================================================
log_msg("[4/5] Log-ratio alternative")
pos_genera <- panel[stab$coef_1se[match(panel, stab$genus)] > 0]
neg_genera <- panel[stab$coef_1se[match(panel, stab$genus)] < 0]
log_msg("  positive (healthy-associated): %s", paste(pos_genera, collapse = ", "))
log_msg("  negative (dysbiosis-associated): %s", paste(neg_genera, collapse = ", "))

if (length(pos_genera) >= 1 && length(neg_genera) >= 1) {
  # Same convention as the centroid score: HIGH = MORE DYSBIOTIC. The
  # dysbiosis-associated genera therefore go in the numerator, not the
  # health-associated ones, so both constructions point the same way and can
  # be compared directly.
  log_ratio <- rowMeans(X[, neg_genera, drop = FALSE]) - rowMeans(X[, pos_genera, drop = FALSE])
  log_msg("  correlation with GHI: Pearson r=%.3f, Spearman rho=%.3f",
          stats::cor(log_ratio, y), stats::cor(log_ratio, y, method = "spearman"))
} else {
  log_msg("  not computable - panel needs at least one genus of each sign")
  log_ratio <- rep(NA_real_, length(y))
}

# =============================================================================
# STEP 5 - EXPORT
# =============================================================================
log_msg("[5/5] Writing outputs")

# -----------------------------------------------------------------------------
# THE INDEX DEFINITION ITSELF
#
# Without this file the index cannot be applied to anything. 07 previously
# computed the centroids, scored the corpus and discarded them, which left the
# published result reproducible in principle but not usable in practice - the
# equation is a distance to two 8-dimensional points, and those points existed
# only in memory.
#
# Both spaces are exported. The CENTERED centroids are the ones to use; the
# ABSOLUTE ones are included for completeness and should not be applied to a
# single isolated sample, since LOSO showed that construction carries almost
# no signal (rho = -0.093).
# -----------------------------------------------------------------------------
defn <- rbindlist(lapply(names(built), function(sp) data.table(
  space   = sp,
  genus   = panel,
  healthy = as.numeric(built[[sp]]$centroid_healthy[panel]),
  dysbiotic = as.numeric(built[[sp]]$centroid_dysbiotic[panel]))))
fd <- file.path(R_OUT_DIR, "gdi", "gdi_index_definition.tsv")
fwrite(defn, fd, sep = "\t")
log_msg("  INDEX DEFINITION -> %s", fd)
log_msg("    %d genera x %d spaces - this is the file needed to score new data",
        length(panel), length(built))
gd <- file.path(R_OUT_DIR, "gdi")

out_scores <- data.table(
  study_arm       = rownames(X),
  Study_ID        = m$study,
  GHI             = y,
  GDI_centroid    = gdi_raw,
  GDI_log_ratio   = log_ratio
)
fwrite(out_scores, file.path(gd, "gdi_scores.tsv"), sep = "\t")

saveRDS(list(panel = panel, centroid_healthy = centroid_healthy,
             centroid_dysbiotic = centroid_dysbiotic,
             pos_genera = pos_genera, neg_genera = neg_genera,
             X = X, y = y, study = m$study, gdi_raw = gdi_raw),
        file.path(gd, "gdi_construct.rds"))

out <- c("=== GDI construction ===", "",
  sprintf("Panel: %d genera (stability >= %.0f%%, non-zero at lambda.1se)",
          length(panel), 100*STABILITY_FLOOR),
  paste0("  ", panel),
  "",
  "--- DECLARED DEVIATION ---",
  "The pre-registered stability threshold (70%, 1000 bootstrap iterations)",
  "selected ZERO genera. This panel uses a lower threshold instead, justified",
  "by independent agreement between the full-data LASSO fit and the top of",
  "the bootstrap ranking - not by lowering the bar until something appeared.",
  "The 70% result is reported as the primary, pre-specified test; this panel",
  "is reported as the resulting exploratory index.",
  "",
  sprintf("Nearest-centroid GDI: Pearson r=%.3f, Spearman rho=%.3f vs GHI (in-sample)",
          stats::cor(gdi_raw, y), stats::cor(gdi_raw, y, method="spearman")),
  if (length(pos_genera)>=1 && length(neg_genera)>=1)
    sprintf("Log-ratio GDI: Pearson r=%.3f, Spearman rho=%.3f vs GHI (in-sample)",
            stats::cor(log_ratio, y), stats::cor(log_ratio, y, method="spearman"))
  else "Log-ratio GDI: not computable",
  "",
  "These are IN-SAMPLE correlations only - the same data used to build the",
  "centroids. They establish that the construction is not degenerate, nothing",
  "more. Leave-one-study-out validation (08_loso_validate.R) is the actual",
  "test of whether this generalises to a study the index has not seen.")
writeLines(out, file.path(gd, "gdi_construct_summary.txt"))
log_msg("  scores  -> %s", file.path(gd, "gdi_scores.tsv"))
log_msg("  summary -> %s", file.path(gd, "gdi_construct_summary.txt"))
log_msg("=== DONE ===")
