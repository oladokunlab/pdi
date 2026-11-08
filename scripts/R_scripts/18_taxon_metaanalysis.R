#!/usr/bin/env Rscript
# =============================================================================
# 18_taxon_metaanalysis.R
#
# Per-taxon random-effects meta-analysis: how consistent is each candidate
# genus's relationship with the GHI, study by study?
#
# RUN:  Rscript 18_taxon_metaanalysis.R     (needs 05_maaslin2.R to have run)
# TIME: under a minute.
#
# -----------------------------------------------------------------------------
# WHY THIS IS A DIFFERENT QUESTION FROM STABILITY SELECTION
# -----------------------------------------------------------------------------
# Stability selection (06_lasso_selection.R) asks: does this genus keep
# getting SELECTED when the model is refit on resampled subsets of studies.
# That is a question about the selection procedure's behaviour.
#
# This script asks a different question directly: for a given genus, is its
# correlation with the GHI the same sign and a similar size in every study, or
# does it point one way in some studies and the opposite way in others. A
# genus can survive stability selection while being wildly inconsistent
# study-to-study, if its AVERAGE relationship is non-zero even though
# individual studies disagree. Conversely a genus with a small but highly
# consistent effect could be dropped by stability selection (which rewards
# strength of the fitted coefficient) while being a more reliable biological
# signal than a stronger, noisier one.
#
# I-squared is what this script adds that stability selection cannot report:
# the fraction of variation in the taxon's effect that is due to real
# between-study differences rather than sampling noise. I2 < 50% is treated as
# a reasonably consistent association; I2 > 75% means the studies substantially
# disagree with each other about the direction or size of the effect, and that
# genus should not be trusted regardless of its average.
#
# -----------------------------------------------------------------------------
# UNIT OF ANALYSIS: ARMS, NOT RUNS
# -----------------------------------------------------------------------------
# The source template for this step (the pipeline guide) computes per-study
# correlations on individual sequencing runs. That repeats the exact
# pseudo-replication problem this project diagnosed and fixed early on: one
# study can contribute 700+ runs from a handful of distinct GHI values, which
# would make its correlation estimate look far more precise than it is and
# let it dominate the pooled result. Every other stage of this pipeline (LASSO,
# LOSO, MaAsLin2 itself) treats the treatment ARM as the unit; this script
# does the same, for consistency and for the same reason.
#
# -----------------------------------------------------------------------------
# THE MINIMUM-N REQUIREMENT, AND WHAT IT COSTS
# -----------------------------------------------------------------------------
# A within-study Spearman correlation needs a minimum number of arms to mean
# anything. Below MIN_ARMS_PER_STUDY, a "correlation" is close to meaningless -
# with 3 points it can only take a handful of possible values and is often
# exactly +1 or -1. Studies below that floor do not contribute to a taxon's
# pooled estimate. Given this corpus averages roughly 3-4 arms per study, this
# is a real and disclosed constraint: many studies will not qualify for many
# genera, and the number of studies contributing to each taxon's pooled
# estimate is reported alongside it rather than assumed.
#
# -----------------------------------------------------------------------------
# CANDIDATE GENERA
# -----------------------------------------------------------------------------
# Restricted to the genera that passed MaAsLin2 screening (05_maaslin2.R),
# matching the source template's design of using this as a second-stage filter
# on a pre-screened candidate list rather than testing all genera and
# incurring a much harsher multiple-testing correction. If a sparse-PLS
# candidate list exists later, it should be unioned in here.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("metafor", quietly = TRUE))
  stop("Package 'metafor' required.\n  install.packages(\"metafor\", lib = Sys.getenv(\"R_LIBS_USER\"), repos = \"https://cloud.r-project.org\")")

MIN_ARMS_PER_STUDY <- 5     # per the source template's implicit floor
MIN_STUDIES_PER_TAXON <- 3  # fewer than this and rma() is not meaningfully fit
I2_THRESHOLD <- 50          # consistent below this, per the template
P_THRESHOLD  <- 0.05

log_msg("[1/5] Loading MaAsLin2 candidates and the GDI model data")
cand_file <- file.path(R_OUT_DIR, "gdi", "maaslin2_selected_genera.txt")
if (!file.exists(cand_file)) stop("Not found: ", cand_file, "\n  Run 05_maaslin2.R first.")
candidates <- readLines(cand_file)
log_msg("  %d candidate genera from MaAsLin2 (FDR <= 0.25)", length(candidates))

mf <- file.path(R_OUT_DIR, "gdi", "gdi_model.rds")
if (!file.exists(mf)) stop("Not found: ", mf, "\n  Run 06_lasso_selection.R first.")
m <- readRDS(mf)
X <- m$X; y <- m$y; study <- m$study
# BUGFIX (systematic audit, 13 Sep 2026): plain intersect() against
# candidates read from maaslin2_selected_genera.txt silently drops any
# genus MaAsLin2 sanitised differently from X's column names (same root
# cause documented and fixed in 06_lasso_selection.R - MaAsLin2 replaces
# non-alphanumeric characters with periods; the CLR-derived matrix keeps
# the original dashes and brackets). Matching on a canonicalised form
# prevents the same silent drop here.
canon <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
match_idx <- match(canon(candidates), canon(colnames(X)))
candidates <- colnames(X)[na.omit(match_idx)]
log_msg("  %d of those present in the arm-level matrix", length(candidates))
log_msg("  %d arms, %d studies total", nrow(X), length(unique(study)))

# =============================================================================
# 2 - PER-STUDY SPEARMAN, PER TAXON
# =============================================================================
log_msg("[2/5] Computing per-study Spearman correlations")

per_taxon <- function(g) {
  abd <- X[, g]
  rows <- lapply(unique(study), function(s) {
    idx <- study == s
    if (sum(idx) < MIN_ARMS_PER_STUDY) return(NULL)
    a <- abd[idx]; h <- y[idx]
    ok <- is.finite(a) & is.finite(h)
    if (sum(ok) < MIN_ARMS_PER_STUDY) return(NULL)
    r <- suppressWarnings(stats::cor(a[ok], h[ok], method = "spearman"))
    if (!is.finite(r)) return(NULL)
    r <- max(min(r, 0.999), -0.999)   # atanh is undefined at exactly +-1
    n <- sum(ok)
    data.frame(study = s, r = r, n = n, ri = atanh(r), sei = 1 / sqrt(n - 3))
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) < MIN_STUDIES_PER_TAXON) return(NULL)
  do.call(rbind, rows)
}

results <- list()
for (g in candidates) {
  sc <- per_taxon(g)
  if (is.null(sc)) { results[[g]] <- NULL; next }
  fit <- tryCatch(
    metafor::rma(yi = sc$ri, sei = sc$sei, method = "REML"),
    error = function(e) NULL)
  if (is.null(fit)) { results[[g]] <- NULL; next }
  results[[g]] <- list(
    taxon = g, n_studies = nrow(sc), n_arms = sum(sc$n),
    pooled_r = tanh(as.numeric(fit$beta)),
    ci_lower = tanh(fit$ci.lb), ci_upper = tanh(fit$ci.ub),
    pval = fit$pval, I2 = fit$I2, studies = sc$study)
}
n_evaluable <- sum(!vapply(results, is.null, logical(1)))
log_msg("  %d of %d candidates had enough studies to fit (>= %d studies, >= %d arms each)",
        n_evaluable, length(candidates), MIN_STUDIES_PER_TAXON, MIN_ARMS_PER_STUDY)

# =============================================================================
# 3 - REPORT AND FILTER
# =============================================================================
log_msg("[3/5] Pooled effects, sorted by |pooled r|")
tab <- Filter(Negate(is.null), results)
if (!length(tab)) {
  log_msg("  No candidate genus had enough per-study data to evaluate.")
  log_msg("  This is itself informative: it means the corpus cannot currently")
  log_msg("  support this consistency check, independent of whether any true")
  log_msg("  association exists.")
} else {
  ord <- order(-abs(vapply(tab, function(z) z$pooled_r, numeric(1))))
  tab <- tab[ord]
  log_msg("")
  log_msg("  %-30s %8s %8s %10s %12s %8s %8s", "taxon", "studies", "arms",
          "pooled r", "95% CI", "p", "I2")
  for (z in tab)
    log_msg("  %-30s %8d %8d %+10.3f [%+.2f,%+.2f] %8.4f %7.0f%%%s",
            gsub("^g__", "", z$taxon), z$n_studies, z$n_arms, z$pooled_r,
            z$ci_lower, z$ci_upper, z$pval, z$I2,
            if (z$I2 < I2_THRESHOLD && z$pval < P_THRESHOLD) "  * stable" else "")
}

stable <- Filter(function(z) z$I2 < I2_THRESHOLD && z$pval < P_THRESHOLD, tab)
log_msg("")
log_msg("[4/5] Stable candidates (I2 < %d%%, p < %.2f): %d of %d evaluable",
        I2_THRESHOLD, P_THRESHOLD, length(stable), length(tab))
if (length(stable)) {
  for (z in stable)
    log_msg("    %-30s pooled r=%+.3f  I2=%.0f%%  (%d studies)",
            gsub("^g__", "", z$taxon), z$pooled_r, z$I2, z$n_studies)
} else {
  log_msg("    None. Either no candidate is consistent across studies at this")
  log_msg("    threshold, or too few studies had enough arms to evaluate it.")
}

# =============================================================================
# 4 - COMPARE AGAINST STABILITY SELECTION
# =============================================================================
log_msg("[5/5] Cross-check against LASSO stability selection")
stab_file <- file.path(R_OUT_DIR, "gdi", "gdi_stability.tsv")
if (file.exists(stab_file) && length(tab)) {
  stabtab <- fread(stab_file)
  log_msg("  %-30s %10s %10s", "taxon", "meta-I2", "LASSO stab")
  for (z in tab) {
    srow <- stabtab[genus == z$taxon]
    stabval <- if (nrow(srow)) sprintf("%.0f%%", 100 * srow$stability[1]) else "n/a"
    log_msg("  %-30s %9.0f%% %10s%s", gsub("^g__", "", z$taxon), z$I2, stabval,
            if (z$I2 < I2_THRESHOLD && nrow(srow) && srow$stability[1] >= 0.20)
              "  agrees with both criteria" else "")
  }
  log_msg("")
  log_msg("  A genus agreeing on both - low I2 here AND high stability there -")
  log_msg("  is supported by two independent kinds of evidence: consistency")
  log_msg("  across studies, and robustness to which studies are included.")
} else {
  log_msg("  gdi_stability.tsv not found or no evaluable taxa - skipping cross-check")
}

# =============================================================================
# WRITE
# =============================================================================
out <- c("PER-TAXON RANDOM-EFFECTS META-ANALYSIS", strrep("=", 74), "",
  sprintf("Candidates from MaAsLin2 (FDR <= 0.25): %d", length(candidates)),
  sprintf("Evaluable (>= %d studies with >= %d arms each): %d",
          MIN_STUDIES_PER_TAXON, MIN_ARMS_PER_STUDY, length(tab)),
  sprintf("Stable (I2 < %d%%, p < %.2f): %d", I2_THRESHOLD, P_THRESHOLD, length(stable)),
  "",
  sprintf("  %-30s %8s %8s %10s %8s %8s", "taxon", "studies", "arms", "pooled r", "p", "I2"))
for (z in tab)
  out <- c(out, sprintf("  %-30s %8d %8d %+10.3f %8.4f %7.0f%%",
                        gsub("^g__", "", z$taxon), z$n_studies, z$n_arms,
                        z$pooled_r, z$pval, z$I2))
out <- c(out, "",
  "This tests consistency of each taxon's GHI association across studies,",
  "independent of and complementary to LASSO stability selection. Arms are the",
  "unit of analysis throughout, matching every other stage of this pipeline.")
writeLines(out, file.path(R_OUT_DIR, "gdi", "taxon_metaanalysis.txt"))
fwrite(rbindlist(lapply(tab, function(z) data.table(
         taxon = z$taxon, n_studies = z$n_studies, n_arms = z$n_arms,
         pooled_r = z$pooled_r, ci_lower = z$ci_lower, ci_upper = z$ci_upper,
         pval = z$pval, I2 = z$I2))),
       file.path(R_OUT_DIR, "gdi", "taxon_metaanalysis.tsv"), sep = "\t")
log_msg("")
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "taxon_metaanalysis.txt"))
log_msg("=== DONE ===")
