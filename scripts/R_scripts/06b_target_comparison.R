#!/usr/bin/env Rscript
# =============================================================================
# 06b_target_comparison.R
#
# Does the genus panel depend on WHICH GHI version it is regressed against?
#
# RUN:  Rscript 06b_target_comparison.R
# TIME: 3x the 05+06 runtime. Expect 15-45 minutes.
#
# -----------------------------------------------------------------------------
# WHY THIS EXISTS
# -----------------------------------------------------------------------------
# 04f found the Morphology and Biomarker arms of the GHI are uncorrelated at
# sample level (r = -0.007, n = 136), and 04g found the domains form three
# groups that do not track each other: Morphology; Barrier/Claudin1/SCFA/
# Inflammation; and the two Immunoglobulin domains.
#
# That does NOT invalidate the GHI. The GHI is a FORMATIVE index - gut health
# is DEFINED as the combination of these facets, in the same sense that the UN
# Human Development Index defines development as income plus longevity plus
# schooling. Formative indices do not require their components to correlate,
# and cannot be validated by internal consistency. They are validated
# EXTERNALLY, which 04e does against known intervention direction.
#
# What the orthogonality DOES raise is a concrete modelling question. A sample
# scored from morphology alone and a sample scored from cytokines alone both
# enter the LASSO as "y", but they were measured on different instruments. If
# those instruments carry different microbial signals, the model is being
# asked to predict a target whose meaning shifts between rows - which weakens
# selection and may be part of why stability topped out at 54% rather than the
# registered 70%.
#
# This script answers that empirically rather than by argument: run the same
# pipeline three times, once per GHI version, and compare the panels.
#
# -----------------------------------------------------------------------------
# THIS IS A COMPARISON, NOT A COMPETITION
# -----------------------------------------------------------------------------
# All three targets are reported, whatever they show. Running three models and
# keeping whichever gives the most flattering panel is selection on the
# outcome - the same error that was correctly avoided when GHI_Combined was
# kept as primary despite losing both validation tests in 04d/04e.
#
# GHI_Combined remains the pre-specified primary. Nothing here changes that.
# A change of primary target is a protocol decision for the PI and
# statistician, and needs a PROSPERO amendment.
#
# -----------------------------------------------------------------------------
# HOW TO READ THE RESULT
# -----------------------------------------------------------------------------
#   PANELS OVERLAP SUBSTANTIALLY
#     The facets converge at the microbiome level even though they do not at
#     the biomarker level. That is a genuine result and it VINDICATES the
#     composite: different measurement instruments, same underlying microbial
#     signature. Report it. Combined stands.
#
#   PANELS ARE DISJOINT
#     Morphology and biomarkers carry different microbial signals. Combined is
#     averaging two targets and diluting both. The honest output is two
#     indices, reported side by side.
#
#   ONE TARGET YIELDS A PANEL AND THE OTHERS YIELD NOTHING
#     Only that facet carries microbiome signal in this corpus. Report the
#     null for the others - a negative result across 44 studies is a finding,
#     not a failure.
#
# Jaccard overlap is reported alongside raw counts because two panels of very
# different size can share genera without meaning the same thing.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
SD <- .locate()
source(file.path(SD, "00_r_config.R"))

TARGETS <- c("GHI_Combined", "GHI_Morphology", "GHI_Biomarker")

# Where the per-target results are parked so the three runs do not overwrite
# each other. The normal single-target outputs in gdi/ are left untouched.
CMP_DIR <- file.path(R_OUT_DIR, "gdi", "target_comparison")
dir.create(CMP_DIR, recursive = TRUE, showWarnings = FALSE)

GDI_DIR <- file.path(R_OUT_DIR, "gdi")
RSCRIPT <- file.path(R.home("bin"), "Rscript")

# -----------------------------------------------------------------------------
# The target is hard-coded in 05_maaslin2.R (GHI_TARGET) and MaAsLin2 filters
# genera per target, so the comparison must start at 05, not at 06. Rather than
# edit the production scripts, each is copied to a temporary file with the
# target line rewritten. The originals are never modified.
# -----------------------------------------------------------------------------
src05 <- file.path(SD, "05_maaslin2.R")
src06 <- file.path(SD, "06_lasso_selection.R")
for (f in c(src05, src06)) if (!file.exists(f)) stop("Not found: ", f)

# ---- preflight ---------------------------------------------------------------
# Confirm the paths this script's own config resolved to actually exist before
# launching three runs. A path error caught here costs a second; caught after
# the fact it costs the whole job and leaves three identical .err files.
log_msg("[preflight] PDI_ROOT   = %s", PDI_ROOT)
log_msg("[preflight] R_OUT_DIR  = %s", R_OUT_DIR)
if (!dir.exists(MERGED_DIR))
  stop("MERGED_DIR does not exist: ", MERGED_DIR,
       "\n  PDI_ROOT resolved to ", PDI_ROOT,
       "\n  Run this script from the real R_scripts directory.")
ghi_check <- file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")
if (!file.exists(ghi_check))
  stop("GHI scores not found: ", ghi_check, "\n  Run 04b_ghi_pca.R first.")
gh <- fread(ghi_check, nrows = 1)
missing_tg <- setdiff(TARGETS, names(gh))
if (length(missing_tg))
  stop("Target column(s) absent from ", basename(ghi_check), ": ",
       paste(missing_tg, collapse = ", "))
log_msg("[preflight] all three target columns present - OK")

# -----------------------------------------------------------------------------
# 00_r_config.R derives PDI_ROOT from the calling script's own location
# (SCRIPT_DIR/../..). The temp copies must therefore live in the REAL
# R_scripts directory - putting them anywhere else makes PDI_ROOT resolve to
# the wrong place and every data path underneath it breaks. They are written
# with a .tmp_ prefix and removed on exit.
# -----------------------------------------------------------------------------
tmpdir <- SD

TMP_FILES <- character(0)
on.exit({
  for (f in TMP_FILES) if (file.exists(f)) file.remove(f)
}, add = TRUE)

# Running 05 and 06 under each target overwrites the primary gdi_model.rds and
# gdi_stability.tsv in gdi/. They are copied aside here and restored at the end
# so this script is non-destructive wherever it runs in the chain.
.BAK_PRIMARY <- character(0)
for (.f in c("gdi_input.rds", "gdi_model.rds", "gdi_stability.tsv")) {
  .pp <- file.path(GDI_DIR, .f)
  if (file.exists(.pp)) {
    .bb <- paste0(.pp, ".06b_bak")
    file.copy(.pp, .bb, overwrite = TRUE)
    .BAK_PRIMARY <- c(.BAK_PRIMARY, .bb)
  }
}

results <- list()

for (tg in TARGETS) {

  log_msg("")
  log_msg("========================================================")
  log_msg("TARGET: %s", tg)
  log_msg("========================================================")

  # ---- rewrite the target line in a copy of 05 ------------------------------
  l05 <- readLines(src05, warn = FALSE)
  hit <- grep('^GHI_TARGET\\s*<-', l05)
  if (length(hit) != 1)
    stop("Expected exactly one GHI_TARGET assignment in 05_maaslin2.R, found ", length(hit))
  l05[hit] <- sprintf('GHI_TARGET <- "%s"', tg)
  t05 <- file.path(tmpdir, ".tmp_05_target.R")
  writeLines(l05, t05)

  t06 <- file.path(tmpdir, ".tmp_06_target.R")
  file.copy(src06, t06, overwrite = TRUE)
  TMP_FILES <- union(TMP_FILES, c(t05, t06))

  # ---- clear the intermediate so 05 cannot silently reuse the last target ---
  # The already_done pattern in this project has propagated stale inputs
  # before; here it would silently compare a target against itself.
  inp <- file.path(GDI_DIR, "gdi_input.rds")
  if (file.exists(inp)) file.remove(inp)

  ok <- TRUE
  for (scr in c(t05, t06)) {
    log_msg("  running %s", basename(scr))
    rc <- system2(RSCRIPT, shQuote(scr),
                  stdout = file.path(CMP_DIR, sprintf("%s_%s.log", tg, basename(scr))),
                  stderr = file.path(CMP_DIR, sprintf("%s_%s.err", tg, basename(scr))))
    if (rc != 0) {
      log_msg("  *** %s exited %d for target %s - see the .err file", basename(scr), rc, tg)
      ok <- FALSE
      break
    }
  }
  if (!ok) { results[[tg]] <- NULL; next }

  # ---- harvest and park this target's outputs -------------------------------
  mdl <- file.path(GDI_DIR, "gdi_model.rds")
  stb <- file.path(GDI_DIR, "gdi_stability.tsv")
  if (!file.exists(mdl) || !file.exists(stb)) {
    log_msg("  *** expected outputs missing for %s", tg)
    results[[tg]] <- NULL
    next
  }
  file.copy(mdl, file.path(CMP_DIR, sprintf("%s_model.rds", tg)), overwrite = TRUE)
  file.copy(stb, file.path(CMP_DIR, sprintf("%s_stability.tsv", tg)), overwrite = TRUE)

  m <- readRDS(mdl)
  st <- fread(stb)
  results[[tg]] <- list(
    panel   = m$panel,
    freq    = m$freq,
    n_arms  = nrow(m$X),
    n_study = length(unique(m$study)),
    stab    = st,
    max_stb = if (nrow(st)) max(st$stability) else NA_real_,
    n_70    = if (nrow(st)) sum(st$stability >= 0.70) else 0L
  )
  log_msg("  %s: %d arms / %d studies, %d genera at >=70%%, max stability %.0f%%",
          tg, results[[tg]]$n_arms, results[[tg]]$n_study,
          results[[tg]]$n_70, 100 * results[[tg]]$max_stb)
}

got <- names(results)[!vapply(results, is.null, logical(1))]
if (length(got) < 2)
  stop("Fewer than two targets completed; nothing to compare. Check the .err files in ",
       CMP_DIR)

# =============================================================================
# COMPARISON
# =============================================================================
log_msg("")
log_msg("=== COMPARISON ===")

log_msg("")
log_msg("Per target:")
log_msg("  %-16s %6s %7s %8s %9s", "target", "arms", "studies", ">=70%", "max stab")
for (tg in got)
  log_msg("  %-16s %6d %7d %8d %8.0f%%", tg, results[[tg]]$n_arms,
          results[[tg]]$n_study, results[[tg]]$n_70, 100 * results[[tg]]$max_stb)

# ---- pairwise panel overlap --------------------------------------------------
# Jaccard = shared / union. Reported with the raw counts because a large
# overlap between a 2-genus and a 20-genus panel is not the same claim as a
# large overlap between two 10-genus panels.
log_msg("")
log_msg("Panel overlap (>=70%% panels):")
pairs <- combn(got, 2, simplify = FALSE)
ov <- list()
for (p in pairs) {
  a <- results[[p[1]]]$panel; b <- results[[p[2]]]$panel
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  sh <- intersect(a, b); un <- union(a, b)
  j <- if (length(un) == 0) NA_real_ else length(sh) / length(un)
  ov[[paste(p, collapse = " vs ")]] <- list(shared = sh, jac = j,
                                            na = length(a), nb = length(b))
  log_msg("  %-34s %2d shared of %2d/%2d  Jaccard %.2f",
          paste(p, collapse = " vs "), length(sh), length(a), length(b),
          ifelse(is.na(j), 0, j))
  if (length(sh)) log_msg("      %s", paste(sh, collapse = ", "))
}

# ---- rank agreement on the full stability vector -----------------------------
# The >=70% panels may all be empty, in which case overlap above is
# uninformative. Correlating the full stability vectors asks the softer
# question: do the targets AGREE ABOUT WHICH GENERA MATTER, even if none
# clear the threshold? Spearman, because only the ordering is meaningful.
log_msg("")
log_msg("Rank agreement across all genera (Spearman of stability):")
for (p in pairs) {
  A <- results[[p[1]]]$stab; B <- results[[p[2]]]$stab
  m <- merge(A[, .(genus, sa = stability)], B[, .(genus, sb = stability)], by = "genus")
  if (nrow(m) >= 5 && sd(m$sa) > 0 && sd(m$sb) > 0) {
    rho <- cor(m$sa, m$sb, method = "spearman")
    log_msg("  %-34s rho = %+.3f  (%d shared genera)",
            paste(p, collapse = " vs "), rho, nrow(m))
  } else {
    log_msg("  %-34s not computable", paste(p, collapse = " vs "))
  }
}

# ---- reading ------------------------------------------------------------------
log_msg("")
log_msg("Reading:")
jacs <- vapply(ov, function(x) if (is.na(x$jac)) 0 else x$jac, numeric(1))
any_panel <- any(vapply(got, function(tg) length(results[[tg]]$panel) > 0, logical(1)))
n70 <- vapply(got, function(tg) results[[tg]]$n_70, integer(1))

verdict <- if (all(n70 == 0)) {
  paste("No target yields a panel at the registered 70% threshold. The corpus",
        "does not support a stable genus panel against ANY GHI version. Report",
        "the strongest genera with their stability values and treat this as a",
        "negative result, not a tuning problem.")
} else if (max(jacs) >= 0.5) {
  paste("Panels overlap substantially. The facets converge at the microbiome",
        "level despite not converging at the biomarker level. This supports",
        "GHI_Combined as the target and is worth reporting in its own right.")
} else if (max(jacs) <= 0.2) {
  paste("Panels are largely disjoint. The GHI versions carry different",
        "microbial signals. Averaging them into Combined dilutes both;",
        "consider reporting separate indices. This is a protocol decision -",
        "take it to the PI and statistician with this output.")
} else {
  paste("Partial overlap. Read the shared genera and the rank agreement",
        "together; neither convergence nor divergence is clean here.")
}
log_msg("  %s", verdict)

# =============================================================================
# WRITE
# =============================================================================
lines <- c(
  "GDI TARGET COMPARISON - which GHI version does the panel depend on?",
  "All three targets reported. GHI_Combined remains the pre-specified primary;",
  "this is a sensitivity analysis, not a selection procedure.",
  strrep("=", 74), "",
  sprintf("%-16s %6s %7s %8s %9s", "target", "arms", "studies", ">=70%", "max stab"))
for (tg in got)
  lines <- c(lines, sprintf("%-16s %6d %7d %8d %8.0f%%", tg, results[[tg]]$n_arms,
                            results[[tg]]$n_study, results[[tg]]$n_70,
                            100 * results[[tg]]$max_stb))
lines <- c(lines, "", "Panel overlap:")
for (nm in names(ov))
  lines <- c(lines,
    sprintf("  %-34s %2d shared of %2d/%2d  Jaccard %.2f", nm,
            length(ov[[nm]]$shared), ov[[nm]]$na, ov[[nm]]$nb,
            ifelse(is.na(ov[[nm]]$jac), 0, ov[[nm]]$jac)),
    if (length(ov[[nm]]$shared)) sprintf("      %s", paste(ov[[nm]]$shared, collapse = ", ")) else "      (none)")
lines <- c(lines, "", paste("VERDICT:", verdict), "",
  "A formative index does not require its components to correlate. This",
  "compares whether they select the same genera, which is a different and",
  "answerable question. Panel membership at a 70% threshold is a coarse",
  "readout - read the rank agreement alongside it.")

outs <- file.path(CMP_DIR, "target_comparison_summary.txt")
writeLines(lines, outs)
log_msg("")
log_msg("summary  -> %s", outs)
log_msg("per-target outputs -> %s", CMP_DIR)
for (.bb in .BAK_PRIMARY) {
  .pp <- sub("\\.06b_bak$", "", .bb)
  file.copy(.bb, .pp, overwrite = TRUE)
  file.remove(.bb)
}
for (.f in TMP_FILES) if (file.exists(.f)) file.remove(.f)
unlink(list.files(CMP_DIR, pattern = "\\.tmp_.*\\.(log|err)$", full.names = TRUE))
log_msg("primary gdi_model.rds and gdi_stability.tsv restored to the %s run", TARGETS[1])
log_msg("=== DONE ===")
