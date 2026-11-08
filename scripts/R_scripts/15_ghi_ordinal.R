#!/usr/bin/env Rscript
# =============================================================================
# 15_ghi_ordinal.R
#
# Does the GHI move monotonically across challenge, control and additive?
#
# RUN:  Rscript 15_ghi_ordinal.R     (needs 04b to have run)
# TIME: under a minute.
#
# -----------------------------------------------------------------------------
# WHY THIS IS A DIFFERENT TEST FROM 04e
# -----------------------------------------------------------------------------
# 04e_mixed_model_validation.R runs two separate pairwise contrasts: additive
# against control, and control against challenge. Each is its own
# meta-analysis with its own k and its own p-value. That answers two questions
# - is additive better than control, is control better than challenge - but
# not the single question of whether the index moves in order across all three
# groups at once.
#
# Here every arm is coded on one ordinal scale:
#
#     challenge = -1        control = 0        additive = +1
#
# and the GHI is regressed on that single scale within study. One slope, one
# p-value, testing whether the index rises monotonically from challenge
# through control to additive - the stronger, unified claim implied by scoring
# three groups on a single health axis.
#
# Arms outside these three categories (Challenged_And_Treated,
# Phenotype_Comparison, or missing a label) are excluded, since they do not
# have a place on this scale by construction.
#
# -----------------------------------------------------------------------------
# HOW THE TEST IS RUN
# -----------------------------------------------------------------------------
# Within each study, both the ordinal code and the GHI are centred on that
# study's own mean before fitting, which is the same within-study logic used
# throughout this project: it removes the study's baseline so that only the
# within-study ordering is tested, and it makes the slope comparable across
# studies measured on different original scales.
#
# The model is a simple within-study linear fit of GHI on the ordinal code,
# with permutation inference: the ordinal code is shuffled within each study
# 5,000 times, preserving the number of arms and the clustering of arms within
# studies, to build a null distribution for the slope.
#
# A second version uses ONLY the ordinal ranks (-1, 0, +1) with ties broken by
# the mean, and reports Spearman rho as a distribution-free alternative to the
# slope.
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
N_PERM <- 5000
set.seed(42)

log_msg("[1/4] Loading")
GHI_TSV <- file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")
if (!file.exists(GHI_TSV)) stop("Not found: ", GHI_TSV, "\n  Run 04b_ghi_pca.R first.")
gh <- fread(GHI_TSV)

cls <- trimws(as.character(gh$Challenge_Type))
code <- rep(NA_real_, length(cls))
code[cls %in% CHALLENGE] <- -1
code[cls %in% HEALTHY]   <-  0
code[cls %in% ADDITIVE]  <- +1

log_msg("  %d rows in the GHI table", nrow(gh))
log_msg("  on the ordinal scale: %d challenge (-1), %d control (0), %d additive (+1)",
        sum(code == -1, na.rm = TRUE), sum(code == 0, na.rm = TRUE),
        sum(code == 1, na.rm = TRUE))
log_msg("  excluded (Challenged_And_Treated, Phenotype_Comparison, or unlabelled): %d",
        sum(is.na(code)))

study <- trimws(as.character(gh$Study_ID))

run_one <- function(y, label) {
  ok <- !is.na(code) & !is.na(y)
  x <- code[ok]; yy <- y[ok]; g <- study[ok]
  tb <- table(g); ok2 <- g %in% names(tb)[tb >= 2]
  x <- x[ok2]; yy <- yy[ok2]; g <- g[ok2]

  # within-study centring
  for (s in unique(g)) {
    i <- g == s
    x[i]  <- x[i]  - mean(x[i])
    yy[i] <- yy[i] - mean(yy[i])
  }
  # a study contributing only one distinct ordinal value has nothing left
  # after centring on that value; drop it
  keep_s <- vapply(unique(g), function(s) length(unique(round(x[g == s], 8))) > 1, logical(1))
  keep_s <- names(keep_s)[keep_s]
  m <- g %in% keep_s
  x <- x[m]; yy <- yy[m]; g <- g[m]

  n_arms <- length(x); n_st <- length(unique(g))
  if (n_arms < 15 || n_st < 5)
    return(list(label = label, n = n_arms, k = n_st, slope = NA_real_,
               p = NA_real_, rho = NA_real_, rho_p = NA_real_))

  fit <- stats::lm(yy ~ x - 1)          # through the origin; both centred
  slope <- unname(coef(fit)[1])

  idx <- split(seq_along(x), g)
  obs_slope <- slope
  cnt <- 0L
  for (i in seq_len(N_PERM)) {
    xp <- x
    for (j in idx) if (length(j) > 1) xp[j] <- sample(x[j])
    s2 <- unname(coef(stats::lm(yy ~ xp - 1))[1])
    if (abs(s2) >= abs(obs_slope)) cnt <- cnt + 1L
  }
  p_perm <- (cnt + 1) / (N_PERM + 1)

  rho <- suppressWarnings(stats::cor(x, yy, method = "spearman"))
  rho_ct <- suppressWarnings(stats::cor.test(x, yy, method = "spearman"))

  list(label = label, n = n_arms, k = n_st, slope = slope, p = p_perm,
       rho = rho, rho_p = rho_ct$p.value)
}

log_msg("[2/4] Fitting the ordinal model, three GHI versions")
targets <- list(Morphology = gh$GHI_Morphology_raw,
                Biomarker  = gh$GHI_Biomarker_raw,
                Combined   = gh$GHI_Combined_raw)
res <- lapply(names(targets), function(nm) run_one(targets[[nm]], nm))
names(res) <- names(targets)

log_msg("")
log_msg("  %-11s %6s %7s %10s %10s %10s %10s", "index", "n", "studies",
        "slope", "perm p", "rho", "rho p")
for (nm in names(res)) {
  z <- res[[nm]]
  log_msg("  %-11s %6d %7d %10s %10s %10s %10s", nm, z$n, z$k,
          ifelse(is.na(z$slope), "n/a", sprintf("%+.3f", z$slope)),
          ifelse(is.na(z$p), "n/a", sprintf("%.4f", z$p)),
          ifelse(is.na(z$rho), "n/a", sprintf("%+.3f", z$rho)),
          ifelse(is.na(z$rho_p), "n/a", sprintf("%.4f", z$rho_p)))
}
log_msg("")
log_msg("  Positive slope/rho is the working direction: challenge (-1) lowest,")
log_msg("  additive (+1) highest, control (0) between. Both x and y are")
log_msg("  centred within study before fitting.")

log_msg("[3/4] Reading")
c1 <- res[["Combined"]]
verdict <- if (is.na(c1$p)) {
  "Not evaluable for Combined."
} else if (c1$slope > 0 && c1$p < 0.05) {
  sprintf(paste("The index moves monotonically across the three groups within",
                "study (slope = %+.3f, permutation p = %.4f). This is a single,",
                "unified ordering test rather than two separate pairwise",
                "contrasts, and it is significant."), c1$slope, c1$p)
} else if (c1$slope > 0) {
  sprintf(paste("The index moves in the expected direction (slope = %+.3f) but",
                "the ordinal test does not reach significance (p = %.4f). Two",
                "pairwise contrasts (04e) can each be significant while the",
                "combined ordering is not, if the additive and challenge effects",
                "are not proportional to the +1/-1/0 coding assumed here."),
          c1$slope, c1$p)
} else {
  sprintf("The index does not move in the expected direction (slope = %+.3f).", c1$slope)
}
log_msg("  %s", verdict)

log_msg("[4/4] Writing")
out <- c("GHI AGAINST A SINGLE CHALLENGE-CONTROL-ADDITIVE SCALE", strrep("=", 74), "",
  "Every arm coded challenge = -1, control = 0, additive = +1. The GHI and the",
  "code are both centred within study, then regressed (GHI ~ code) through the",
  "origin. Significance by permutation, shuffling the code within study 5,000",
  "times. Spearman rho on the same centred quantities is reported alongside as",
  "a distribution-free check.", "",
  sprintf("  %-11s %6s %7s %10s %10s %10s %10s", "index", "n", "studies",
          "slope", "perm p", "rho", "rho p"))
for (nm in names(res)) {
  z <- res[[nm]]
  out <- c(out, sprintf("  %-11s %6d %7d %10s %10s %10s %10s", nm, z$n, z$k,
                        ifelse(is.na(z$slope), "n/a", sprintf("%+.3f", z$slope)),
                        ifelse(is.na(z$p), "n/a", sprintf("%.4f", z$p)),
                        ifelse(is.na(z$rho), "n/a", sprintf("%+.3f", z$rho)),
                        ifelse(is.na(z$rho_p), "n/a", sprintf("%.4f", z$rho_p))))
}
out <- c(out, "", verdict, "",
  "This test asks whether the index orders challenge < control < additive on",
  "one scale. 04e_mixed_model_validation.R asks the related but distinct",
  "question of whether each pairwise contrast individually differs from zero.",
  "The two can disagree if the additive and challenge effects are unequal in",
  "size, since the ordinal test assumes they are symmetric around control.")
writeLines(out, file.path(R_OUT_DIR, "ghi", "ghi_ordinal.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "ghi", "ghi_ordinal.txt"))
log_msg("=== DONE ===")
