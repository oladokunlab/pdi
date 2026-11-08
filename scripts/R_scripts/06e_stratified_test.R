#!/usr/bin/env Rscript
# =============================================================================
# 06e_stratified_test.R
#
# Does the genus signal appear when the contrast is restricted, or when age is
# adjusted for?
#
# RUN:  Rscript 06e_stratified_test.R    (needs 05 to have run)
# TIME: 5-15 minutes.
#
# -----------------------------------------------------------------------------
# WHY THIS IS NOT THRESHOLD TUNING
# -----------------------------------------------------------------------------
# After the morphology correction no genus exceeded 15% selection stability
# against any GHI version. The response to that must not be to lower the
# threshold until something passes. It can legitimately be to ask whether the
# model was being fitted where the outcome is meaningful.
#
# TWO REASONS TO THINK IT WAS NOT.
#
# The GHI detects the additive contrast (+0.766, p < 0.0001) and does not
# detect the challenge contrast (+0.157, p = 0.51); the difference between them
# is significant at p = 0.009. For challenged arms the outcome therefore
# carries little information about gut health. Roughly a third of the arms
# entering the LASSO are of that kind. Fitting a predictor against a target
# that is uninformative on part of the data adds noise to selection, and
# restricting to the contrast where the target demonstrably works is a
# principled restriction rather than a search for a favourable subset.
#
# Second, 18 of 48 studies sampled at more than one age, and caecal community
# composition changes substantially between day 7 and day 42. Arms average over
# ages, which absorbs much of this, but where age varies within a study the
# residual within-study contrast carries developmental change alongside
# treatment effect.
#
# -----------------------------------------------------------------------------
# WHAT IS COMPARED
# -----------------------------------------------------------------------------
#   A  all arms                      the current analysis, as the baseline
#   B  additive and control only     the contrast the GHI is validated for
#   C  challenge and control only    the contrast it is not - reported so the
#                                    comparison is symmetric rather than
#                                    selective
#   D  all arms, age-adjusted        age regressed out of the outcome within
#                                    study before fitting
#
# Every configuration uses the same selection rule and the same registered
# threshold. Only the data entering the fit changes, and all four are reported
# whatever they show.
#
# -----------------------------------------------------------------------------
# WHAT WOULD COUNT AS A RESULT
# -----------------------------------------------------------------------------
# B reaching a stability clearly above A would indicate the signal was being
# diluted by arms whose outcome is uninformative. That is a finding, and it
# would need reporting as a restriction with its own justification and its own
# leave-one-study-out validation.
#
# B resembling A means the dilution explanation is wrong and the negative
# result stands. Expect this: the collapse was from 52% to 15%, which is large,
# and B has roughly 40% fewer arms, so it loses power as it gains precision.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

N_BOOT   <- 1000
SUBFRAC  <- 0.50
FLOOR    <- 0.20
REGISTER <- 0.70
set.seed(42)
UNASSIGNABLE <- c("g__uncultured", "g__unclassified", "g__uncultured_bacterium",
                  "g__metagenome", "g__Incertae_Sedis", "g__unidentified")
HEALTHY   <- "Control"
CHALLENGE <- c("Pathogen_Challenge", "Antibiotic_Challenge", "Dietary_Feed_Challenge")
ADDITIVE  <- "Non_Disruptive_Additive"

log_msg("[1/3] Loading")
mf <- file.path(R_OUT_DIR, "gdi", "gdi_model.rds")
if (!file.exists(mf)) stop("Not found: ", mf, "\n  Run 06_lasso_selection.R first.")
m <- readRDS(mf)
X <- m$X; y <- m$y; study <- m$study; arm <- rownames(X)
log_msg("  %d arms, %d genera, %d studies", nrow(X), ncol(X), length(unique(study)))

# ---- intervention class and age, per arm ------------------------------------
gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
k_gh <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
lut <- function(col) {
  v <- tapply(gh[[col]], k_gh, function(z) {
    z <- z[!is.na(z)]
    if (!length(z)) return(NA)
    if (is.numeric(z)) mean(z) else names(sort(table(z), decreasing = TRUE))[1]
  })
  unname(v[arm])
}
cls <- if ("Challenge_Type" %in% names(gh)) trimws(as.character(lut("Challenge_Type"))) else rep(NA_character_, length(arm))
age <- if ("Age_d" %in% names(gh)) suppressWarnings(as.numeric(lut("Age_d"))) else rep(NA_real_, length(arm))
log_msg("  class: %d control, %d additive, %d challenged, %d other",
        sum(cls %in% HEALTHY), sum(cls %in% ADDITIVE),
        sum(cls %in% CHALLENGE), sum(is.na(cls)))
log_msg("  age available for %d of %d arms", sum(!is.na(age)), length(age))

get_coefs <- function(fit, l) {
  co <- as.matrix(coef(fit, s = l))
  co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  co[co[, 1] != 0, 1]
}

# ---- stability selection on an arbitrary subset ------------------------------
run_sel <- function(keep, yv, label) {
  Xs <- X[keep, , drop = FALSE]; ys <- yv[keep]; ss <- study[keep]
  tb <- table(ss); ok <- ss %in% names(tb)[tb >= 2]
  Xs <- Xs[ok, , drop = FALSE]; ys <- ys[ok]; ss <- ss[ok]
  if (length(unique(ss)) < 6 || nrow(Xs) < 25)
    return(list(label = label, n = nrow(Xs), k = length(unique(ss)),
                max = NA_real_, n70 = 0L, n20 = 0L, top = character(0), r2 = NA_real_))
  Xc <- Xs; yc <- ys
  for (s in unique(ss)) { i <- ss == s
    Xc[i, ] <- sweep(Xs[i, , drop = FALSE], 2, colMeans(Xs[i, , drop = FALSE]), "-")
    yc[i]   <- ys[i] - mean(ys[i]) }
  st <- unique(ss); ntake <- max(3, floor(length(st) * SUBFRAC))
  cnt <- setNames(numeric(ncol(Xc)), colnames(Xc)); nok <- 0
  for (b in seq_len(N_BOOT)) {
    tk <- sample(st, ntake); idx <- which(ss %in% tk)
    if (length(idx) < 10 || length(unique(yc[idx])) < 5) next
    f <- tryCatch(cv.glmnet(Xc[idx, , drop = FALSE], yc[idx], alpha = 1,
                            nfolds = min(5, length(idx)), standardize = TRUE,
                            intercept = FALSE), error = function(e) NULL)
    if (is.null(f)) next
    sel <- get_coefs(f, f$lambda.1se)
    cnt[names(sel)] <- cnt[names(sel)] + 1; nok <- nok + 1
  }
  if (nok < 50) return(list(label = label, n = nrow(Xc), k = length(unique(ss)),
                            max = NA_real_, n70 = 0L, n20 = 0L, top = character(0), r2 = NA_real_))
  fr <- cnt / nok
  fr <- fr[!(names(fr) %in% UNASSIGNABLE)]
  r2 <- tryCatch({
    cvf <- cv.glmnet(Xc, yc, alpha = 1, standardize = TRUE, intercept = FALSE)
    pr <- as.numeric(predict(cvf, newx = Xc, s = cvf$lambda.1se))
    if (stats::sd(pr) == 0) 0 else max(0, 1 - sum((yc - pr)^2) / sum((yc - mean(yc))^2))
  }, error = function(e) NA_real_)
  list(label = label, n = nrow(Xc), k = length(unique(ss)),
       max = max(fr), n70 = sum(fr >= REGISTER), n20 = sum(fr >= FLOOR),
       top = names(sort(fr, decreasing = TRUE))[1:min(6, length(fr))],
       topv = sort(fr, decreasing = TRUE)[1:min(6, length(fr))], r2 = r2)
}

# ---- age-adjusted outcome ----------------------------------------------------
# Age is regressed out of the outcome within each study, so that what remains
# is variation in gut health not attributable to the birds being older. Only
# studies with more than one sampling age are affected; elsewhere age is
# constant and the residual equals the centred outcome.
y_adj <- y
if (sum(!is.na(age)) > 20) {
  for (s in unique(study)) {
    i <- study == s & !is.na(age) & !is.na(y)
    if (sum(i) >= 3 && length(unique(age[i])) > 1) {
      fit <- stats::lm(y[i] ~ age[i])
      y_adj[i] <- stats::residuals(fit) + mean(y[i])
    }
  }
  n_adj <- sum(vapply(unique(study), function(s) {
    i <- study == s & !is.na(age); sum(i) >= 3 && length(unique(age[i])) > 1 }, logical(1)))
  log_msg("  age adjustment applied within %d studies sampling at multiple ages", n_adj)
}

log_msg("[2/3] Running four configurations")
cfg <- list(
  run_sel(rep(TRUE, nrow(X)),                         y,     "A  all arms"),
  run_sel(cls %in% c(HEALTHY, ADDITIVE),              y,     "B  additive + control"),
  run_sel(cls %in% c(HEALTHY, CHALLENGE),             y,     "C  challenge + control"),
  run_sel(rep(TRUE, nrow(X)),                         y_adj, "D  all arms, age-adjusted"))

log_msg("[3/3] Comparison")
log_msg("")
log_msg("  %-26s %5s %6s %9s %7s %7s %7s", "configuration", "arms", "studies",
        "max stab", ">=20%", ">=70%", "R2")
for (z in cfg)
  log_msg("  %-26s %5d %6d %8s %7d %7d %7s", z$label, z$n, z$k,
          if (is.na(z$max)) "n/a" else sprintf("%.0f%%", 100*z$max),
          z$n20, z$n70, if (is.na(z$r2)) "n/a" else sprintf("%.3f", z$r2))
for (z in cfg) if (length(z$top)) {
  log_msg("")
  log_msg("  %s - strongest genera:", z$label)
  for (i in seq_along(z$top))
    log_msg("      %-40s %.0f%%", gsub("^g__", "", z$top[i]), 100*z$topv[i])
}

base <- cfg[[1]]$max; bt <- cfg[[2]]$max
log_msg("")
verdict <- if (is.na(bt) || is.na(base)) {
  "Not comparable - one configuration had too few studies to fit."
} else if (bt >= REGISTER) {
  "B reaches the registered threshold. Restricting to the validated contrast recovers a panel; it must now be validated by leave-one-study-out and reported as a restriction with its own justification."
} else if (bt - base >= 0.15) {
  sprintf(paste("B is markedly stronger than A (%.0f%% against %.0f%%) without reaching the",
                "registered threshold. Consistent with dilution by arms whose outcome is",
                "uninformative, and worth reporting, but not a panel."), 100*bt, 100*base)
} else {
  sprintf(paste("B resembles A (%.0f%% against %.0f%%). Restricting to the validated contrast",
                "does not recover the signal, so dilution by challenge arms does not explain",
                "the negative result. The finding stands."), 100*bt, 100*base)
}
log_msg("  %s", verdict)

out <- c("STRATIFIED AND AGE-ADJUSTED SELECTION", strrep("=", 74), "",
  "Same selection rule and same registered threshold throughout; only the data",
  "entering the fit differs. All four configurations reported.", "",
  sprintf("  %-26s %5s %6s %9s %7s %7s %7s", "configuration", "arms", "studies",
          "max stab", ">=20%", ">=70%", "R2"))
for (z in cfg)
  out <- c(out, sprintf("  %-26s %5d %6d %8s %7d %7d %7s", z$label, z$n, z$k,
           ifelse(is.na(z$max), "n/a", sprintf("%.0f%%", 100*z$max)), z$n20, z$n70,
           ifelse(is.na(z$r2), "n/a", sprintf("%.3f", z$r2))))
out <- c(out, "", verdict, "",
  "Configuration B restricts to the contrast the GHI is validated for; the",
  "index detects the additive contrast and does not detect the challenge one.",
  "Configuration D removes within-study age variation from the outcome, which",
  "affects only the studies sampling at more than one age.")
writeLines(out, file.path(R_OUT_DIR, "gdi", "stratified_test.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "stratified_test.txt"))
log_msg("=== DONE ===")
