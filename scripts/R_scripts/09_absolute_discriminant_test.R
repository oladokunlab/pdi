#!/usr/bin/env Rscript
# =============================================================================
# 09_absolute_discriminant_test.R
#
# Can a canine-DI-style ABSOLUTE formula be built from this corpus?
#
# RUN:  Rscript 09_absolute_discriminant_test.R   (needs 05, 06 first)
# TIME: a few minutes. Reads only; writes a summary.
#
# -----------------------------------------------------------------------------
# THE QUESTION
# -----------------------------------------------------------------------------
# The GDI as built needs within-study centring: a sample is scored relative to
# the mean of its own study. That makes it a within-trial ranking tool. LOSO
# gave rho = -0.309 centred but only -0.093 absolute, so the current
# construction does not support scoring one isolated sample.
#
# The canine and feline Dysbiosis Indices do exactly that - one qPCR panel, a
# fixed formula, a cutoff of 0 - and they manage it for two reasons this
# script tests separately:
#
#   1. They are fitted against a CLASS LABEL (healthy vs chronic enteropathy),
#      not against a continuous index. A label is study-independent in a way
#      that a biomarker composite is not.
#   2. They use a LINEAR DISCRIMINANT on absolute abundances, not a distance
#      to a centroid.
#
# So the failure of the absolute GDI might be caused by the target (a
# study-dependent GHI), by the construction (centroid distance), or by the
# data genuinely lacking a study-independent signal. Those imply very
# different next steps, and guessing between them is not necessary.
#
# -----------------------------------------------------------------------------
# WHAT IS TESTED, AND WHY EACH ONE
# -----------------------------------------------------------------------------
#   A. CENTROID + GHI, absolute          the current construction, as baseline
#   B. LOGISTIC + LABEL, absolute        the canine DI design: fit control vs
#                                        challenged directly on absolute CLR
#   C. LOGISTIC + LABEL, centred         same target, centred inputs - isolates
#                                        whether centring or the target matters
#
# All three are evaluated by leave-one-study-out, because in-sample AUC on 135
# arms will look good regardless and would be meaningless here.
#
# -----------------------------------------------------------------------------
# HOW TO READ THE RESULT
# -----------------------------------------------------------------------------
#   B clearly beats A   -> the TARGET was the problem. A label-trained absolute
#                          index is worth building properly; that is a genuine
#                          route to a canine-DI-style assay.
#   B is no better       -> the corpus has no study-independent signal at this
#                          taxonomic resolution. Report the within-trial index
#                          and state the absolute limitation. No amount of
#                          reformulation fixes it.
#   C good, B poor       -> the signal is real but only within study. Same
#                          conclusion as above for deployment purposes, but it
#                          confirms the biology is there.
#
# A negative result here is informative and worth reporting: it defines what a
# 16S-derived poultry index can and cannot do, which nobody has established.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

set.seed(42)
UNASSIGNABLE <- c("g__uncultured", "g__unclassified", "g__uncultured_bacterium",
                  "g__metagenome", "g__Incertae_Sedis", "g__unidentified")
HEALTHY   <- c("Control")
CHALLENGE <- c("Pathogen_Challenge", "Antibiotic_Challenge", "Dietary_Feed_Challenge")

log_msg("[1/5] Loading")
mf <- file.path(R_OUT_DIR, "gdi", "gdi_model.rds")
if (!file.exists(mf)) stop("Not found: ", mf, "\n  Run 06_lasso_selection.R first.")
m <- readRDS(mf)
X <- m$X; y <- m$y; study <- m$study
arm <- rownames(X)
log_msg("  %d arms, %d genera, %d studies", nrow(X), ncol(X), length(unique(study)))

# ---- class label from the curated column ------------------------------------
cls <- rep(NA_character_, length(arm))
gt <- file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")
if (file.exists(gt)) {
  gh <- fread(gt)
  if ("Challenge_Type" %in% names(gh)) {
    k <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
    lut <- tapply(trimws(as.character(gh$Challenge_Type)), k, function(v) {
      v <- v[!is.na(v) & nzchar(v)]
      if (!length(v)) NA_character_ else names(sort(table(v), decreasing = TRUE))[1] })
    cls <- unname(lut[arm])
  }
}
lab <- ifelse(cls %in% HEALTHY, 0L, ifelse(cls %in% CHALLENGE, 1L, NA_integer_))
log_msg("  labelled arms: %d control, %d challenged", sum(lab %in% 0L), sum(lab %in% 1L))
if (sum(lab %in% 0L) < 15 || sum(lab %in% 1L) < 15)
  stop("Too few labelled arms for a label-trained model.")

# Studies that carry BOTH classes are the only ones informative for a
# within-study version, and are also the only ones where a fold can be scored
# on both sides. Reported so the reduction is visible.
lab_ok <- !is.na(lab)
both <- names(which(tapply(lab[lab_ok], study[lab_ok],
                           function(v) length(unique(v)) == 2)))
log_msg("  studies with both classes: %d of %d",
        length(both), length(unique(study[lab_ok])))

center_within <- function(M, g) {
  out <- M
  for (s in unique(g)) { i <- g == s
    out[i, ] <- sweep(M[i, , drop = FALSE], 2, colMeans(M[i, , drop = FALSE]), "-") }
  out
}
auc_of <- function(sc, y1) {
  s1 <- sc[y1 == 1]; s0 <- sc[y1 == 0]
  if (length(s1) < 3 || length(s0) < 3) return(NA_real_)
  r <- rank(c(s1, s0))
  (sum(r[seq_along(s1)]) - length(s1)*(length(s1)+1)/2) / (length(s1)*length(s0))
}
dist_to <- function(r, c0) sqrt(sum((r - c0)^2))

# =============================================================================
# A - CURRENT CONSTRUCTION: centroid on GHI, scored absolute
# =============================================================================
log_msg("[2/5] A: centroid trained on GHI, scored on ABSOLUTE values (baseline)")
panel_pub <- setdiff(m$panel, UNASSIGNABLE)
scoreA <- rep(NA_real_, length(arm))
for (s in unique(study)) {
  tr <- study != s; te <- !tr
  Xtr <- X[tr, panel_pub, drop = FALSE]; ytr <- y[tr]
  hh <- ytr > stats::median(ytr)
  if (sum(hh) < 3 || sum(!hh) < 3) next
  ch <- colMeans(Xtr[hh, , drop = FALSE]); cd <- colMeans(Xtr[!hh, , drop = FALSE])
  scoreA[te] <- apply(X[te, panel_pub, drop = FALSE], 1,
                      function(r) dist_to(r, ch) - dist_to(r, cd))
}
kA <- !is.na(scoreA) & !is.na(lab)
aucA <- auc_of(scoreA[kA], lab[kA])
log_msg("  LOSO AUC (challenged vs control) = %s",
        if (is.na(aucA)) "n/a" else sprintf("%.3f", aucA))

# =============================================================================
# B - CANINE-DI DESIGN: penalised logistic on the LABEL, absolute inputs
#
# L1-penalised so the panel stays small enough to be a qPCR assay. Fitted only
# on labelled arms, and the held-out study is never in the training set.
# =============================================================================
log_msg("[3/5] B: logistic on the LABEL, ABSOLUTE inputs (the canine-DI design)")
fit_logit_loso <- function(Xin) {
  sc <- rep(NA_real_, length(arm))
  sel_count <- setNames(numeric(ncol(Xin)), colnames(Xin)); nfold <- 0
  for (s in unique(study[!is.na(lab)])) {
    tr <- !is.na(lab) & study != s
    te <- !is.na(lab) & study == s
    if (sum(te) < 2) next
    if (length(unique(lab[tr])) < 2) next
    if (sum(lab[tr] == 1) < 8 || sum(lab[tr] == 0) < 8) next
    f <- tryCatch(cv.glmnet(Xin[tr, , drop = FALSE], lab[tr], family = "binomial",
                            alpha = 1, nfolds = 5, standardize = TRUE),
                  error = function(e) NULL)
    if (is.null(f)) next
    sc[te] <- as.numeric(predict(f, newx = Xin[te, , drop = FALSE], s = f$lambda.1se))
    co <- as.matrix(coef(f, s = f$lambda.1se))
    nz <- rownames(co)[co[, 1] != 0]; nz <- setdiff(nz, "(Intercept)")
    sel_count[nz] <- sel_count[nz] + 1; nfold <- nfold + 1
  }
  list(score = sc, sel = sel_count, nfold = nfold)
}
Xa <- X[, setdiff(colnames(X), UNASSIGNABLE), drop = FALSE]
B <- fit_logit_loso(Xa)
kB <- !is.na(B$score) & !is.na(lab)
aucB <- auc_of(B$score[kB], lab[kB])
log_msg("  %d folds fitted; LOSO AUC = %s", B$nfold,
        if (is.na(aucB)) "n/a" else sprintf("%.3f", aucB))
if (B$nfold > 0) {
  top <- sort(B$sel[B$sel > 0], decreasing = TRUE)
  log_msg("  genera selected in >=50%% of folds:")
  shown <- FALSE
  for (g in names(top)) if (top[[g]] / B$nfold >= 0.5) {
    log_msg("    %-40s %d/%d", g, top[[g]], B$nfold); shown <- TRUE }
  if (!shown) log_msg("    none")
}

# =============================================================================
# C - SAME TARGET, CENTRED INPUTS
# Isolates whether any shortfall in B comes from the target or from centring.
# =============================================================================
log_msg("[4/5] C: logistic on the LABEL, CENTRED inputs")
Xc <- center_within(Xa, study)
Cc <- fit_logit_loso(Xc)
kC <- !is.na(Cc$score) & !is.na(lab)
aucC <- auc_of(Cc$score[kC], lab[kC])
log_msg("  %d folds fitted; LOSO AUC = %s", Cc$nfold,
        if (is.na(aucC)) "n/a" else sprintf("%.3f", aucC))

# =============================================================================
# VERDICT
# =============================================================================
log_msg("[5/5] Reading")
log_msg("  A  centroid / GHI  / absolute : %s", if (is.na(aucA)) "n/a" else sprintf("%.3f", aucA))
log_msg("  B  logistic / label / absolute : %s", if (is.na(aucB)) "n/a" else sprintf("%.3f", aucB))
log_msg("  C  logistic / label / centred  : %s", if (is.na(aucC)) "n/a" else sprintf("%.3f", aucC))

verdict <- if (is.na(aucB)) {
  "Model B was not evaluable; too few labelled arms per fold."
} else if (aucB >= 0.70) {
  paste("Model B reaches conventional discrimination. A label-trained model on",
        "absolute abundances is achievable in this corpus, indicating that the",
        "weaker absolute performance of the centroid index follows from",
        "regression against a study-dependent continuous outcome rather than",
        "from the abundance data themselves.")
} else if (!is.na(aucA) && aucB - aucA >= 0.05) {
  paste("Model B improves on model A but remains below 0.70. The choice of",
        "outcome accounts for part of the shortfall in absolute performance,",
        "though discrimination is insufficient for single-sample application.")
} else {
  paste("Model B does not improve on model A. No study-independent",
        "discrimination at genus level was recovered by either an alternative",
        "outcome or an alternative model form. The limitation therefore",
        "reflects the resolution of the abundance data rather than the",
        "construction of the index.")
}
log_msg("  %s", verdict)
if (!is.na(aucC) && !is.na(aucB) && aucC - aucB >= 0.10)
  log_msg("  C beats B by %.2f: the signal is real but lives within study.", aucC - aucB)

out <- c("ABSOLUTE DISCRIMINANT FEASIBILITY",
  strrep("=", 74), "",
  "OBJECTIVE", strrep("-", 74),
  "Assessment of whether a fixed-coefficient index applicable to a single",
  "sample, comparable in form to the canine and feline Dysbiosis Indices, can",
  "be derived from this corpus. Three constructions were compared, differing in",
  "outcome and in model form.", "",
  sprintf("Labelled arms                   : %d control, %d challenged",
          sum(lab %in% 0L), sum(lab %in% 1L)),
  sprintf("Studies containing both classes : %d", length(both)), "",
  "RESULTS", strrep("-", 74),
  "Area under the ROC curve, challenged against control, leave-one-study-out.", "",
  sprintf("  A  centroid, continuous outcome, absolute abundances : %s",
          ifelse(is.na(aucA), "n/a", sprintf("%.3f", aucA))),
  sprintf("  B  penalised logistic, class label, absolute         : %s",
          ifelse(is.na(aucB), "n/a", sprintf("%.3f", aucB))),
  sprintf("  C  penalised logistic, class label, within-study     : %s",
          ifelse(is.na(aucC), "n/a", sprintf("%.3f", aucC))),
  "", strwrap(verdict, width = 74), "",
  "NOTES", strrep("-", 74),
  "1. All estimates are leave-one-study-out. Apparent performance on the full",
  sprintf("   %d arms is not reported, as it does not bear on transfer to an", nrow(X)),
  "   independent study.",
  "",
  "2. Models B and C were fitted only on arms carrying a control or challenge",
  "   label, and only in folds where the training set retained both classes,",
  "   which reduces the number of evaluable folds relative to model A.",
  "",
  "3. Placeholder taxonomic labels were excluded from all three models.")
writeLines(out, file.path(R_OUT_DIR, "gdi", "absolute_discriminant_test.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "absolute_discriminant_test.txt"))
log_msg("=== DONE ===")
