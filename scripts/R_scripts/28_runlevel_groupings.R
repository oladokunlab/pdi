#!/usr/bin/env Rscript
# =============================================================================
# 28_runlevel_groupings.R
#
# Run-level GDI under two additional groupings of the treatment label. No GHI
# used anywhere. Companion to 27_runlevel_gdi.R, which tested control vs
# challenge alone.
#
# RUN:  Rscript 28_runlevel_groupings.R
# TIME: 20-40 minutes (two full analyses, 1000 resamples each).
#
# -----------------------------------------------------------------------------
# WHY THE GROUPING IS A REAL QUESTION AND NOT A DETAIL
# -----------------------------------------------------------------------------
# 27_runlevel_gdi.R used control vs challenge only, on the reasoning that a
# supplemented bird is not the same state as an untreated one. That is
# defensible, and it excluded the largest group in the corpus: only 13 studies
# contain both a control and a challenged arm.
#
# There is an equally defensible argument the other way. In the field the
# question is "is this gut healthy", not "was this bird challenged". A control
# bird and a supplemented bird are both NOT dysbiotic. Grouping them together
# is arguably closer to the clinical question, and it lets any study with an
# additive arm and a challenge arm contribute.
#
# Rather than pick one, both are run here so the choice becomes a reported
# result instead of an unexamined assumption.
#
#   B  healthy (control + additive) vs challenged     - binary
#   C  control vs additive vs challenge               - three-class, run as
#      three one-vs-rest classifiers, since a single ordered model would
#      assume the three classes lie on one line, which is exactly what is in
#      question
#
# -----------------------------------------------------------------------------
# WHAT TO EXPECT, STATED IN ADVANCE
# -----------------------------------------------------------------------------
# Grouping B should raise the study count, because the qualifying condition
# becomes easier to meet. It will NOT fix the underlying constraint: the
# corpus has roughly twenty studies with a challenge arm at all, and no
# grouping creates challenge studies that do not exist.
#
# A related caution: 20_control_vs_treated.R already pooled in the opposite
# direction - control vs everything treated - and returned AUC 0.450, below
# chance. Different grouping, but a reason not to expect much.
#
# The number that decides anything is LOSO AUC on held-out studies. In-sample
# AUC is reported alongside it purely so the gap is visible; a large gap means
# the model is fitting study-specific structure rather than biology, which is
# what 27 found (0.775 in-sample against 0.471 held out).
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

N_BOOT <- 1000; SUBFRAC <- 0.50; STABILITY_THR <- 0.70
set.seed(42)
UNASSIGNABLE <- c("g__uncultured","g__unclassified","g__uncultured_bacterium",
                  "g__metagenome","g__Incertae_Sedis","g__unidentified")
HEALTHY   <- "Control"
ADDITIVE  <- "Non_Disruptive_Additive"
CHALLENGE <- c("Pathogen_Challenge","Antibiotic_Challenge","Dietary_Feed_Challenge")

log_msg("[1/3] Loading run-level data and treatment labels")
clr <- fread(file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv"))
M <- as.matrix(clr[, -1, with = FALSE]); rownames(M) <- clr[[names(clr)[1]]]
md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2s <- setNames(trimws(as.character(md$Study_ID)), rid)
r2t <- setNames(trimws(as.character(md$Treatment_Group)), rid)
runs <- intersect(colnames(M), names(r2s))
M <- M[, runs, drop = FALSE]

bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
bk <- paste(trimws(bio$Study_ID), trimws(bio$Treatment_Group), sep = "||")
cls_lut <- tapply(trimws(as.character(bio$Challenge_Type)), bk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1] })
cls <- unname(cls_lut[paste(r2s[runs], r2t[runs], sep = "||")])
study_all <- unname(r2s[runs])
log_msg("  %d runs; class counts:", length(runs))
print(table(cls, useNA = "ifany"))

auc_of <- function(p, y) {
  p1 <- p[y == 1]; p0 <- p[y == 0]
  if (length(p1) < 3 || length(p0) < 3) return(NA_real_)
  r <- rank(c(p1, p0))
  (sum(r[seq_along(p1)]) - length(p1)*(length(p1)+1)/2) / (length(p1)*length(p0))
}

# -----------------------------------------------------------------------------
# one binary analysis: pos = class(es) coded 1, neg = class(es) coded 0
# -----------------------------------------------------------------------------
run_binary <- function(pos_set, neg_set, label) {
  keep <- cls %in% c(pos_set, neg_set)
  X <- t(M[, keep, drop = FALSE])
  y <- as.integer(cls[keep] %in% pos_set)
  s <- study_all[keep]

  tb <- table(s, y)
  if (!all(c("0","1") %in% colnames(tb))) {
    log_msg("  %s: one class absent entirely - skipped", label); return(NULL) }
  both <- rownames(tb)[tb[,"0"] > 0 & tb[,"1"] > 0]
  ok <- s %in% both
  X <- X[ok,,drop=FALSE]; y <- y[ok]; s <- s[ok]
  if (length(both) < 4) {
    log_msg("  %s: only %d studies contain both classes - too few", label, length(both)); return(NULL) }

  log_msg("")
  log_msg("=== %s ===", label)
  log_msg("  %d runs (%d positive, %d negative) across %d studies with both classes",
          nrow(X), sum(y==1), sum(y==0), length(both))

  Xc <- X
  for (st in unique(s)) { i <- s == st
    Xc[i,] <- sweep(X[i,,drop=FALSE], 2, colMeans(X[i,,drop=FALSE]), "-") }
  Xc <- Xc[, apply(Xc, 2, function(v) stats::sd(v) > 0), drop = FALSE]

  studies <- unique(s); ntake <- max(3, floor(length(studies)*SUBFRAC))
  cnt <- setNames(numeric(ncol(Xc)), colnames(Xc)); nok <- 0
  for (b in seq_len(N_BOOT)) {
    tk <- sample(studies, ntake); idx <- which(s %in% tk)
    if (length(idx) < 20 || length(unique(y[idx])) < 2) next
    f <- tryCatch(cv.glmnet(Xc[idx,,drop=FALSE], y[idx], family="binomial",
                            alpha=1, nfolds=5, standardize=TRUE), error=function(e) NULL)
    if (is.null(f)) next
    co <- as.matrix(coef(f, s=f$lambda.1se)); co <- co[rownames(co)!="(Intercept)",,drop=FALSE]
    sel <- rownames(co)[co[,1] != 0]
    cnt[sel] <- cnt[sel] + 1; nok <- nok + 1
  }
  freq <- if (nok>0) cnt/nok else cnt
  freq <- freq[!(names(freq) %in% UNASSIGNABLE)]
  log_msg("  %d of %d resamples usable; genera at %.0f%%: %d",
          nok, N_BOOT, 100*STABILITY_THR, sum(freq >= STABILITY_THR))
  top <- sort(freq, decreasing=TRUE)[1:min(8,length(freq))]
  for (i in seq_along(top)) log_msg("    %-30s %3.0f%%", gsub("^g__","",names(top)[i]), 100*top[i])

  pred <- rep(NA_real_, length(y))
  for (st in unique(s)) {
    tr <- s != st; te <- s == st
    if (length(unique(y[tr])) < 2 || sum(te) < 2) next
    f <- tryCatch(cv.glmnet(Xc[tr,,drop=FALSE], y[tr], family="binomial",
                            alpha=1, nfolds=5, standardize=TRUE), error=function(e) NULL)
    if (is.null(f)) next
    pred[te] <- as.numeric(predict(f, newx=Xc[te,,drop=FALSE], s=f$lambda.1se, type="response"))
  }
  okp <- !is.na(pred)
  a_loso <- auc_of(pred[okp], y[okp])
  ff <- cv.glmnet(Xc, y, family="binomial", alpha=1, nfolds=5, standardize=TRUE)
  a_in <- auc_of(as.numeric(predict(ff, newx=Xc, s=ff$lambda.1se, type="response")), y)
  n1se <- sum(as.matrix(coef(ff, s=ff$lambda.1se))[-1,1] != 0)

  log_msg("  in-sample AUC = %.3f (%d genera)", a_in, n1se)
  log_msg("  LOSO AUC      = %.3f   <- the honest figure", a_loso)
  if (!is.na(a_in) && !is.na(a_loso) && a_in - a_loso > 0.15)
    log_msg("  *** large gap: fitting study-specific structure, not transferable signal.")

  list(label=label, runs=nrow(Xc), studies=length(both), pos=sum(y==1), neg=sum(y==0),
       max_stab=if(length(freq)) max(freq) else NA, n70=sum(freq>=STABILITY_THR),
       auc_in=a_in, auc_loso=a_loso, n1se=n1se, top=top)
}

log_msg("[2/3] Grouping B - healthy (control + additive) vs challenged")
B <- run_binary(CHALLENGE, c(HEALTHY, ADDITIVE), "B  healthy vs challenged")

log_msg("")
log_msg("[3/3] Grouping C - three classes, one-vs-rest")
C1 <- run_binary(CHALLENGE, c(HEALTHY, ADDITIVE), "C1 challenge vs rest")
C2 <- run_binary(ADDITIVE,  c(HEALTHY, CHALLENGE), "C2 additive vs rest")
C3 <- run_binary(HEALTHY,   c(ADDITIVE, CHALLENGE), "C3 control vs rest")

log_msg("")
log_msg("=========================================================")
log_msg("  %-26s %6s %8s %9s %9s %7s", "grouping", "runs", "studies", "in-samp", "LOSO", ">=70%")
res <- Filter(Negate(is.null), list(B, C2, C3))
for (r in res)
  log_msg("  %-26s %6d %8d %9.3f %9.3f %7d",
          r$label, r$runs, r$studies, r$auc_in, r$auc_loso, r$n70)
log_msg("")
log_msg("  (C1 is identical to B by construction and is not repeated above.)")
# Read 27's figures from its own output rather than repeating them here; a
# hardcoded copy goes stale the moment the arm classes change.
# 28 sources 00_r_config.R, which defines R_OUT_DIR; DIR_GDI is a 50_setup.R
# variable and does not exist here.
ref_file <- file.path(R_OUT_DIR, "gdi", "runlevel_gdi.txt")
ref <- if (file.exists(ref_file)) readLines(ref_file, warn = FALSE) else character(0)
grab <- function(pat) {
  hit <- grep(pat, ref, perl = TRUE, value = TRUE)
  if (!length(hit)) return(NA_character_)
  sub(pat, "\\1", hit[1], perl = TRUE)
}
# Patterns match runlevel_gdi.txt as written by 27, not its console log:
#   Runs: 268 (134 control, 134 challenged)
#   Studies contributing both classes: 13
#   In-sample AUC: 0.775 (12 genera at lambda.1se)
#   LOSO AUC: 0.471   <- the honest figure
r_runs <- grab("^Runs: ([0-9]+).*$")
r_stud <- grab("^Studies contributing both classes: ([0-9]+).*$")
r_in   <- grab("^In-sample AUC: ([0-9.]+).*$")
r_loso <- grab("^LOSO AUC: ([0-9.]+).*$")
if (!any(is.na(c(r_runs, r_stud, r_in, r_loso)))) {
  log_msg("  27_runlevel_gdi.R, control vs challenge only: %s runs, %s studies,",
          r_runs, r_stud)
  log_msg("  in-sample %s, LOSO %s.", r_in, r_loso)
} else {
  log_msg("  For the control-vs-challenge-only model see runlevel_gdi.txt.")
}

out <- c("RUN-LEVEL GDI UNDER ALTERNATIVE TREATMENT GROUPINGS", strrep("=",74), "",
  "No GHI used. Target is the treatment label under three alternative groupings, shown against the reference fit below.", "",
  sprintf("  %-26s %6s %8s %9s %9s %7s", "grouping","runs","studies","in-samp","LOSO",">=70%"))
for (r in res)
  out <- c(out, sprintf("  %-26s %6d %8d %9.3f %9.3f %7d",
           r$label, r$runs, r$studies, r$auc_in, r$auc_loso, r$n70))
out <- c(out, "",
  "Reference - control vs challenge only (27_runlevel_gdi.R):",
  if (!any(is.na(c(r_runs, r_stud, r_in, r_loso))))
    sprintf("  %s runs, %s studies, in-sample AUC %s, LOSO AUC %s",
            r_runs, r_stud, r_in, r_loso)
  else "  see runlevel_gdi.txt", "",
  "LOSO AUC is the honest figure. A large in-sample/LOSO gap indicates the",
  "model fits study-specific structure rather than a transferable signature.")
for (r in res) {
  out <- c(out, "", sprintf("%s - strongest candidates:", r$label))
  for (i in seq_along(r$top))
    out <- c(out, sprintf("  %-30s %3.0f%%", gsub("^g__","",names(r$top)[i]), 100*r$top[i]))
}
writeLines(out, file.path(R_OUT_DIR, "gdi", "runlevel_groupings.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "runlevel_groupings.txt"))
log_msg("=== DONE ===")
