#!/usr/bin/env Rscript
# =============================================================================
# 21_target_comparison_domains.R
#
# Does the microbiome predict Morphology, Biomarker, or Performance
# separately - not just Combined? Both supervised and unsupervised, all four
# targets side by side.
#
# RUN:  Rscript 21_target_comparison_domains.R   (needs 06_lasso_selection.R
#       and 04b_ghi_pca.R with the Performance split to have run)
# TIME: 5-15 minutes (four supervised fits, each with 1000 resamples).
#
# -----------------------------------------------------------------------------
# WHY THIS IS WORTH RUNNING NOW, SPECIFICALLY
# -----------------------------------------------------------------------------
# Every prior GDI attempt targeted GHI_Combined - the average of all ten
# domains. Averaging can hide a real relationship: if the microbiome
# genuinely tracks one facet of gut health strongly and another not at all,
# blending them into one number dilutes the signal that exists.
#
# Performance did not exist as its own target until today (it was folded into
# Biomarker) and is now the strongest, cleanest-validating component of the
# whole GHI - the only one with a formally significant additive-vs-challenge
# difference (p=0.0004). That is new, specific motivation to ask whether the
# microbiome predicts it separately, alongside the same question for
# Morphology and the now-corrected Biomarker.
#
# -----------------------------------------------------------------------------
# SUPERVISED vs UNSUPERVISED, BOTH REPORTED
# -----------------------------------------------------------------------------
# SUPERVISED: LASSO with stability selection, fitted separately against each
# target. Identical rule throughout - within-study centering, study-level 50%
# resampling without replacement, 1000 iterations, registered 70% threshold.
# Only the target column changes between runs.
#
# UNSUPERVISED: the signed control-to-challenge axis (built from genus
# abundances and treatment labels only, no GHI anywhere in its construction)
# is built ONCE, then correlated against each of the four targets in turn.
# This is not four different axes - it is one axis, tested against four
# different outcomes, which is the fair way to ask whether the axis is
# capturing something specific to one facet of health.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

N_BOOT <- 1000; SUBFRAC <- 0.50; FLOOR <- 0.20; REGISTER <- 0.70
N_PERM <- 5000
set.seed(42)
UNASSIGNABLE <- c("g__uncultured","g__unclassified","g__uncultured_bacterium",
                  "g__metagenome","g__Incertae_Sedis","g__unidentified")
HEALTHY <- "Control"
CHALLENGE <- c("Pathogen_Challenge","Antibiotic_Challenge","Dietary_Feed_Challenge")

log_msg("[1/6] Loading genus matrix and all four GHI targets")
m <- readRDS(file.path(R_OUT_DIR, "gdi", "gdi_model.rds"))
X_full <- m$X; study_full <- m$study; arm <- rownames(X_full)

gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
gk <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
TARGETS <- c("Morphology", "Biomarker", "Performance", "Combined")
Y <- sapply(TARGETS, function(t) {
  col <- paste0("GHI_", t, "_raw")
  if (!col %in% names(gh)) return(rep(NA_real_, length(arm)))
  unname(tapply(gh[[col]], gk, mean, na.rm = TRUE)[arm])
})
colnames(Y) <- TARGETS
for (t in TARGETS) log_msg("  %-12s available for %d of %d arms", t, sum(!is.na(Y[,t])), length(arm))

cls <- unname(tapply(trimws(as.character(gh$Challenge_Type)), gk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1] })[arm])

# =============================================================================
# 2 - SUPERVISED: LASSO STABILITY SELECTION, ONE RUN PER TARGET
# =============================================================================
log_msg("[2/6] Supervised: LASSO stability selection per target")

get_coefs <- function(f, l) {
  co <- as.matrix(coef(f, s = l)); co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  co[co[,1] != 0, 1]
}
run_lasso <- function(y) {
  ok <- !is.na(y)
  X <- X_full[ok, , drop = FALSE]; yv <- y[ok]; st <- study_full[ok]
  tb <- table(st); keep <- st %in% names(tb)[tb >= 2]
  X <- X[keep, , drop = FALSE]; yv <- yv[keep]; st <- st[keep]
  if (nrow(X) < 20 || length(unique(st)) < 6)
    return(list(n = nrow(X), k = length(unique(st)), max_stab = NA, n70 = 0, r2 = NA))
  Xc <- X; yc <- yv
  for (s in unique(st)) { i <- st == s
    Xc[i,] <- sweep(X[i,,drop=FALSE], 2, colMeans(X[i,,drop=FALSE]), "-")
    yc[i] <- yv[i] - mean(yv[i]) }
  studies <- unique(st); ntake <- max(3, floor(length(studies) * SUBFRAC))
  cnt <- setNames(numeric(ncol(Xc)), colnames(Xc)); nok <- 0
  for (b in seq_len(N_BOOT)) {
    tk <- sample(studies, ntake); idx <- which(st %in% tk)
    if (length(idx) < 10 || length(unique(yc[idx])) < 5) next
    f <- tryCatch(cv.glmnet(Xc[idx,,drop=FALSE], yc[idx], alpha=1,
                            nfolds=min(5,length(idx)), standardize=TRUE, intercept=FALSE),
                 error=function(e) NULL)
    if (is.null(f)) next
    sel <- get_coefs(f, f$lambda.1se); cnt[names(sel)] <- cnt[names(sel)] + 1; nok <- nok + 1
  }
  freq <- if (nok > 0) cnt / nok else cnt
  freq <- freq[!(names(freq) %in% UNASSIGNABLE)]
  r2 <- tryCatch({ cvf <- cv.glmnet(Xc, yc, alpha=1, standardize=TRUE, intercept=FALSE)
                   pr <- as.numeric(predict(cvf, newx=Xc, s=cvf$lambda.1se))
                   if (sd(pr)==0) 0 else max(0, 1 - sum((yc-pr)^2)/sum((yc-mean(yc))^2)) },
                error=function(e) NA_real_)
  list(n = nrow(Xc), k = length(unique(st)), max_stab = if(length(freq)) max(freq) else NA,
       n70 = sum(freq >= REGISTER), r2 = r2,
       top = sort(freq, decreasing=TRUE)[1:min(5,length(freq))])
}

sup_res <- list()
for (t in TARGETS) { log_msg("  fitting: %s", t); sup_res[[t]] <- run_lasso(Y[,t]) }

log_msg("")
log_msg("  %-12s %6s %7s %10s %6s %8s", "target", "arms", "studies", "max stab", ">=70%", "R2")
for (t in TARGETS) { z <- sup_res[[t]]
  log_msg("  %-12s %6d %7d %9s %6d %8s", t, z$n, z$k,
          ifelse(is.na(z$max_stab), "n/a", sprintf("%.0f%%", 100*z$max_stab)),
          z$n70, ifelse(is.na(z$r2), "n/a", sprintf("%.3f", z$r2))) }
log_msg("")
for (t in TARGETS) { z <- sup_res[[t]]
  if (is.null(z$top) || !length(z$top)) next
  log_msg("  %s strongest: %s", t, paste(sprintf("%s=%.0f%%", gsub("^g__","",names(z$top)), 100*z$top), collapse=", ")) }

# =============================================================================
# 3 - UNSUPERVISED: ONE SIGNED AXIS, BUILT ONCE
# =============================================================================
log_msg("[3/6] Unsupervised: building the signed control->challenge axis (no GHI used)")

isC <- cls %in% HEALTHY; isX <- cls %in% CHALLENGE
donors <- unique(study_full)[vapply(unique(study_full), function(s) {
  i <- study_full == s; any(isC[i]) && any(isX[i]) }, logical(1))]
log_msg("  %d studies contribute both a control and a challenged arm", length(donors))

dir_from <- function(s) {
  i <- which(study_full == s); ic <- i[isC[i]]; ix <- i[isX[i]]
  if (!length(ic) || !length(ix)) return(NULL)
  d <- colMeans(X_full[ix,,drop=FALSE]) - colMeans(X_full[ic,,drop=FALSE])
  L <- sqrt(sum(d^2)); if (!is.finite(L) || L==0) return(NULL); d/L
}
axis_score <- rep(NA_real_, length(arm))
for (s in unique(study_full)) {
  others <- setdiff(donors, s)
  dirs <- Filter(Negate(is.null), lapply(others, dir_from))
  if (!length(dirs)) next
  v <- Reduce(`+`, dirs) / length(dirs); L <- sqrt(sum(v^2))
  if (!is.finite(L) || L == 0) next
  i <- which(study_full == s); if (length(i) < 2) next
  Xi <- sweep(X_full[i,,drop=FALSE], 2, colMeans(X_full[i,,drop=FALSE]), "-")
  axis_score[i] <- as.numeric(Xi %*% (v/L))
}
log_msg("  axis scored for %d of %d arms", sum(!is.na(axis_score)), length(arm))

# =============================================================================
# 4 - CORRELATE THE SAME AXIS AGAINST ALL FOUR TARGETS
# =============================================================================
log_msg("[4/6] Correlating the axis against each target")
perm_p <- function(x, y, g, n = N_PERM) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; g <- g[ok]
  if (length(x) < 10) return(c(NA_real_, NA_real_, length(x)))
  obs <- suppressWarnings(cor(x, y, method="spearman"))
  idx <- split(seq_along(x), g); cnt <- 0L
  for (i in seq_len(n)) { yp <- y
    for (j in idx) if (length(j) > 1) yp[j] <- y[sample(j)]
    if (suppressWarnings(cor(x, yp, method="spearman")) <= obs) cnt <- cnt + 1L }
  c(obs, (cnt+1)/(n+1), length(x))
}
within_center <- function(x, y, g) {
  for (s in unique(g)) { i <- g==s & !is.na(x) & !is.na(y)
    if (sum(i) >= 2) { x[i] <- x[i]-mean(x[i]); y[i] <- y[i]-mean(y[i]) } else { x[g==s] <- NA; y[g==s] <- NA } }
  list(x=x, y=y)
}

unsup_res <- list()
for (t in TARGETS) {
  cc <- within_center(axis_score, Y[,t], study_full)
  r <- perm_p(cc$x, cc$y, study_full)
  unsup_res[[t]] <- r
}
log_msg("")
log_msg("  %-12s %10s %10s %6s", "target", "rho", "perm p", "n")
for (t in TARGETS) { r <- unsup_res[[t]]
  log_msg("  %-12s %+10.3f %10s %6d", t, r[1], ifelse(is.na(r[2]),"n/a",sprintf("%.4f",r[2])), r[3]) }

# =============================================================================
# 5 - READING
# =============================================================================
log_msg("[5/6] Reading")
best_sup <- names(which.max(sapply(TARGETS, function(t) ifelse(is.na(sup_res[[t]]$max_stab),0,sup_res[[t]]$max_stab))))
best_unsup <- names(which.min(sapply(TARGETS, function(t) ifelse(is.na(unsup_res[[t]][2]),1,unsup_res[[t]][2]))))
log_msg("  Supervised: %s shows the highest stability (%.0f%%)", best_sup, 100*sup_res[[best_sup]]$max_stab)
log_msg("  Unsupervised: %s shows the strongest axis correlation (p=%.4f)", best_unsup, unsup_res[[best_unsup]][2])
log_msg("  Neither being a formal 'winner' matters unless it clears the pre-registered")
log_msg("  70%% threshold (supervised) or p<0.05 (unsupervised) - check the tables above,")
log_msg("  not just which target ranked highest among four.")

# =============================================================================
# 6 - WRITE
# =============================================================================
log_msg("[6/6] Writing")
out <- c("TARGET COMPARISON: MORPHOLOGY / BIOMARKER / PERFORMANCE / COMBINED",
  strrep("=", 74), "",
  "SUPERVISED (LASSO stability selection, identical rule per target):", "",
  sprintf("  %-12s %6s %7s %10s %6s %8s", "target","arms","studies","max stab",">=70%","R2"))
for (t in TARGETS) { z <- sup_res[[t]]
  out <- c(out, sprintf("  %-12s %6d %7d %9s %6d %8s", t, z$n, z$k,
           ifelse(is.na(z$max_stab),"n/a",sprintf("%.0f%%",100*z$max_stab)), z$n70,
           ifelse(is.na(z$r2),"n/a",sprintf("%.3f",z$r2)))) }
out <- c(out, "", "UNSUPERVISED (one signed control->challenge axis, no GHI in its construction):", "",
  sprintf("  %-12s %10s %10s %6s", "target","rho","perm p","n"))
for (t in TARGETS) { r <- unsup_res[[t]]
  out <- c(out, sprintf("  %-12s %+10.3f %10s %6d", t, r[1],
           ifelse(is.na(r[2]),"n/a",sprintf("%.4f",r[2])), r[3])) }
writeLines(out, file.path(R_OUT_DIR, "gdi", "target_comparison_domains.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "target_comparison_domains.txt"))
log_msg("=== DONE ===")
