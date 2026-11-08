#!/usr/bin/env Rscript
# =============================================================================
# 20_control_vs_treated.R
#
# Control vs everything else (additive and challenge pooled into one
# "treated" group). A genuinely new configuration, not yet tested.
#
# RUN:  Rscript 20_control_vs_treated.R    (needs 06_lasso_selection.R to have run)
#
# -----------------------------------------------------------------------------
# WHY THIS IS DIFFERENT FROM EVERYTHING ELSE TESTED
# -----------------------------------------------------------------------------
# Every earlier configuration split arms by TYPE of intervention (control vs
# additive, or control vs challenge, separately). This pools any intervention
# together and asks the simplest possible question: does anything happening
# to the bird show up in the microbiome at all, regardless of what it is.
#
# Same selection rule, same 70% registered threshold, same resampling as
# every other configuration in this project. Reported whatever it shows.
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
set.seed(42)
UNASSIGNABLE <- c("g__uncultured","g__unclassified","g__uncultured_bacterium",
                  "g__metagenome","g__Incertae_Sedis","g__unidentified")
HEALTHY <- "Control"

log_msg("[1/4] Loading and building the two-group label")
m <- readRDS(file.path(R_OUT_DIR, "gdi", "gdi_model.rds"))
X <- m$X; y <- m$y; study <- m$study; arm <- rownames(X)

gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
gk <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
cls <- unname(tapply(trimws(as.character(gh$Challenge_Type)), gk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1] })[arm])

# TREATED = anything that is not Control and not the excluded categories
group <- ifelse(cls %in% HEALTHY, "Control",
         ifelse(cls %in% c("Challenged_And_Treated","Phenotype_Comparison") | is.na(cls), NA, "Treated"))
keep <- !is.na(group)
tb <- table(study[keep]); keep <- keep & study %in% names(tb)[tb >= 2]
Xk <- X[keep, , drop = FALSE]; yk <- y[keep]; sk <- study[keep]; gk2 <- group[keep]
log_msg("  %d arms (%d control, %d treated) across %d studies",
        nrow(Xk), sum(gk2 == "Control"), sum(gk2 == "Treated"), length(unique(sk)))

Xc <- Xk; yc <- yk
for (s in unique(sk)) { i <- sk == s
  Xc[i, ] <- sweep(Xk[i, , drop = FALSE], 2, colMeans(Xk[i, , drop = FALSE]), "-")
  yc[i] <- yk[i] - mean(yk[i]) }

log_msg("[2/4] Stability selection, %d resamples", N_BOOT)
st <- unique(sk); ntake <- max(3, floor(length(st) * SUBFRAC))
cnt <- setNames(numeric(ncol(Xc)), colnames(Xc)); nok <- 0
for (b in seq_len(N_BOOT)) {
  tkeep <- sample(st, ntake); idx <- which(sk %in% tkeep)
  if (length(idx) < 10 || length(unique(yc[idx])) < 5) next
  f <- tryCatch(cv.glmnet(Xc[idx, , drop = FALSE], yc[idx], alpha = 1,
                          nfolds = min(5, length(idx)), standardize = TRUE, intercept = FALSE),
               error = function(e) NULL)
  if (is.null(f)) next
  co <- as.matrix(coef(f, s = f$lambda.1se)); co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  sel <- rownames(co)[co[,1] != 0]
  cnt[sel] <- cnt[sel] + 1; nok <- nok + 1
}
freq <- cnt / nok; freq <- freq[!(names(freq) %in% UNASSIGNABLE)]
log_msg("  %d of %d resamples usable", nok, N_BOOT)
log_msg("  genera at the registered %.0f%% threshold: %d", 100*REGISTER, sum(freq >= REGISTER))
top <- sort(freq, decreasing = TRUE)[1:min(10, length(freq))]
log_msg("  strongest:")
for (i in seq_along(top)) log_msg("    %-32s %3.0f%%", gsub("^g__","",names(top)[i]), 100*top[i])

log_msg("[3/4] In-sample fit")
full <- tryCatch({ f <- cv.glmnet(Xc, yc, alpha=1, standardize=TRUE, intercept=FALSE)
                   pr <- as.numeric(predict(f, newx=Xc, s=f$lambda.1se))
                   r2 <- if (sd(pr)==0) 0 else max(0, 1 - sum((yc-pr)^2)/sum((yc-mean(yc))^2))
                   list(r2=r2) }, error=function(e) list(r2=NA_real_))
log_msg("  R2 = %.3f", full$r2)

log_msg("[4/4] AUC: treated vs control, on the strongest single genus")
if (length(top) > 0) {
  g1 <- names(top)[1]
  s1 <- Xk[gk2=="Treated", g1]; s0 <- Xk[gk2=="Control", g1]
  r <- rank(c(s1,s0))
  auc <- (sum(r[seq_along(s1)]) - length(s1)*(length(s1)+1)/2) / (length(s1)*length(s0))
  log_msg("  %s alone: AUC = %.3f", gsub("^g__","",g1), auc)
}

out <- c("CONTROL VS EVERYTHING ELSE (ADDITIVE + CHALLENGE POOLED)", strrep("=", 74), "",
  sprintf("Arms: %d (%d control, %d treated), %d studies", nrow(Xk), sum(gk2=="Control"),
          sum(gk2=="Treated"), length(unique(sk))),
  sprintf("Genera at registered 70%% threshold: %d", sum(freq >= REGISTER)),
  sprintf("In-sample R2: %.3f", full$r2), "",
  "Strongest candidates:")
for (i in seq_along(top)) out <- c(out, sprintf("  %-32s %3.0f%%", gsub("^g__","",names(top)[i]), 100*top[i]))
if (length(top) > 0) out <- c(out, "", sprintf("AUC (%s alone), treated vs control: %.3f", gsub("^g__","",names(top)[1]), auc))
writeLines(out, file.path(R_OUT_DIR, "gdi", "control_vs_treated.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "control_vs_treated.txt"))
log_msg("=== DONE ===")
