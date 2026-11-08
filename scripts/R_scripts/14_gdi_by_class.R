#!/usr/bin/env Rscript
# =============================================================================
# 14_gdi_by_class.R
#
# Build and evaluate the GDI separately within control, additive and challenge
# arms, each against the GHI of the same arms.
#
# RUN:  Rscript 14_gdi_by_class.R     (needs 05 and 06 to have run)
# TIME: 5-15 minutes.
#
# -----------------------------------------------------------------------------
# THE QUESTION
# -----------------------------------------------------------------------------
# Every previous fit pooled all arms and asked whether genus abundance tracks
# the GHI. That assumes one relationship holds across every treatment context.
# It may not. The microbiome of a challenged bird and of a bird on a feed
# additive are perturbed in different ways, and a relationship present inside
# one group can be cancelled by an opposite or absent relationship in another
# when the two are pooled.
#
# Here selection and evaluation are run separately within each class.
#
# -----------------------------------------------------------------------------
# HOW CENTRING IS HANDLED, AND WHY IT MATTERS HERE
# -----------------------------------------------------------------------------
# Most studies contribute a single control arm. Centring within study AFTER
# splitting by class would therefore leave the control group with nothing - a
# single value minus its own mean is zero.
#
# So both abundances and the GHI are centred within study using ALL of that
# study's arms, which is what removes the study baseline, and only then split
# by class. Each class analysis is thus a within-study contrast, and a class
# represented by one arm in a study still contributes its deviation from that
# study's mean.
#
# -----------------------------------------------------------------------------
# WHAT IS REPORTED
# -----------------------------------------------------------------------------
# For each class: how many arms and studies are available, whether stability
# selection yields a panel under the same registered rule used throughout, the
# in-sample fit, and the correlation between the resulting score and the GHI.
#
# Sample sizes per class are small - roughly 35 to 60 arms - so a panel is
# unlikely and an in-sample correlation would be optimistic. The informative
# output is the comparison between classes: if one class shows markedly more
# signal than the others, that is worth pursuing; if all three are flat, the
# pooled null was not caused by mixing contexts.
# =============================================================================
source("00_r_config.R")
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

N_BOOT <- 1000; SUBFRAC <- 0.50; FLOOR <- 0.20; REGISTER <- 0.70
set.seed(42)
UNASSIGNABLE <- c("g__uncultured","g__unclassified","g__uncultured_bacterium",
                  "g__metagenome","g__Incertae_Sedis","g__unidentified")
HEALTHY <- "Control"
CHALLENGE <- c("Pathogen_Challenge","Antibiotic_Challenge","Dietary_Feed_Challenge")
ADDITIVE <- "Non_Disruptive_Additive"

log_msg("[1/3] Loading")
m <- readRDS(file.path(R_OUT_DIR, "gdi", "gdi_model.rds"))
X <- m$X; y <- m$y; study <- m$study; arm <- rownames(X)

gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
gk <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
cls <- unname(tapply(trimws(as.character(gh$Challenge_Type)), gk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1] })[arm])

# Centre within study using ALL arms, then split. See the note above.
Xc <- X; yc <- y
for (s in unique(study)) {
  i <- study == s
  if (sum(i) >= 2) {
    Xc[i, ] <- sweep(X[i, , drop = FALSE], 2, colMeans(X[i, , drop = FALSE]), "-")
    yc[i]   <- y[i] - mean(y[i])
  } else { Xc[i, ] <- NA; yc[i] <- NA }
}
log_msg("  %d arms, %d studies", nrow(X), length(unique(study)))
log_msg("  control %d, additive %d, challenge %d, other/unlabelled %d",
        sum(cls %in% HEALTHY), sum(cls %in% ADDITIVE),
        sum(cls %in% CHALLENGE), sum(is.na(cls)))

get_coefs <- function(f, l) {
  co <- as.matrix(coef(f, s = l)); co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  co[co[,1] != 0, 1]
}

run_class <- function(sel, label) {
  ok <- sel & !is.na(yc) & stats::complete.cases(Xc)
  Xs <- Xc[ok, , drop = FALSE]; ys <- yc[ok]; ss <- study[ok]
  nst <- length(unique(ss))
  log_msg("")
  log_msg("  === %s: %d arms, %d studies ===", label, nrow(Xs), nst)
  if (nrow(Xs) < 20 || nst < 6) {
    log_msg("    too few arms or studies to fit"); return(NULL) }

  st <- unique(ss); nt <- max(3, floor(nst * SUBFRAC))
  cnt <- setNames(numeric(ncol(Xs)), colnames(Xs)); nok <- 0
  for (b in seq_len(N_BOOT)) {
    tk <- sample(st, nt); idx <- which(ss %in% tk)
    if (length(idx) < 10 || length(unique(ys[idx])) < 5) next
    f <- tryCatch(cv.glmnet(Xs[idx, , drop = FALSE], ys[idx], alpha = 1,
                            nfolds = min(5, length(idx)), standardize = TRUE,
                            intercept = FALSE), error = function(e) NULL)
    if (is.null(f)) next
    cnt[names(get_coefs(f, f$lambda.1se))] <- cnt[names(get_coefs(f, f$lambda.1se))] + 1
    nok <- nok + 1
  }
  if (nok < 50) { log_msg("    only %d usable resamples", nok); return(NULL) }
  fr <- cnt / nok; fr <- fr[!(names(fr) %in% UNASSIGNABLE)]

  full <- tryCatch({ cv <- cv.glmnet(Xs, ys, alpha = 1, standardize = TRUE, intercept = FALSE)
                     list(sel = get_coefs(cv, cv$lambda.1se), fit = cv) },
                   error = function(e) NULL)
  panel <- if (is.null(full)) character(0) else
    setdiff(names(fr)[fr >= FLOOR & names(fr) %in% names(full$sel)], UNASSIGNABLE)

  r2 <- NA_real_
  if (!is.null(full)) {
    pr <- as.numeric(predict(full$fit, newx = Xs, s = full$fit$lambda.1se))
    r2 <- if (stats::sd(pr) == 0) 0 else max(0, 1 - sum((ys - pr)^2) / sum((ys - mean(ys))^2))
  }
  top <- sort(fr, decreasing = TRUE)[seq_len(min(6, length(fr)))]
  log_msg("    max stability %.0f%%, %d genera at >=%.0f%%, %d at the registered %.0f%%",
          100*max(fr), sum(fr >= FLOOR), 100*FLOOR, sum(fr >= REGISTER), 100*REGISTER)
  log_msg("    panel under the full rule: %d genera", length(panel))
  log_msg("    in-sample R2 = %.3f", r2)
  log_msg("    strongest:")
  for (i in seq_along(top))
    log_msg("      %-38s %3.0f%%%s", gsub("^g__","",names(top)[i]), 100*top[i],
            if (names(top)[i] %in% panel) "  [panel]" else "")
  list(label = label, n = nrow(Xs), k = nst, max = max(fr),
       n_panel = length(panel), r2 = r2, panel = panel)
}

log_msg("[2/3] Fitting within each class")
res <- list(
  run_class(cls %in% HEALTHY,   "CONTROL arms"),
  run_class(cls %in% ADDITIVE,  "ADDITIVE arms"),
  run_class(cls %in% CHALLENGE, "CHALLENGE arms"),
  run_class(!is.na(cls),        "ALL arms (reference)"))
res <- Filter(Negate(is.null), res)

log_msg("[3/3] Comparison")
log_msg("")
log_msg("  %-22s %6s %8s %10s %8s %7s", "class", "arms", "studies", "max stab", "panel", "R2")
for (z in res)
  log_msg("  %-22s %6d %8d %9.0f%% %8d %7.3f", z$label, z$n, z$k, 100*z$max, z$n_panel, z$r2)

best <- res[[which.max(vapply(res, function(z) z$max, numeric(1)))]]
verdict <- if (best$n_panel > 0 && best$r2 > 0.05) {
  sprintf(paste("%s shows the most signal (max stability %.0f%%, R2 %.3f, %d-genus",
                "panel). Worth validating by leave-one-study-out before it is",
                "treated as a result - the in-sample figure is optimistic on a",
                "subset this size."), best$label, 100*best$max, best$r2, best$n_panel)
} else {
  paste("No class yields a panel with any in-sample fit. Splitting by treatment",
        "context does not recover signal, so the pooled null was not caused by",
        "mixing contexts with different relationships.")
}
log_msg("  %s", verdict)

out <- c("GDI FITTED SEPARATELY WITHIN EACH TREATMENT CLASS", strrep("=",74), "",
  "Abundances and the GHI are centred within study using all arms of that",
  "study, then split by class, so that a class represented by a single arm in a",
  "study still contributes its deviation from that study's mean.", "",
  sprintf("  %-22s %6s %8s %10s %8s %7s", "class","arms","studies","max stab","panel","R2"))
for (z in res)
  out <- c(out, sprintf("  %-22s %6d %8d %9.0f%% %8d %7.3f",
                        z$label, z$n, z$k, 100*z$max, z$n_panel, z$r2))
out <- c(out, "", verdict, "",
  "Selection rule, stability floor and registered threshold are identical to",
  "the primary analysis; only the arms entering each fit differ. Per-class",
  "sample sizes are small, so in-sample figures are optimistic and a panel",
  "would require leave-one-study-out before being reported as a result.")
writeLines(out, file.path(R_OUT_DIR, "gdi", "gdi_by_class.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "gdi_by_class.txt"))
log_msg("=== DONE ===")
