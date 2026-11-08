# =============================================================================
# 66_variance_components.R
#
# Between-study variance in the GHI and what it does to inference.
#
#   1. Variance components of the unchallenged control groups: a random-
#      intercept model (REML) with study as the grouping factor. Reports the
#      between- and within-study variances, the study share (ICC), the mean
#      number of control groups per study, and the design effect
#      1 + (m - 1) * ICC.
#   2. Paired within-study contrast: for every study holding both an
#      unchallenged control and a challenged group, the mean control score
#      minus the mean challenged score, tested with a one-sample t-test.
#   3. Pooled group-level contrast: OLS of the score on a control indicator,
#      with a naive standard error and a CR1 cluster-robust standard error
#      clustered on study (the small-sample correction used by Stata and by
#      sandwich::vcovCL(type = "HC1")).
#
# The antibiotic-challenge class holds a single group in a single study and
# is withheld from the challenged set, as in the Results.
#
# Input : r_analysis/ghi/ghi_scores_pca.tsv
# Output: r_analysis/verification/66_variance_components.txt
#
# Usage : export PDI_ROOT=/scratch/user/$USER/pdi && Rscript 66_variance_components.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

PDI_ROOT <- Sys.getenv("PDI_ROOT",
                       file.path("/scratch/user", Sys.getenv("USER"), "pdi"))
IN  <- file.path(PDI_ROOT, "r_analysis", "ghi", "ghi_scores_pca.tsv")
OUT <- file.path(PDI_ROOT, "r_analysis", "verification", "66_variance_components.txt")
stopifnot(file.exists(IN))

CONTROL    <- "Unchallenged control"
CHALLENGED <- c("Infectious challenge", "Toxin challenge", "Dietary challenge",
                "Environmental challenge")          # antibiotic withheld (1 group)
SCORE      <- "GHI_Combined_raw"

g <- fread(IN)
g <- g[!is.na(get(SCORE))]

out <- character()
say <- function(...) out <<- c(out, sprintf(...))

say("=== Variance components and clustered inference for the GHI ===")
say("Input: %s", sub(PDI_ROOT, "$PDI_ROOT", IN, fixed = TRUE))
say("Scored groups: %d in %d studies", nrow(g), uniqueN(g$Study_ID))
say("")

# --- 1. variance components of the control groups ---------------------------
ctl <- g[Intervention_Class == CONTROL]
fit <- lmer(as.formula(paste(SCORE, "~ 1 + (1 | Study_ID)")), data = ctl, REML = TRUE)
vc  <- as.data.frame(VarCorr(fit))
v_between <- vc$vcov[vc$grp == "Study_ID"]
v_within  <- vc$vcov[vc$grp == "Residual"]
icc <- v_between / (v_between + v_within)
m   <- nrow(ctl) / uniqueN(ctl$Study_ID)
de  <- 1 + (m - 1) * icc

say("--- 1. Control groups, random intercept for study (REML) ---")
say("control groups          : %d in %d studies", nrow(ctl), uniqueN(ctl$Study_ID))
say("control mean (SD)       : %.3f (%.3f)", mean(ctl[[SCORE]]), sd(ctl[[SCORE]]))
say("between-study variance  : %.3f", v_between)
say("within-study variance   : %.3f", v_within)
say("study share (ICC)       : %.1f%%", 100 * icc)
say("control groups per study: %.2f", m)
say("design effect           : %.2f", de)
say("")

# --- 2. paired within-study contrast ---------------------------------------
per_study <- g[, .(ctl = mean(get(SCORE)[Intervention_Class == CONTROL]),
                   chl = mean(get(SCORE)[Intervention_Class %in% CHALLENGED])),
               by = Study_ID][is.finite(ctl) & is.finite(chl)]
per_study[, diff := ctl - chl]
tt <- t.test(per_study$diff)
say("--- 2. Paired within-study contrast, control minus challenged ---")
say("studies                 : %d", nrow(per_study))
say("mean difference         : %+.3f", unname(tt$estimate))
say("SE                      : %.3f", unname(tt$stderr))
say("t (df)                  : %.2f (%d)", unname(tt$statistic), as.integer(tt$parameter))
say("P (two-sided)           : %.3f", tt$p.value)
say("")

# --- 3. pooled group-level contrast, naive and clustered SE -----------------
d <- g[Intervention_Class %in% c(CONTROL, CHALLENGED)]
d[, ctrl := as.numeric(Intervention_Class == CONTROL)]
ols <- lm(as.formula(paste(SCORE, "~ ctrl")), data = d)
X <- model.matrix(ols); e <- residuals(ols)
n <- nrow(X); k <- ncol(X)
cl <- d$Study_ID; G <- length(unique(cl))
bread <- solve(crossprod(X))
meat  <- Reduce(`+`, lapply(split(seq_len(n), cl), function(i) {
  s <- crossprod(X[i, , drop = FALSE], e[i]); s %*% t(s)
}))
V_cr1 <- (G / (G - 1)) * ((n - 1) / (n - k)) * bread %*% meat %*% bread
est   <- coef(ols)[["ctrl"]]
se_nv <- summary(ols)$coefficients["ctrl", "Std. Error"]
se_cl <- sqrt(V_cr1["ctrl", "ctrl"])
z     <- est / se_cl
say("--- 3. Pooled group-level contrast, control minus challenged ---")
say("groups                  : %d (%d control, %d challenged) in %d studies",
    n, sum(d$ctrl), sum(1 - d$ctrl), G)
say("estimate                : %+.3f", est)
say("naive SE                : %.3f", se_nv)
say("clustered SE (CR1)      : %.3f", se_cl)
say("z, P (clustered)        : %.2f, %.3f", z, 2 * pnorm(-abs(z)))

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
writeLines(out, OUT)
cat(out, sep = "\n")
