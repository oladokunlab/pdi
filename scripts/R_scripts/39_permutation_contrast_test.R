#!/usr/bin/env Rscript
# =============================================================================
# 39_permutation_contrast_test.R
#
# Fully nonparametric significance test for the 13 genera flagged by
# 37_per_taxon_contrast.R - the test owed since Jesus's original question
# and not yet built.
#
# RUN:  Rscript 39_permutation_contrast_test.R
#
# -----------------------------------------------------------------------------
# WHAT THIS CLOSES
# -----------------------------------------------------------------------------
# 37's pooled p-values are correct but only half nonparametric. The per-study
# statistic (a mean difference) makes no distributional assumption, but the
# POOLING step does: each per-study difference is combined by a random-effects
# formula and the final p-value comes from treating the combined result as
# approximately normal. That normal-theory step is what Jesus was asking
# about, and it was never actually tested against a distribution-free
# alternative.
#
# This script replaces ONLY that pooling/significance step with a
# permutation test, while keeping everything else identical: same per-study
# difference definition, same random-effects weighting (weight = n_ctrl +
# n_trt), same three contrasts. Only the source of the p-value changes.
#
# -----------------------------------------------------------------------------
# HOW THE PERMUTATION WORKS
# -----------------------------------------------------------------------------
# For a given genus and contrast, each contributing study supplies some
# control arms and some treated arms. The real pooled estimate is computed
# from the true control/treated labels. To build a null distribution: within
# each study, randomly reassign which arms are "control" and which are
# "treated" (keeping the same NUMBER of each, since that only shuffles
# labels, not the data), recompute that study's difference under the fake
# labels, then pool across studies exactly as before. Repeat this thousands
# of times. The empirical p-value is the fraction of these relabelled,
# pooled results at least as extreme as the real one.
#
# This makes no assumption about the shape of any distribution anywhere.
# Its only real limitation is granularity: a study with one control arm and
# one treated arm has only two possible relabellings, so it cannot
# contribute much resolution to the null distribution on its own - the
# permutation is still valid, it is just coarse-grained for small studies,
# which is a property of the data, not a flaw in the method.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

CONTROL   <- "Control"
ADDITIVE  <- "Non_Disruptive_Additive"
CHALLENGE <- c("Pathogen_Challenge", "Antibiotic_Challenge", "Dietary_Feed_Challenge")
N_PERM <- 5000
set.seed(42)

log_msg("[1/4] Loading the same candidate genera and arm data as script 37")
clr <- as.data.frame(fread(file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv")))
gname_col <- names(clr)[1]
genera <- clr[[gname_col]]
X <- as.matrix(clr[, -1, drop = FALSE]); rownames(X) <- genera
runs <- colnames(X)

md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2study <- setNames(trimws(as.character(md$Study_ID)), rid)
r2arm   <- setNames(trimws(as.character(md$Treatment_Group)), rid)

ghi <- as.data.frame(fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")))
arm2cls <- setNames(trimws(as.character(ghi$Challenge_Type)),
                     paste(trimws(as.character(ghi$Study_ID)),
                           trimws(as.character(ghi$Treatment_Group)), sep = "||"))

common <- intersect(runs, names(r2study))
arm_key <- paste(r2study[common], r2arm[common], sep = "||")
uniq_arms <- unique(arm_key)
A <- matrix(NA_real_, nrow = nrow(X), ncol = length(uniq_arms),
            dimnames = list(rownames(X), uniq_arms))
for (a in uniq_arms) A[, a] <- rowMeans(X[, common[arm_key == a], drop = FALSE], na.rm = TRUE)
arm_study <- sub("\\|\\|.*$", "", uniq_arms)
arm_cls   <- unname(arm2cls[uniq_arms])

candidates <- c("g__UCG-005","g__Phascolarctobacterium","g__Escherichia-Shigella","g__CHKCI001",
                "g__Merdibacter","g__Paludicola","g__RF39","g__Streptococcus",
                "g__Clostridia_vadinBB60_group","g__Tyzzerella","g__Ligilactobacillus",
                "g__Lactobacillus","g__Anaerostipes")
candidates <- intersect(candidates, rownames(A))
log_msg("  %d candidate genera loaded", length(candidates))

pool_re <- function(diffs, weights) {
  k <- length(diffs)
  if (k < 2) return(NA_real_)
  mu_f <- sum(weights * diffs) / sum(weights)
  Q <- sum(weights * (diffs - mu_f)^2)
  C <- sum(weights) - sum(weights^2) / sum(weights)
  tau2 <- max(0, (Q - (k - 1)) / C)
  wi <- 1 / (1 / weights + tau2)
  sum(wi * diffs) / sum(wi)
}

# returns, per study contributing to this contrast: the pooled arm values
# (control + treated together), and how many belong to each group - the
# raw material for both the real statistic and every permutation.
get_study_pools <- function(genus, treated_classes) {
  vals <- A[genus, ]
  out <- list()
  for (s in unique(arm_study)) {
    idx <- arm_study == s
    ctrl_idx <- idx & arm_cls == CONTROL
    trt_idx  <- idx & arm_cls %in% treated_classes
    ctrl_idx[is.na(ctrl_idx)] <- FALSE; trt_idx[is.na(trt_idx)] <- FALSE
    v_c <- vals[ctrl_idx & !is.na(vals)]; v_t <- vals[trt_idx & !is.na(vals)]
    if (length(v_c) >= 1 && length(v_t) >= 1) out[[s]] <- list(pool = c(v_c, v_t), nc = length(v_c), nt = length(v_t))
  }
  out
}

real_stat <- function(pools) {
  diffs <- vapply(pools, function(p) mean(p$pool[(p$nc+1):(p$nc+p$nt)]) - mean(p$pool[1:p$nc]), numeric(1))
  w <- vapply(pools, function(p) p$nc + p$nt, numeric(1))
  pool_re(diffs, w)
}

perm_stat <- function(pools) {
  diffs <- vapply(pools, function(p) {
    perm <- sample(p$pool)
    mean(perm[(p$nc+1):(p$nc+p$nt)]) - mean(perm[1:p$nc])
  }, numeric(1))
  w <- vapply(pools, function(p) p$nc + p$nt, numeric(1))
  pool_re(diffs, w)
}

run_contrast <- function(label, treated_classes) {
  log_msg("  --- %s ---", label)
  log_msg("  %-28s %5s %10s %14s %14s", "genus", "k", "obs.est", "parametric p", "permutation p")
  res <- list()
  for (g in candidates) {
    pools <- get_study_pools(g, treated_classes)
    if (length(pools) < 2) next
    obs <- real_stat(pools)
    null_dist <- replicate(N_PERM, perm_stat(pools))
    null_dist <- null_dist[!is.na(null_dist)]
    p_perm <- (1 + sum(abs(null_dist) >= abs(obs))) / (length(null_dist) + 1)
    res[[g]] <- list(k = length(pools), est = obs, p_perm = p_perm)
    log_msg("  %-28s %5d %10.3f %14s %14.4f%s",
            g, length(pools), obs, "(see 37's table)", p_perm,
            if (p_perm <= 0.05) "  *" else "")
  }
  res
}

log_msg("")
log_msg("[2/4] ADDITIVE vs control, permutation test, %d permutations per genus", N_PERM)
res_add <- run_contrast("ADDITIVE", ADDITIVE)

log_msg("")
log_msg("[3/4] CHALLENGE vs control")
res_chal <- run_contrast("CHALLENGE", CHALLENGE)

log_msg("")
log_msg("[4/4] COMBINED")
res_comb <- run_contrast("COMBINED", c(ADDITIVE, CHALLENGE))

to_df <- function(res, label) {
  if (!length(res)) return(NULL)
  data.frame(genus = names(res), contrast = label,
             k = vapply(res, function(x) x$k, numeric(1)),
             estimate = vapply(res, function(x) x$est, numeric(1)),
             permutation_p = vapply(res, function(x) x$p_perm, numeric(1)))
}
out <- do.call(rbind, list(to_df(res_add,"ADDITIVE"), to_df(res_chal,"CHALLENGE"), to_df(res_comb,"COMBINED")))
fwrite(out, file.path(R_OUT_DIR, "gdi", "taxon_contrast_permutation.tsv"), sep = "\t")

log_msg("")
log_msg("=== DONE ===")
log_msg("Full table -> %s", file.path(R_OUT_DIR, "gdi", "taxon_contrast_permutation.tsv"))
log_msg("")
log_msg("READING: compare the permutation_p column here against the parametric p in")
log_msg("taxon_contrast_meta.tsv (script 37's output), genus by genus. Close agreement")
log_msg("means the normal-theory approximation was fine for this data. Large")
log_msg("disagreement, especially for genera with few contributing studies, means the")
log_msg("parametric p-value should not be trusted on its own for that genus.")
