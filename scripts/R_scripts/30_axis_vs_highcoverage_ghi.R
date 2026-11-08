#!/usr/bin/env Rscript
# =============================================================================
# 30_axis_vs_highcoverage_ghi.R
#
# Does the unsupervised axis correlate better with a HIGH-COVERAGE GHI than
# with the full one?
#
# RUN:  Rscript 30_axis_vs_highcoverage_ghi.R
# TIME: 2-4 minutes.
#
# -----------------------------------------------------------------------------
# THE REASONING
# -----------------------------------------------------------------------------
# 29_ghi_diagnostics.R showed GHI discrimination rises from AUC 0.653 to 0.852
# once samples with fewer than three measured domains are excluded. Those
# low-coverage samples are essentially single noisy measurements carrying the
# label of a composite index.
#
# Every previous correlation between the microbiome axis and the GHI used the
# full, diluted GHI. If a share of that weak correlation was the TARGET being
# noisy rather than the axis being uninformative, restricting to a
# high-coverage GHI should raise it.
#
# The axis is the right construction to test this with. It has no penalty
# parameter and performs no feature selection, so it has very little to overfit
# with - which is why it is the only construction in this project that reached
# significance out-of-study. A LASSO or random forest at this sample size would
# simply memorise.
#
# -----------------------------------------------------------------------------
# WHAT THIS CANNOT ESTABLISH
# -----------------------------------------------------------------------------
# Restricting coverage costs arms: roughly 141 down to something near 49. That
# is a severe power penalty, and a non-significant result here would be
# genuinely ambiguous - underpowered rather than absent.
#
# So the informative comparison is not significance but EFFECT SIZE. If rho
# rises materially as coverage increases, target noise was part of the problem.
# If rho is flat or falls, it was not, and the weak correlation is a property
# of the microbiome rather than of the GHI.
#
# Both the correlation and the number of arms behind it are reported at every
# coverage level so the trade is visible.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

HEALTHY   <- "Control"
CHALLENGE <- c("Pathogen_Challenge","Antibiotic_Challenge","Dietary_Feed_Challenge")
N_PERM <- 5000
set.seed(42)

# =============================================================================
log_msg("[1/4] Building the unsupervised axis (no GHI used, leave-one-study-out)")
m <- readRDS(file.path(R_OUT_DIR, "gdi", "gdi_model.rds"))
X <- m$X; study <- m$study; arm <- rownames(X)

gh <- as.data.frame(fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")))
gh$Study_ID <- trimws(gh$Study_ID)
gk <- paste(gh$Study_ID, trimws(gh$Treatment_Group), sep = "||")
cls <- unname(tapply(trimws(as.character(gh$Challenge_Type)), gk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1] })[arm])
isC <- cls %in% HEALTHY; isX <- cls %in% CHALLENGE

donors <- unique(study)[vapply(unique(study), function(s) {
  i <- study == s; any(isC[i]) && any(isX[i]) }, logical(1))]
log_msg("  %d donor studies (contain both a control and a challenged arm)", length(donors))

dir_from <- function(s) {
  i <- which(study == s); ic <- i[isC[i]]; ix <- i[isX[i]]
  if (!length(ic) || !length(ix)) return(NULL)
  d <- colMeans(X[ix,,drop=FALSE]) - colMeans(X[ic,,drop=FALSE])
  L <- sqrt(sum(d^2)); if (!is.finite(L) || L == 0) return(NULL); d/L
}
axis <- rep(NA_real_, length(arm))
for (s in unique(study)) {
  dirs <- Filter(Negate(is.null), lapply(setdiff(donors, s), dir_from))
  if (!length(dirs)) next
  v <- Reduce(`+`, dirs)/length(dirs); L <- sqrt(sum(v^2))
  if (!is.finite(L) || L == 0) next
  i <- which(study == s); if (length(i) < 2) next
  Xi <- sweep(X[i,,drop=FALSE], 2, colMeans(X[i,,drop=FALSE]), "-")
  axis[i] <- as.numeric(Xi %*% (v/L))
}
log_msg("  axis scored for %d of %d arms", sum(!is.na(axis)), length(arm))

# =============================================================================
log_msg("[2/4] Building GHI targets at each coverage level")
dom_cols <- grep("^GHI_dom_", names(gh), value = TRUE)
gh$ndom <- rowSums(!is.na(gh[, dom_cols, drop = FALSE]))

# collapse the per-sample GHI to arm level at each coverage threshold
arm_ghi <- function(min_dom, target) {
  sel <- gh$ndom >= min_dom & is.finite(gh[[target]])
  if (!any(sel)) return(setNames(rep(NA_real_, length(arm)), arm))
  k <- paste(gh$Study_ID[sel], trimws(gh$Treatment_Group[sel]), sep = "||")
  v <- tapply(gh[[target]][sel], k, mean)
  unname(v[arm])
}

# =============================================================================
log_msg("[3/4] Correlating the axis against each target")
within_center <- function(x, y, g) {
  for (s in unique(g)) {
    i <- g == s & !is.na(x) & !is.na(y)
    if (sum(i) >= 2) { x[i] <- x[i] - mean(x[i]); y[i] <- y[i] - mean(y[i]) }
    else { x[g == s] <- NA; y[g == s] <- NA }
  }
  list(x = x, y = y)
}
perm_p <- function(x, y, g, n = N_PERM) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; g <- g[ok]
  if (length(x) < 10) return(c(NA_real_, NA_real_, length(x), length(unique(g))))
  obs <- suppressWarnings(cor(x, y, method = "spearman"))
  idx <- split(seq_along(x), g); cnt <- 0L
  for (i in seq_len(n)) { yp <- y
    for (j in idx) if (length(j) > 1) yp[j] <- y[sample(j)]
    if (suppressWarnings(cor(x, yp, method = "spearman")) <= obs) cnt <- cnt + 1L }
  c(obs, (cnt+1)/(n+1), length(x), length(unique(g)))
}

TARGETS <- c(Combined = "GHI_Combined_raw", Biomarker = "GHI_Biomarker_raw")
res <- list()
for (tn in names(TARGETS)) {
  log_msg("")
  log_msg("  --- target: GHI %s ---", tn)
  log_msg("  %-16s %10s %10s %7s %8s", "min domains", "rho", "perm p", "arms", "studies")
  for (k in 1:4) {
    y <- arm_ghi(k, TARGETS[[tn]])
    cc <- within_center(axis, y, study)
    r <- perm_p(cc$x, cc$y, study)
    res[[paste(tn, k)]] <- c(list(target = tn, mindom = k), as.list(r))
    log_msg("  >= %-13d %+10.3f %10s %7d %8d", k, r[1],
            ifelse(is.na(r[2]), "n/a", sprintf("%.4f", r[2])), r[3], r[4])
  }
}

# =============================================================================
log_msg("")
log_msg("[4/4] Reading")
get <- function(tn, k) res[[paste(tn, k)]]
for (tn in names(TARGETS)) {
  a <- get(tn, 1); b <- get(tn, 3)
  if (is.null(a) || is.null(b) || is.na(a[[3]]) || is.na(b[[3]])) next
  d <- abs(b[[3]]) - abs(a[[3]])
  log_msg("  %s: rho %+.3f at all samples (%d arms) -> %+.3f at >=3 domains (%d arms)",
          tn, a[[3]], a[[5]], b[[3]], b[[5]])
  if (d > 0.05)
    log_msg("      correlation STRENGTHENED as target noise fell - part of the weak")
    else if (d < -0.05)
    log_msg("      correlation WEAKENED - the restriction cost more in power than it")
    else
    log_msg("      correlation essentially unchanged - target noise was not the")
}
log_msg("")
log_msg("  Interpretation depends on effect size, not significance: restricting")
log_msg("  coverage costs arms, so a non-significant result at >=3 domains is")
log_msg("  ambiguous rather than negative. A rise in rho would indicate the noisy")
log_msg("  target was part of the problem; a flat or falling rho would indicate")
log_msg("  the weak relationship is a property of the microbiome itself.")

out <- c("UNSUPERVISED AXIS vs GHI AT INCREASING BIOMARKER COVERAGE",
  strrep("=", 74), "",
  "The axis uses genus abundances and treatment labels only, computed",
  "leave-one-study-out. No GHI enters its construction. Only the TARGET",
  "changes across rows: the GHI restricted to samples with at least N domains.", "",
  sprintf("  %-12s %-14s %10s %10s %7s %8s", "target", "min domains", "rho", "perm p", "arms", "studies"))
for (nm in names(res)) { r <- res[[nm]]
  out <- c(out, sprintf("  %-12s >= %-11d %+10.3f %10s %7d %8d",
           r$target, r$mindom, r[[3]],
           ifelse(is.na(r[[4]]), "n/a", sprintf("%.4f", r[[4]])), r[[5]], r[[6]])) }
out <- c(out, "",
  "Restricting coverage costs arms, so significance is not the right test here.",
  "A rise in |rho| with coverage would show that GHI noise contributed to the",
  "weak correlation; a flat or falling |rho| would show it did not.")
writeLines(out, file.path(R_OUT_DIR, "gdi", "axis_vs_coverage.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "axis_vs_coverage.txt"))
log_msg("=== DONE ===")
