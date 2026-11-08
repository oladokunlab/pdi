#!/usr/bin/env Rscript
# =============================================================================
# 12_gdi_three_ways.R
#
# Build three dysbiosis indices from 16S data alone and compare them.
#
# RUN:  Rscript 12_gdi_three_ways.R
# TIME: 2-5 minutes.
#
# -----------------------------------------------------------------------------
# THE THREE CONSTRUCTIONS
# -----------------------------------------------------------------------------
# A. DEVIATION FROM CONTROL
#    Aitchison distance from the study's own control centroid. The literal
#    reading of dysbiosis: departure from a healthy reference.
#    Inherently within-study - it needs that study's own controls - so it
#    cannot be validated out-of-study. It can still be tested against outcomes.
#
# B. SIGNED CHALLENGE AXIS
#    Position along the direction running from control arms to challenged arms.
#    Adds direction to A, so that moving towards challenged scores positive
#    while a beneficial shift scores negative. Built by pooling the direction
#    from OTHER studies, so it is validated out-of-study.
#
# C. FIXED REFERENCE PANEL
#    The canine and feline Dysbiosis Index design: identify the genera that
#    differ most between control and challenged birds, fix them as a panel, and
#    score any sample as a weighted sum. Panel and weights are derived from the
#    training studies only, so this is validated out-of-study.
#
# -----------------------------------------------------------------------------
# HOW EACH IS JUDGED, AND WHY THE TESTS DIFFER
# -----------------------------------------------------------------------------
# Two outcomes are available: the GHI, and the control/challenge label. Neither
# is used to BUILD A or B; C uses the label and is therefore evaluated only
# out-of-study.
#
#   correlation with the GHI, within study
#       The main test. Negative is the working direction: further towards
#       dysbiosis means lower gut health.
#
#   AUC, challenged versus control
#       Meaningless for A, where control arms sit at distance zero by
#       construction and would separate perfectly for a trivial reason.
#       Reported for B and C, and for C only out-of-study.
#
#   AUC, challenged versus additive
#       The fair test for A. Both groups differ from control, so neither is
#       privileged; the question is whether challenge deviates MORE. Reported
#       for all three.
#
# Where a measure is inherently within-study, that is stated rather than
# papered over with a number that cannot mean what it appears to.
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
PANEL_N   <- 7          # following the canine and feline indices
N_PERM    <- 5000
set.seed(42)

# =============================================================================
# 1 - ARM-LEVEL COMMUNITY MATRIX
# =============================================================================
log_msg("[1/5] Loading community data")
clr <- fread(file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv"))
M <- as.matrix(clr[, -1, with = FALSE]); rownames(M) <- clr[[names(clr)[1]]]
md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2s <- setNames(trimws(as.character(md$Study_ID)), rid)
r2t <- setNames(trimws(as.character(md$Treatment_Group)), rid)
M <- M[, colnames(M) %in% names(r2s), drop = FALSE]
ak <- paste(r2s[colnames(M)], r2t[colnames(M)], sep = "||")
A  <- t(vapply(unique(ak), function(a) rowMeans(M[, ak == a, drop = FALSE]), numeric(nrow(M))))
st <- sub("\\|\\|.*$", "", rownames(A))

bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
bk <- paste(trimws(bio$Study_ID), trimws(bio$Treatment_Group), sep = "||")
cls <- unname(tapply(trimws(as.character(bio$Challenge_Type)), bk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1] })[rownames(A)])

gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
gk <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
GHI <- unname(tapply(gh$GHI_Combined_raw, gk, function(z) mean(z, na.rm = TRUE))[rownames(A)])

isC <- cls %in% HEALTHY; isX <- cls %in% CHALLENGE; isA <- cls %in% ADDITIVE
log_msg("  %d arms, %d studies (%d control, %d challenged, %d additive)",
        nrow(A), length(unique(st)), sum(isC), sum(isX), sum(isA))
donors <- names(which(vapply(unique(st), function(s)
  any(isC[st == s]) && any(isX[st == s]), logical(1))))
log_msg("  %d studies contain both a control and a challenged arm", length(donors))

# =============================================================================
# 2 - THE THREE INDICES
# =============================================================================
log_msg("[2/5] Building")
n <- nrow(A)

# ---- A: deviation from the study's own control centroid ---------------------
A_dev <- rep(NA_real_, n)
for (s in unique(st)) {
  i <- which(st == s); ic <- i[isC[i]]
  if (!length(ic)) next
  cen <- colMeans(A[ic, , drop = FALSE])
  A_dev[i] <- sqrt(rowSums(sweep(A[i, , drop = FALSE], 2, cen, "-")^2))
}
log_msg("  A  deviation from control : %d arms (within-study by construction)", sum(!is.na(A_dev)))

# ---- B: signed challenge axis, built from other studies ---------------------
unit_dir <- function(s) {
  i <- which(st == s); ic <- i[isC[i]]; ix <- i[isX[i]]
  if (!length(ic) || !length(ix)) return(NULL)
  d <- colMeans(A[ix, , drop = FALSE]) - colMeans(A[ic, , drop = FALSE])
  L <- sqrt(sum(d^2)); if (!is.finite(L) || L == 0) return(NULL); d / L
}
# Only studies containing both arm types can DEFINE a direction, but once a
# direction exists it can score any arm. A study that contributed nothing to
# building the axis carries no leakage, so the additive-only and control-only
# studies are scored using the axis built from all donors. Donor studies are
# scored leave-one-out, from the other donors.
mean_dir <- function(studies) {
  dirs <- Filter(Negate(is.null), lapply(studies, unit_dir))
  if (!length(dirs)) return(NULL)
  v <- Reduce(`+`, dirs) / length(dirs); L <- sqrt(sum(v^2))
  if (!is.finite(L) || L == 0) return(NULL); v / L
}
B_axis <- rep(NA_real_, n)
ax_all <- mean_dir(donors)
for (s in unique(st)) {
  ax <- if (s %in% donors) mean_dir(setdiff(donors, s)) else ax_all
  if (is.null(ax)) next
  i <- which(st == s); if (length(i) < 2) next
  B_axis[i] <- as.numeric(sweep(A[i, , drop = FALSE], 2, colMeans(A[i, , drop = FALSE]), "-") %*% ax)
}
log_msg("  B  signed challenge axis  : %d arms across %d studies", sum(!is.na(B_axis)),
        length(unique(st[!is.na(B_axis)])))
log_msg("       donor studies scored leave-one-out; all others scored by the")
log_msg("       axis from all %d donors, to which they contributed nothing", length(donors))

# ---- C: fixed reference panel, canine-DI style ------------------------------
# Per genus, the mean within-study difference between challenged and control
# arms across training studies. The largest absolute differences define the
# panel; each genus contributes with the sign of its difference.
panel_from <- function(studies, k = PANEL_N) {
  diffs <- list()
  for (s in studies) {
    i <- which(st == s); ic <- i[isC[i]]; ix <- i[isX[i]]
    if (!length(ic) || !length(ix)) next
    diffs[[s]] <- colMeans(A[ix, , drop = FALSE]) - colMeans(A[ic, , drop = FALSE])
  }
  if (length(diffs) < 3) return(NULL)
  D <- Reduce(`+`, diffs) / length(diffs)
  sel <- names(sort(abs(D), decreasing = TRUE))[seq_len(min(k, length(D)))]
  list(genera = sel, w = sign(D[sel]))
}
# Same logic as B: donors scored leave-one-out, everyone else scored by the
# panel derived from all donors.
C_panel <- rep(NA_real_, n); panels <- list()
p_all <- panel_from(donors)
for (s in unique(st)) {
  p <- if (s %in% donors) panel_from(setdiff(donors, s)) else p_all
  if (is.null(p)) next
  i <- which(st == s); if (length(i) < 2) next
  Xi <- sweep(A[i, p$genera, drop = FALSE], 2, colMeans(A[i, p$genera, drop = FALSE]), "-")
  C_panel[i] <- as.numeric(Xi %*% p$w) / length(p$genera)
  if (s %in% donors) panels[[s]] <- p$genera
}
log_msg("  C  fixed reference panel  : %d arms across %d studies (%d-genus panel)",
        sum(!is.na(C_panel)), length(unique(st[!is.na(C_panel)])), PANEL_N)
if (!is.null(p_all))
  log_msg("       panel from all donors: %s",
          paste(gsub("^g__", "", p_all$genera), collapse = ", "))
if (length(panels)) {
  tb <- sort(table(unlist(panels)), decreasing = TRUE)
  log_msg("     genera appearing across folds:")
  for (g in names(tb)[seq_len(min(10, length(tb)))])
    log_msg("       %-40s %2d/%2d", gsub("^g__", "", g), tb[[g]], length(panels))
}

# =============================================================================
# 3 - EVALUATION
# =============================================================================
log_msg("[3/5] Evaluation")

within_center <- function(x, y) {
  for (s in unique(st)) {
    i <- st == s & !is.na(x) & !is.na(y)
    if (sum(i) >= 2) { x[i] <- x[i] - mean(x[i]); y[i] <- y[i] - mean(y[i]) }
    else { x[st == s] <- NA; y[st == s] <- NA }
  }
  list(x = x, y = y)
}
perm_p <- function(x, y, g, n = N_PERM) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; g <- g[ok]
  if (length(x) < 10) return(c(NA_real_, NA_real_, length(x)))
  obs <- suppressWarnings(cor(x, y, method = "spearman"))
  idx <- split(seq_along(x), g); cnt <- 0L
  for (i in seq_len(n)) { yp <- y
    for (j in idx) if (length(j) > 1) yp[j] <- y[sample(j)]
    if (suppressWarnings(cor(x, yp, method = "spearman")) <= obs) cnt <- cnt + 1L }
  c(obs, (cnt + 1) / (n + 1), length(x))
}
auc_of <- function(x, pos, neg) {
  s1 <- x[pos & !is.na(x)]; s0 <- x[neg & !is.na(x)]
  if (length(s1) < 5 || length(s0) < 5) return(NA_real_)
  r <- rank(c(s1, s0))
  (sum(r[seq_along(s1)]) - length(s1) * (length(s1) + 1) / 2) / (length(s1) * length(s0))
}

MEAS <- list(A = A_dev, B = B_axis, C = C_panel)
LAB  <- c(A = "deviation from control", B = "signed challenge axis", C = "fixed reference panel")
OOS  <- c(A = FALSE, B = TRUE, C = TRUE)

res <- list()
for (k in names(MEAS)) {
  x <- MEAS[[k]]
  cc <- within_center(x, GHI)
  r  <- perm_p(cc$x, cc$y, st)
  a_xc <- if (OOS[[k]]) auc_of(x, isX, isC) else NA_real_    # meaningless for A
  a_xa <- auc_of(x, isX, isA)   # now computable for all three
  res[[k]] <- list(rho = r[1], p = r[2], n = r[3], auc_xc = a_xc, auc_xa = a_xa)
}

log_msg("")
log_msg("  %-24s %8s %9s %6s %10s %10s", "index", "rho GHI", "perm p", "n",
        "AUC X/ctrl", "AUC X/add")
for (k in names(MEAS)) {
  z <- res[[k]]
  log_msg("  %-24s %+8.3f %9s %6d %10s %10s", LAB[[k]],
          ifelse(is.na(z$rho), NA, z$rho),
          ifelse(is.na(z$p), "n/a", sprintf("%.4f", z$p)), z$n,
          ifelse(is.na(z$auc_xc), "n/a *", sprintf("%.3f", z$auc_xc)),
          ifelse(is.na(z$auc_xa), "n/a", sprintf("%.3f", z$auc_xa)))
}
log_msg("")
log_msg("  rho GHI    within-study Spearman against the GHI; NEGATIVE is working")
log_msg("  AUC X/ctrl challenged vs control; not computed for A, where controls")
log_msg("             sit at zero by construction and would separate trivially")
log_msg("  AUC X/add  challenged vs additive; the fair test for all three, since")
log_msg("             both groups differ from control")

# =============================================================================
# 4 - WHICH, IF ANY
# =============================================================================
log_msg("[4/5] Reading")
best <- NULL; bestv <- 0
for (k in names(res)) {
  z <- res[[k]]
  if (!is.na(z$rho) && z$rho < 0 && !is.na(z$p) && z$p < 0.05 && abs(z$rho) > bestv) {
    best <- k; bestv <- abs(z$rho)
  }
}
verdict <- if (is.null(best)) {
  paste("None of the three correlates with the GHI in the working direction at",
        "p < 0.05. Note that A is inherently within-study and cannot be applied",
        "to a new study without that study's own controls, so even a positive",
        "result for A would not give a transferable index.")
} else if (OOS[[best]]) {
  sprintf(paste("%s performs best and is validated out-of-study (rho = %+.3f,",
                "p = %.4f). It can be applied to a study it was not built on,",
                "which makes it a usable index."), LAB[[best]], res[[best]]$rho, res[[best]]$p)
} else {
  sprintf(paste("%s performs best (rho = %+.3f, p = %.4f) but is inherently",
                "within-study: scoring a new sample requires control arms from",
                "that same study. Usable for comparing arms inside one trial,",
                "not as a standalone index."), LAB[[best]], res[[best]]$rho, res[[best]]$p)
}
log_msg("  %s", verdict)

# =============================================================================
# 5 - WRITE
# =============================================================================
log_msg("[5/5] Writing")
out <- c("THREE 16S-ONLY DYSBIOSIS INDICES", strrep("=", 74), "",
  "A  deviation from the study's own control centroid (within-study only)",
  "B  signed control-to-challenge axis, built from other studies",
  sprintf("C  fixed %d-genus reference panel, built from other studies", PANEL_N), "",
  sprintf("Arms: %d across %d studies; %d studies contain both arm types",
          nrow(A), length(unique(st)), length(donors)), "",
  sprintf("  %-24s %8s %9s %6s %10s %10s", "index", "rho GHI", "perm p", "n",
          "AUC X/ctrl", "AUC X/add"))
for (k in names(MEAS)) {
  z <- res[[k]]
  out <- c(out, sprintf("  %-24s %+8.3f %9s %6d %10s %10s", LAB[[k]],
           ifelse(is.na(z$rho), NA, z$rho),
           ifelse(is.na(z$p), "n/a", sprintf("%.4f", z$p)), z$n,
           ifelse(is.na(z$auc_xc), "n/a", sprintf("%.3f", z$auc_xc)),
           ifelse(is.na(z$auc_xa), "n/a", sprintf("%.3f", z$auc_xa))))
}
out <- c(out, "", verdict, "",
  "A cannot be validated out-of-study because it requires control arms from",
  "the study being scored. B and C are built from the training studies only",
  "and applied to the held-out study, so their figures are out-of-sample.",
  "AUC challenged-versus-control is not computed for A: its control arms sit",
  "at distance zero by construction and would separate perfectly for a trivial",
  "reason.")
writeLines(out, file.path(R_OUT_DIR, "gdi", "gdi_three_ways.txt"))
fwrite(data.table(arm = rownames(A), study = st, class = cls,
                  A_deviation = A_dev, B_axis = B_axis, C_panel = C_panel,
                  GHI_Combined = GHI),
       file.path(R_OUT_DIR, "gdi", "gdi_three_ways_scores.tsv"), sep = "\t")
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "gdi_three_ways.txt"))
log_msg("=== DONE ===")
