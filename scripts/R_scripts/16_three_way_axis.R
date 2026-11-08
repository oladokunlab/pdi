#!/usr/bin/env Rscript
# =============================================================================
# 16_three_way_axis.R
#
# A GHI-free GDI built from all three groups: control, additive, challenge.
#
# RUN:  Rscript 16_three_way_axis.R
# TIME: under a minute.
#
# -----------------------------------------------------------------------------
# WHY THIS IS DIFFERENT FROM THE TWO-GROUP AXIS
# -----------------------------------------------------------------------------
# The earlier signed axis used only control and challenge arms, needing a
# study to contain both - which only 13 of 42 studies do, since most are
# additive-only trials with no challenge condition.
#
# This version uses all three groups. Within each study, the control arms
# define the anchor point. Two directions are then computed from that anchor:
#   control -> challenge   (the "getting worse" direction)
#   control -> additive    (the "getting better" direction)
# Every arm - including the additive ones that the two-group axis could not
# use - is scored by its position along both directions.
#
# Nothing here is fitted to the GHI. The GHI is used only afterward, to check
# whether the resulting scores line up with it.
#
# -----------------------------------------------------------------------------
# WHY TWO DIRECTIONS RATHER THAN ONE -1/0/+1 AXIS
# -----------------------------------------------------------------------------
# Forcing challenge and additive onto a single straight line assumes the
# microbiome moves the same way, just in opposite directions, for both. There
# is no reason to expect that: a pathogen challenge and a beneficial
# supplement are different biological processes and there is no reason their
# effects on the community should be exact opposites. Reporting two directions
# says what each one does without assuming they are mirror images.
#
# A study needs a control arm plus AT LEAST ONE of challenge or additive to
# contribute a direction; it need not have both.
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
N_PERM <- 5000
set.seed(42)

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

isC <- cls %in% HEALTHY; isX <- cls %in% CHALLENGE; isA <- cls %in% ADDITIVE
log_msg("  %d arms, %d studies (%d control, %d challenge, %d additive)",
        nrow(A), length(unique(st)), sum(isC), sum(isX), sum(isA))

donors <- unique(st)[vapply(unique(st), function(s) {
  i <- st == s
  any(isC[i]) && (any(isX[i]) || any(isA[i]))
}, logical(1))]
log_msg("  %d studies contain a control arm plus challenge and/or additive", length(donors))

# =============================================================================
# 2 - TWO DIRECTIONS PER STUDY, POOLED ACROSS STUDIES
# =============================================================================
log_msg("[2/5] Computing control->challenge and control->additive directions")

dir_from <- function(s, target_set) {
  i <- which(st == s); ic <- i[isC[i]]
  it <- i[cls[i] %in% target_set]
  if (!length(ic) || !length(it)) return(NULL)
  d <- colMeans(A[it, , drop = FALSE]) - colMeans(A[ic, , drop = FALSE])
  L <- sqrt(sum(d^2)); if (!is.finite(L) || L == 0) return(NULL); d / L
}
pooled_dir <- function(studies, target_set) {
  ds <- Filter(Negate(is.null), lapply(studies, dir_from, target_set = target_set))
  if (!length(ds)) return(NULL)
  v <- Reduce(`+`, ds) / length(ds); L <- sqrt(sum(v^2))
  if (!is.finite(L) || L == 0) return(NULL); v / L
}

n_x_donors <- sum(vapply(donors, function(s) !is.null(dir_from(s, CHALLENGE)), logical(1)))
n_a_donors <- sum(vapply(donors, function(s) !is.null(dir_from(s, ADDITIVE)), logical(1)))
log_msg("  studies contributing a control->challenge direction: %d", n_x_donors)
log_msg("  studies contributing a control->additive direction: %d", n_a_donors)

score_challenge <- rep(NA_real_, nrow(A))
score_additive  <- rep(NA_real_, nrow(A))
for (s in unique(st)) {
  others <- setdiff(donors, s)
  i <- which(st == s); if (length(i) < 2) next
  Xi <- sweep(A[i, , drop = FALSE], 2, colMeans(A[i, , drop = FALSE]), "-")

  vx <- pooled_dir(others, CHALLENGE)
  if (!is.null(vx)) score_challenge[i] <- as.numeric(Xi %*% vx)

  va <- pooled_dir(others, ADDITIVE)
  if (!is.null(va)) score_additive[i] <- as.numeric(Xi %*% va)
}
log_msg("  challenge-axis scores: %d arms", sum(!is.na(score_challenge)))
log_msg("  additive-axis scores: %d arms", sum(!is.na(score_additive)))

# =============================================================================
# 3 - ATTACH GHI, EVALUATION ONLY
# =============================================================================
log_msg("[3/5] Attaching GHI (evaluation only)")
gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
gk <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
GHI <- unname(tapply(gh$GHI_Combined_raw, gk, function(z) mean(z, na.rm = TRUE))[rownames(A)])

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
  if (length(x) < 10) return(c(NA_real_, NA_real_, length(x)))
  obs <- suppressWarnings(cor(x, y, method = "spearman"))
  idx <- split(seq_along(x), g); cnt <- 0L
  for (i in seq_len(n)) { yp <- y
    for (j in idx) if (length(j) > 1) yp[j] <- y[sample(j)]
    if (suppressWarnings(cor(x, yp, method = "spearman")) <= obs) cnt <- cnt + 1L }
  c(obs, (cnt + 1) / (n + 1), length(x))
}

log_msg("[4/5] Correlation with the GHI, and separation by group")
res <- list()
for (nm in c("challenge", "additive")) {
  x <- if (nm == "challenge") score_challenge else score_additive
  cc <- within_center(x, GHI, st)
  r <- perm_p(cc$x, cc$y, st)
  res[[nm]] <- r
  log_msg("  %-10s axis vs GHI: rho = %+.3f  perm p = %s  n = %d",
          nm, r[1], ifelse(is.na(r[2]), "n/a", sprintf("%.4f", r[2])), r[3])
}

auc_of <- function(x, pos, neg) {
  s1 <- x[pos & !is.na(x)]; s0 <- x[neg & !is.na(x)]
  if (length(s1) < 5 || length(s0) < 5) return(NA_real_)
  r <- rank(c(s1, s0))
  (sum(r[seq_along(s1)]) - length(s1)*(length(s1)+1)/2) / (length(s1)*length(s0))
}
auc_xc <- auc_of(score_challenge, isX, isC)
auc_ac <- auc_of(score_additive, isA, isC)
log_msg("")
log_msg("  AUC, challenge-axis: challenge vs control = %.3f", auc_xc)
log_msg("  AUC, additive-axis: additive vs control = %.3f", auc_ac)

log_msg("[5/5] Writing")
out <- c("THREE-WAY GHI-FREE AXIS: CONTROL, CHALLENGE, ADDITIVE", strrep("=", 74), "",
  "Two directions computed per study from control arms as the anchor:",
  "  control -> challenge, and control -> additive. Pooled across donor",
  "studies, applied out-of-study (each study scored by a direction built",
  "without it). No GHI used in construction.", "",
  sprintf("Studies with control + (challenge and/or additive): %d", length(donors)),
  sprintf("  contributing control->challenge: %d", n_x_donors),
  sprintf("  contributing control->additive: %d", n_a_donors), "",
  sprintf("  %-24s %10s %10s %6s", "axis", "rho vs GHI", "perm p", "n"),
  sprintf("  %-24s %+10.3f %10s %6d", "control->challenge", res$challenge[1],
          ifelse(is.na(res$challenge[2]), "n/a", sprintf("%.4f", res$challenge[2])), res$challenge[3]),
  sprintf("  %-24s %+10.3f %10s %6d", "control->additive", res$additive[1],
          ifelse(is.na(res$additive[2]), "n/a", sprintf("%.4f", res$additive[2])), res$additive[3]),
  "",
  sprintf("AUC, challenge-axis (challenge vs control): %.3f", auc_xc),
  sprintf("AUC, additive-axis (additive vs control): %.3f", auc_ac))
writeLines(out, file.path(R_OUT_DIR, "gdi", "three_way_axis.txt"))
fwrite(data.table(arm = rownames(A), study = st, class = cls,
                  score_challenge = score_challenge, score_additive = score_additive,
                  GHI_Combined = GHI),
       file.path(R_OUT_DIR, "gdi", "three_way_axis_scores.tsv"), sep = "\t")
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "three_way_axis.txt"))
log_msg("=== DONE ===")
