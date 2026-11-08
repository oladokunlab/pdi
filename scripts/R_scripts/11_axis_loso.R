#!/usr/bin/env Rscript
# =============================================================================
# 11_axis_loso.R
#
# Does the control-to-challenge axis transfer to a study it was not built on?
#
# RUN:  Rscript 11_axis_loso.R
# TIME: under a minute.
#
# -----------------------------------------------------------------------------
# THE PROBLEM THIS FIXES
# -----------------------------------------------------------------------------
# In 10_unsupervised_gdi.R the axis was defined inside each study, using that
# study's own control and challenged arms, and then used to score the same
# arms. The correlation with the GHI (-0.504 within study) is therefore
# partly circular: the direction and the thing being measured come from the
# same data.
#
# Here the axis is built by POOLING THE OTHER STUDIES and then applied to the
# held-out one. The held-out study contributes nothing to defining what
# "towards dysbiosis" means. A correlation that survives that is a property of
# the biology rather than of the arithmetic.
#
# -----------------------------------------------------------------------------
# HOW THE AXIS IS BUILT AND APPLIED
# -----------------------------------------------------------------------------
# Within each training study, control arms are averaged and challenged arms are
# averaged, and the difference between those two centroids is that study's
# local direction of deterioration. Differencing within study first is what
# removes the study's own baseline, so pooling the directions across studies
# does not simply recover between-study variation. The pooled direction is the
# mean of the per-study directions, normalised to unit length.
#
# The held-out study's arms are then centred on their own abundance means -
# using abundances only, never the GHI - and projected onto that pooled
# direction. Positive means shifted towards the challenged end.
#
# Nothing is fitted, no taxa are selected, and no threshold is applied. The
# only quantity carried across the fold boundary is a direction in genus space.
#
# -----------------------------------------------------------------------------
# WHAT COUNTS AS A RESULT
# -----------------------------------------------------------------------------
# A working index correlates NEGATIVELY with the GHI: further towards the
# challenged end means lower gut health. Significance is assessed by permuting
# the GHI within each held-out study, which preserves the clustering of arms
# and tests the statistic actually reported.
#
# Only studies containing both a control and a challenged arm can contribute a
# direction, so the training set is smaller than the full corpus and the
# per-fold estimates will be noisy. A weaker figure than the within-study
# -0.504 is expected even if the signal is genuine; the question is whether it
# survives at all.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

HEALTHY   <- "Control"
CHALLENGE <- c("Pathogen_Challenge", "Antibiotic_Challenge", "Dietary_Feed_Challenge")
N_PERM    <- 5000
set.seed(42)

# =============================================================================
# 1 - REBUILD THE ARM MATRIX (same construction as 10)
# =============================================================================
log_msg("[1/4] Loading")
clr <- fread(file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv"))
M <- as.matrix(clr[, -1, with = FALSE]); rownames(M) <- clr[[names(clr)[1]]]

md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2s <- setNames(trimws(as.character(md$Study_ID)), rid)
r2t <- setNames(trimws(as.character(md$Treatment_Group)), rid)

M <- M[, colnames(M) %in% names(r2s), drop = FALSE]
runs <- colnames(M)
ak <- paste(r2s[runs], r2t[runs], sep = "||")
A <- t(vapply(unique(ak), function(a) rowMeans(M[, ak == a, drop = FALSE]), numeric(nrow(M))))
st <- sub("\\|\\|.*$", "", rownames(A))

bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
bk <- paste(trimws(bio$Study_ID), trimws(bio$Treatment_Group), sep = "||")
cls <- unname(tapply(trimws(as.character(bio$Challenge_Type)), bk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1] })[rownames(A)])

gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
gk <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
pull <- function(col) if (col %in% names(gh))
  unname(tapply(gh[[col]], gk, function(z) mean(z, na.rm = TRUE))[rownames(A)]) else rep(NA_real_, nrow(A))
GHI <- list(Combined = pull("GHI_Combined_raw"),
            Morphology = pull("GHI_Morphology_raw"),
            Biomarker  = pull("GHI_Biomarker_raw"))

log_msg("  %d arms, %d studies", nrow(A), length(unique(st)))

# studies able to contribute a direction
donor <- vapply(unique(st), function(s)
  any(cls[st == s] %in% HEALTHY) && any(cls[st == s] %in% CHALLENGE), logical(1))
donors <- names(donor)[donor]
log_msg("  %d studies contain BOTH a control and a challenged arm", length(donors))
if (length(donors) < 4) stop("Too few studies with both arm types to hold one out.")

# =============================================================================
# 2 - LEAVE ONE STUDY OUT
# =============================================================================
log_msg("[2/4] Building the axis from other studies, scoring the held-out one")

axis_from <- function(studies) {
  dirs <- list()
  for (s in studies) {
    i  <- which(st == s)
    ic <- i[cls[i] %in% HEALTHY]; ix <- i[cls[i] %in% CHALLENGE]
    if (!length(ic) || !length(ix)) next
    d <- colMeans(A[ix, , drop = FALSE]) - colMeans(A[ic, , drop = FALSE])
    L <- sqrt(sum(d^2))
    if (is.finite(L) && L > 0) dirs[[s]] <- d / L      # unit direction per study
  }
  if (!length(dirs)) return(NULL)
  v <- Reduce(`+`, dirs) / length(dirs)
  L <- sqrt(sum(v^2))
  if (!is.finite(L) || L == 0) return(NULL)
  v / L
}

score <- rep(NA_real_, nrow(A))
n_used <- 0L
for (s_out in donors) {
  ax <- axis_from(setdiff(donors, s_out))
  if (is.null(ax)) next
  i <- which(st == s_out)
  if (length(i) < 2) next
  Xi <- sweep(A[i, , drop = FALSE], 2, colMeans(A[i, , drop = FALSE]), "-")
  score[i] <- as.numeric(Xi %*% ax)
  n_used <- n_used + 1L
}
log_msg("  %d studies scored by an axis built without them", n_used)
log_msg("  %d arms received a score", sum(!is.na(score)))

# =============================================================================
# 3 - CORRELATION WITH THE GHI
# =============================================================================
log_msg("[3/4] Held-out correlation (negative is the working direction)")

perm_p <- function(x, y, g, n = N_PERM) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]; g <- g[ok]
  if (length(x) < 10) return(c(NA_real_, NA_real_, length(x)))
  obs <- suppressWarnings(cor(x, y, method = "spearman"))
  idx <- split(seq_along(x), g); cnt <- 0L
  for (i in seq_len(n)) {
    yp <- y
    for (j in idx) if (length(j) > 1) yp[j] <- y[sample(j)]
    if (suppressWarnings(cor(x, yp, method = "spearman")) <= obs) cnt <- cnt + 1L
  }
  c(obs, (cnt + 1) / (n + 1), length(x))
}

res <- list()
log_msg("")
log_msg("  %-12s %9s %10s %6s", "GHI", "rho", "perm p", "n")
for (nm in names(GHI)) {
  # centre both within the held-out study, so the comparison is a within-study
  # contrast rather than a between-study one
  x <- score; y <- GHI[[nm]]
  for (s in unique(st)) {
    i <- st == s & !is.na(x) & !is.na(y)
    if (sum(i) >= 2) { x[i] <- x[i] - mean(x[i]); y[i] <- y[i] - mean(y[i]) }
    else { x[st == s] <- NA; y[st == s] <- NA }
  }
  r <- perm_p(x, y, st)
  res[[nm]] <- r
  log_msg("  %-12s %+9.3f %10s %6d", nm, r[1],
          ifelse(is.na(r[2]), "n/a", sprintf("%.4f", r[2])), r[3])
}

# per-study, for inspection
log_msg("")
log_msg("  per-study (inspection only, small n per study):")
for (s in donors) {
  i <- st == s & !is.na(score) & !is.na(GHI$Combined)
  if (sum(i) < 3) next
  r <- suppressWarnings(cor(score[i], GHI$Combined[i], method = "spearman"))
  log_msg("    %-24s n=%2d rho=%+.3f", s, sum(i), r)
}

# =============================================================================
# 4 - READING
# =============================================================================
log_msg("[4/4] Reading")
c1 <- res[["Combined"]]
verdict <- if (is.na(c1[1])) {
  "Not evaluable."
} else if (c1[1] <= -0.25 && !is.na(c1[2]) && c1[2] < 0.05) {
  sprintf(paste("The axis transfers to studies it was not built on (rho = %+.3f,",
                "p = %.4f). This is an unsupervised index - no taxa selected, no",
                "threshold applied, no model fitted - so the correlation with the",
                "GHI is external validation. Worth reporting as a community-level",
                "dysbiosis measure."), c1[1], c1[2])
} else if (c1[1] < 0 && !is.na(c1[2]) && c1[2] < 0.05) {
  sprintf(paste("A weak but significant transfer (rho = %+.3f, p = %.4f). Report",
                "with its limitations rather than as a deployable index."), c1[1], c1[2])
} else {
  sprintf(paste("The axis does not transfer (rho = %+.3f, p = %s). The within-study",
                "correlation of -0.504 reported earlier was therefore circular: the",
                "direction and the arms it scored came from the same study. The",
                "negative result stands."), c1[1], ifelse(is.na(c1[2]), "n/a", sprintf("%.4f", c1[2])))
}
log_msg("  %s", verdict)

out <- c("CONTROL-TO-CHALLENGE AXIS, LEAVE-ONE-STUDY-OUT", strrep("=", 74), "",
  "The axis is built by pooling within-study control-to-challenge directions",
  "from the OTHER studies and applied to the held-out one. Nothing is fitted,",
  "no taxa are selected, no threshold is applied.", "",
  sprintf("Studies contributing a direction : %d", length(donors)),
  sprintf("Studies scored out-of-sample     : %d", n_used),
  sprintf("Arms scored                      : %d", sum(!is.na(score))), "",
  sprintf("  %-12s %9s %10s %6s", "GHI", "rho", "perm p", "n"))
for (nm in names(res))
  out <- c(out, sprintf("  %-12s %+9.3f %10s %6d", nm, res[[nm]][1],
                        ifelse(is.na(res[[nm]][2]), "n/a", sprintf("%.4f", res[[nm]][2])),
                        res[[nm]][3]))
out <- c(out, "", verdict, "",
  "Both the score and the GHI are centred within the held-out study before",
  "correlating, so this is a within-study contrast. Permutation shuffles the",
  "GHI within study, preserving the clustering of arms.")
writeLines(out, file.path(R_OUT_DIR, "gdi", "axis_loso.txt"))
fwrite(data.table(arm = rownames(A), study = st, class = cls,
                  axis_score = score, GHI_Combined = GHI$Combined),
       file.path(R_OUT_DIR, "gdi", "axis_loso_scores.tsv"), sep = "\t")
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "axis_loso.txt"))
log_msg("=== DONE ===")
