# 50_setup.R ------------------------------------------------------------------
# Shared paths, helpers and data loaders for the PDI verification/robustness
# scripts (50-66). Source this at the top of every other script.
#
#   PDI_ROOT  environment variable, or edit the default below.
#   On Grace:  export PDI_ROOT=/scratch/user/$USER/pdi
# -----------------------------------------------------------------------------

PDI_ROOT <- Sys.getenv("PDI_ROOT", unset = "")
if (!nzchar(PDI_ROOT)) stop(paste(
  "PDI_ROOT is not set. The whole 50-series resolves its inputs and outputs",
  "from it, and the previous default of '.' silently built a parallel output",
  "tree under the working directory instead of failing. Run:",
  "  export PDI_ROOT=/path/to/pdi"), call. = FALSE)
if (!dir.exists(file.path(PDI_ROOT, "r_analysis")))
  stop(sprintf("PDI_ROOT=%s has no r_analysis/ directory -- wrong root?", PDI_ROOT),
       call. = FALSE)
DIR_GENUS <- file.path(PDI_ROOT, "r_analysis", "genus")
DIR_GHI   <- file.path(PDI_ROOT, "r_analysis", "ghi")
DIR_GDI   <- file.path(PDI_ROOT, "r_analysis", "gdi")
DIR_META  <- file.path(PDI_ROOT, "scripts")
DIR_OUT   <- file.path(PDI_ROOT, "r_analysis", "verification")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

SEED <- 42
N_PERM <- 9999

# Tee console output to a file, so every analysis leaves a permanent record on
# disk. Three scripts in the original pipeline (31, 40, 42) lacked this and
# produced no output file. Set PDI_NO_SINK=1 to disable when sourcing
# interactively inside an IDE or notebook kernel.
open_log <- function(name) {
  if (identical(Sys.getenv("PDI_NO_SINK"), "1")) return(invisible(NULL))
  con <- file(file.path(DIR_OUT, paste0(name, ".txt")), open = "wt")
  sink(con, split = TRUE)
  invisible(con)
}
close_log <- function(con) {
  if (is.null(con)) return(invisible(NULL))
  while (sink.number() > 0) sink()
  close(con)
  invisible(NULL)
}

hdr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

# --- loaders -----------------------------------------------------------------

load_counts <- function() {
  as.matrix(read.delim(file.path(DIR_GENUS, "genus_counts.tsv"), row.names = 1, check.names = FALSE))
}
load_prop_adj <- function() {
  as.matrix(read.delim(file.path(DIR_GENUS, "genus_proportions_batch_adjusted.tsv"),
                       row.names = 1, check.names = FALSE))
}
load_clr_published <- function() {
  as.matrix(read.delim(file.path(DIR_GENUS, "genus_clr_batch_adjusted.tsv"),
                       row.names = 1, check.names = FALSE))
}
load_seq_meta <- function() {
  m <- read.delim(file.path(DIR_META, "GDI_metadata_qiime2.tsv"), check.names = FALSE,
                  colClasses = "character", quote = "", comment.char = "",
                  na.strings = c("", "NA"))
  m <- m[!startsWith(as.character(m[[1]]), "#"), , drop = FALSE]
  names(m)[1] <- "sample_id"
  m
}
load_ghi <- function() {
  read.delim(file.path(DIR_GHI, "ghi_scores_pca.tsv"), check.names = FALSE)
}

GHI_DOMAINS <- c("Morphology", "Barrier", "Inflammation", "Anti_Inflammatory",
                 "SCFA", "Performance", "Immunoglobulin_Systemic")
GHI_DOM_COLS <- paste0("GHI_dom_", GHI_DOMAINS)
CHALLENGE_CLASSES <- c("Pathogen_Challenge", "Dietary_Feed_Challenge", "Antibiotic_Challenge")

# ---------------------------------------------------------------------------
# Eight-class intervention scheme (applied 2026-10). Supersedes Challenge_Type
# for the additive and challenge contrasts only; Challenge_Type is retained in
# the workbook and in ghi_scores_pca.tsv, and every other script still reads it.
# Classes excluded from the primary challenge contrast on purpose:
#   "Challenged and treated" a third condition, not a version of either parent
# ---------------------------------------------------------------------------
CONTROL_CLASS_V2     <- "Unchallenged control"
ADDITIVE_CLASS_V2    <- "Non-disruptive additive"
CHALLENGE_CLASSES_V2 <- c("Dietary challenge", "Infectious challenge",
                          "Toxin challenge", "Antibiotic challenge",
                          "Environmental challenge")


# arm-level GHI (mean over samples within Study_ID x Treatment_Group)
ghi_by_arm <- function(ghi = load_ghi(), target = "GHI_Combined_raw") {
  a <- aggregate(list(GHI = ghi[[target]]),
                 by = list(Study_ID = ghi$Study_ID, Treatment_Group = ghi$Treatment_Group),
                 FUN = function(z) mean(z, na.rm = TRUE))
  cls <- aggregate(list(Class = ghi$Challenge_Type),
                   by = list(Study_ID = ghi$Study_ID, Treatment_Group = ghi$Treatment_Group),
                   FUN = function(z) z[1])
  age <- aggregate(list(Age_d = suppressWarnings(as.numeric(ghi$Age_d))),
                   by = list(Study_ID = ghi$Study_ID, Treatment_Group = ghi$Treatment_Group),
                   FUN = function(z) mean(z, na.rm = TRUE))
  a <- merge(merge(a, cls), age)
  a[is.finite(a$GHI), , drop = FALSE]
}

# feature matrix (features x samples) averaged to arms, returned as arms x features
matrix_by_arm <- function(M, meta) {
  keep <- intersect(colnames(M), meta$sample_id)
  M <- M[, keep, drop = FALSE]
  mm <- meta[match(keep, meta$sample_id), ]
  arm <- paste(mm$Study_ID, mm$Treatment_Group, sep = "||")
  t(apply(M, 1, function(v) tapply(v, arm, mean)))
}

# within-study centring
center_within <- function(x, grp) x - ave(x, grp, FUN = function(z) mean(z, na.rm = TRUE))

# within-study restricted permutation test for a Spearman correlation
perm_spearman_within <- function(x, y, grp, B = N_PERM, seed = SEED) {
  set.seed(seed)
  obs <- suppressWarnings(cor(x, y, method = "spearman"))
  idx <- split(seq_along(y), grp)
  null <- numeric(B)
  for (b in seq_len(B)) {
    yy <- y
    for (ii in idx) if (length(ii) > 1) yy[ii] <- sample(y[ii])
    null[b] <- suppressWarnings(cor(x, yy, method = "spearman"))
  }
  list(rho = obs, p = (sum(abs(null) >= abs(obs)) + 1) / (B + 1), n = length(y))
}

# paired within-study contrast: treated arms minus control arms, per study
# n_ctrl / n_trt are the observations contributing to THIS contrast, and
# n_arms = n_ctrl + n_trt is the precision proxy used for funnel asymmetry.
# Do not substitute the study's total arm count -- it gives different results.
study_contrast <- function(df, value, treated_classes, invert = FALSE,
                           control_class = "Control") {
  studies <- unique(df$Study_ID); rows <- list()
  for (s in studies) {
    sub <- df[df$Study_ID == s, ]
    ctl <- sub[[value]][sub$Class %in% control_class]
    trt <- sub[[value]][sub$Class %in% treated_classes]
    ctl <- ctl[is.finite(ctl)]; trt <- trt[is.finite(trt)]
    if (length(ctl) > 0 && length(trt) > 0) {
      e <- mean(trt) - mean(ctl)
      rows[[s]] <- data.frame(study = s, effect = if (invert) -e else e,
                              n_ctrl = length(ctl), n_trt = length(trt),
                              n_arms = length(ctl) + length(trt),
                              sd_trt = if (length(trt) > 1) sd(trt) else NA_real_,
                              stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
