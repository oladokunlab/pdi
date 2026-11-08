#!/usr/bin/env Rscript
# =============================================================================
# 17_maaslin_class.R
#
# MaAsLin2, GHI-free: tests genera against the CONTROL/ADDITIVE/CHALLENGE
# label directly, rather than against the GHI.
#
# RUN:  Rscript 17_maaslin_class.R
# TIME: 2-5 minutes.
#
# -----------------------------------------------------------------------------
# WHY THIS COMPLETES THE GHI-FREE SIDE
# -----------------------------------------------------------------------------
# Every GHI-free result so far (the signed axis, the fixed panel, the
# three-way axis) was built from simple mean differences between groups -
# no statistical model, no random effect for study, no correction for
# multiple testing across genera. MaAsLin2 was used only on the GHI-based
# side, as the screening step feeding LASSO.
#
# This runs MaAsLin2 itself as the GHI-free test: does each genus associate
# with the CHALLENGE/ADDITIVE label, treating study as a random effect (as
# in the registered protocol) and controlling the false discovery rate
# across all 84 genera at once. This is a properly specified test where the
# earlier GHI-free measures were descriptive.
#
# Two fixed-effect terms are fitted together (challenge, additive), both
# against control as the reference level, in one model per genus - so each
# genus gets a coefficient and a q-value for "how it differs under
# challenge" and a separate one for "how it differs under additive",
# adjusted for the same study-level random effect throughout.
#
# Genera surviving this screen can then be checked against the genera that
# drive the unsupervised axes, to see whether the two independent methods
# implicate the same taxa.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("Maaslin2", quietly = TRUE)) stop("Package 'Maaslin2' required.")

HEALTHY   <- "Control"
CHALLENGE <- c("Pathogen_Challenge", "Antibiotic_Challenge", "Dietary_Feed_Challenge")
ADDITIVE  <- "Non_Disruptive_Additive"
FDR_THRESHOLD <- 0.25   # same permissive screening threshold as 05_maaslin2.R

log_msg("[1/5] Loading community data and labels")
clr <- fread(file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv"))
gcol <- names(clr)[1]
clr_mat <- as.matrix(clr[, -1, with = FALSE]); rownames(clr_mat) <- clr[[gcol]]

md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2s <- setNames(trimws(as.character(md$Study_ID)), rid)
r2t <- setNames(trimws(as.character(md$Treatment_Group)), rid)

bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
bk <- paste(trimws(bio$Study_ID), trimws(bio$Treatment_Group), sep = "||")
cls_lut <- tapply(trimws(as.character(bio$Challenge_Type)), bk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing = TRUE))[1]
})

runs <- colnames(clr_mat)
runs <- runs[runs %in% names(r2s)]
clr_mat <- clr_mat[, runs, drop = FALSE]
rk <- paste(r2s[runs], r2t[runs], sep = "||")
cls <- unname(cls_lut[rk])

grp <- rep(NA_character_, length(cls))
grp[cls %in% HEALTHY]   <- "Control"
grp[cls %in% CHALLENGE] <- "Challenge"
grp[cls %in% ADDITIVE]  <- "Additive"

ok <- !is.na(grp)
clr_mat <- clr_mat[, ok, drop = FALSE]
meta <- data.frame(
  row.names = colnames(clr_mat),
  Group     = factor(grp[ok], levels = c("Control", "Challenge", "Additive")),
  Study_ID  = factor(r2s[colnames(clr_mat)])
)
log_msg("  %d runs across %d studies (%d control, %d challenge, %d additive)",
        nrow(meta), nlevels(droplevels(meta$Study_ID)),
        sum(meta$Group == "Control"), sum(meta$Group == "Challenge"),
        sum(meta$Group == "Additive"))

# =============================================================================
# 2 - MAASLIN2, GROUP AS FIXED EFFECT, STUDY AS RANDOM EFFECT
#
# Same normalization/transform settings as the registered GHI-based run:
# input is already CLR, so both are NONE.
# =============================================================================
log_msg("[2/5] Running MaAsLin2 (Group as fixed effect, Study_ID as random effect)")
out_dir <- file.path(R_OUT_DIR, "gdi", "maaslin2_class")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit <- Maaslin2::Maaslin2(
  input_data      = as.data.frame(t(clr_mat)),
  input_metadata  = meta,
  output          = out_dir,
  fixed_effects   = c("Group"),
  random_effects  = c("Study_ID"),
  reference       = c("Group,Control"),
  normalization   = "NONE",
  transform       = "NONE",
  min_prevalence  = 0,
  min_abundance   = 0,
  max_significance = FDR_THRESHOLD,
  plot_heatmap    = FALSE,
  plot_scatter    = FALSE,
  standardize     = FALSE
)

log_msg("[3/5] Extracting results")
res_file <- file.path(out_dir, "all_results.tsv")
if (!file.exists(res_file)) stop("MaAsLin2 produced no results file.")
res <- fread(res_file, sep = "\t", header = TRUE)

for (term in c("GroupChallenge", "GroupAdditive")) {
  sub <- res[res$metadata == "Group" & res$value == sub("Group", "", term), ]
  if (!nrow(sub)) next
  setorder(sub, qval)
  n_sig <- sum(sub$qval <= FDR_THRESHOLD, na.rm = TRUE)
  log_msg("")
  log_msg("  %s vs Control: %d of %d genera pass FDR <= %.2f", term, n_sig, nrow(sub), FDR_THRESHOLD)
  top <- head(sub, 10)
  for (i in seq_len(nrow(top)))
    log_msg("    %-32s coef=%+.4f  q=%.4g%s", top$feature[i], top$coef[i], top$qval[i],
            if (top$qval[i] <= FDR_THRESHOLD) "  *" else "")
}

# =============================================================================
# 4 - COMPARE AGAINST THE UNSUPERVISED AXIS GENERA
# =============================================================================
log_msg("[4/5] Comparing against the genera behind the unsupervised axes")
axis_genera <- c("g__Lactobacillus", "g__Streptococcus", "g__Anaerostipes",
                 "g__Colidextribacter", "g__Limosilactobacillus", "g__Tyzzerella")
challenge_res <- res[res$metadata == "Group" & res$value == "Challenge", ]
setorder(challenge_res, qval)
match_tbl <- challenge_res[challenge_res$feature %in% axis_genera, ]
if (nrow(match_tbl)) {
  log_msg("  Of the genera behind the signed control-to-challenge axis,")
  log_msg("  here is what MaAsLin2 independently found for each, vs Challenge:")
  for (i in seq_len(nrow(match_tbl)))
    log_msg("    %-32s coef=%+.4f  q=%.4g%s", match_tbl$feature[i], match_tbl$coef[i],
            match_tbl$qval[i], if (match_tbl$qval[i] <= FDR_THRESHOLD) "  * agrees" else "  not significant here")
} else {
  log_msg("  None of the axis genera appear in the MaAsLin2 results table.")
}

log_msg("[5/5] Done")
log_msg("full results -> %s", res_file)
log_msg("=== DONE ===")
