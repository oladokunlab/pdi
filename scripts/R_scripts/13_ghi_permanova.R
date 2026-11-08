#!/usr/bin/env Rscript
# =============================================================================
# 13_ghi_permanova.R
#
# Across every study, does the Gut Health Index explain ANY of the variation in
# community composition?
#
# RUN:  Rscript 13_ghi_permanova.R
# TIME: 2-4 minutes.
#
# -----------------------------------------------------------------------------
# WHY THIS IS THE CLEANEST VERSION OF THE QUESTION
# -----------------------------------------------------------------------------
# Every previous test asked whether some CONSTRUCTED index tracks gut health -
# a selected panel, a centroid distance, a projection onto an axis. Each
# depends on how the index was built, so a null result can always be blamed on
# the construction.
#
# This asks the question directly and without any construction. PERMANOVA takes
# the community distance matrix and asks how much of its variation is
# attributable to the GHI, entered as a continuous term. Nothing is selected,
# nothing is fitted to an outcome, no threshold is applied, and no treatment
# labels are used.
#
# Age is included alongside it, because age is known to dominate cecal
# composition and because the marginal test then reports what the GHI explains
# INDEPENDENTLY of age rather than confounded with it. Permutations are
# restricted within study, so the test asks whether the GHI tracks composition
# within studies rather than between them.
#
# -----------------------------------------------------------------------------
# WHAT THE ANSWER MEANS
# -----------------------------------------------------------------------------
# A non-negligible, significant R2 for the GHI would mean the microbiome does
# carry gut health information and that previous failures were failures of
# construction - worth returning to with a different estimator.
#
# A negligible or non-significant R2 closes the question. It would mean no
# index built from these genera can predict this outcome, because the
# relationship is not present in the community data to begin with, and the
# negative result is a property of the data rather than of any method.
# =============================================================================
source("00_r_config.R")
suppressPackageStartupMessages(library(vegan))

log_msg("[1/3] Building the arm-level matrix")
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

gh <- fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"))
gk <- paste(trimws(gh$Study_ID), trimws(gh$Treatment_Group), sep = "||")
pull <- function(col) {
  if (!col %in% names(gh)) return(rep(NA_real_, nrow(A)))
  v <- tapply(gh[[col]], gk, function(z) mean(suppressWarnings(as.numeric(z)), na.rm = TRUE))
  unname(v[rownames(A)])
}
GHI <- pull("GHI_Combined_raw"); AGE <- pull("Age_d")

log_msg("  %d arms across %d studies", nrow(A), length(unique(st)))
log_msg("  GHI available for %d arms, age for %d", sum(!is.na(GHI)), sum(!is.na(AGE)))

ok <- !is.na(GHI) & !is.na(AGE)
tb <- table(st[ok]); ok <- ok & st %in% names(tb)[tb >= 2]
log_msg("  %d arms retained across %d studies (GHI, age, and >=2 arms per study)",
        sum(ok), length(unique(st[ok])))
if (sum(ok) < 30) stop("Too few arms.")

log_msg("[2/3] PERMANOVA on Aitchison distance, permutations restricted within study")
D <- dist(A[ok, , drop = FALSE], method = "euclidean")
d <- data.frame(GHI = GHI[ok], Age_d = AGE[ok], Study_ID = factor(st[ok]))
set.seed(42)
fit <- adonis2(D ~ GHI + Age_d, data = d, by = "margin",
               permutations = 9999, strata = d$Study_ID)
print(fit)

r2  <- setNames(fit$R2, rownames(fit))
pv  <- setNames(fit[["Pr(>F)"]], rownames(fit))
log_msg("")
log_msg("  GHI    R2 = %.4f   p = %s  (9999 blocked permutations)", r2[["GHI"]],
        ifelse(is.na(pv[["GHI"]]), "n/a", sprintf("%.3f", pv[["GHI"]])))
log_msg("  Age_d  R2 = %.4f   p = %s  (9999 blocked permutations)", r2[["Age_d"]],
        ifelse(is.na(pv[["Age_d"]]), "n/a", sprintf("%.3f", pv[["Age_d"]])))

log_msg("[3/3] Reading")
v <- if (r2[["GHI"]] >= 0.03 && !is.na(pv[["GHI"]]) && pv[["GHI"]] < 0.05) {
  sprintf(paste("The GHI explains %.1f%% of community variation independently of",
                "age, significantly. The microbiome does carry gut health",
                "information, and previous failures were failures of construction",
                "rather than absence of signal."), 100*r2[["GHI"]])
} else {
  sprintf(paste("The GHI explains %.1f%% of community variation (p = %s), which is",
                "negligible. No index built from these genera can predict this",
                "outcome, because the relationship is not present in the community",
                "data. The negative result is a property of the data, not of any",
                "particular method."), 100*r2[["GHI"]],
          ifelse(is.na(pv[["GHI"]]), "n/a", sprintf("%.3f", pv[["GHI"]])))
}
log_msg("  %s", v)

writeLines(c("DOES THE GHI EXPLAIN COMMUNITY VARIATION?", strrep("=", 74), "",
  "PERMANOVA on Aitchison distance, GHI entered as a continuous term alongside",
  "age, marginal tests, permutations restricted within study. No taxa selected,",
  "nothing fitted to an outcome, no treatment labels used.", "",
  sprintf("Arms    : %d", sum(ok)),
  sprintf("Studies : %d", length(unique(st[ok]))), "",
  sprintf("GHI    R2 = %.4f   p = %s", r2[["GHI"]],
          ifelse(is.na(pv[["GHI"]]), "n/a", sprintf("%.3f", pv[["GHI"]]))),
  sprintf("Age_d  R2 = %.4f   p = %s", r2[["Age_d"]],
          ifelse(is.na(pv[["Age_d"]]), "n/a", sprintf("%.3f", pv[["Age_d"]]))), "",
  "p-values from 9999 permutations restricted within study; resolution 1e-4.",
  "This is a BLOCKED test and is not comparable with the free-permutation",
  "age result in 54_permanova_five_configs.R.", "", v),
  file.path(R_OUT_DIR, "gdi", "ghi_permanova.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "ghi_permanova.txt"))
log_msg("=== DONE ===")
