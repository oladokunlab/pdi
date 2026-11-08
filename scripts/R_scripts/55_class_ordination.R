#!/usr/bin/env Rscript
# =============================================================================
# 55_class_ordination.R
#
# Replaces the QIIME ordinations in merged/ for the class contrast.
#
# Why they are replaced rather than relabelled. The deposited artefacts
# (braycurtis_Intervention_Class_permanova.qzv, braycurtis_PCoA_Group_permanova.qzv)
# are wrong on two counts:
#
#   1. LABELS. Intervention_Class has six levels that do not correspond to the
#      eight-class scheme and conflate insults with additives: its
#      Dietary_Intervention level contains 466 additive runs but also 30 toxin
#      challenge and 22 dietary challenge runs, and its Infectious_Challenge
#      level contains 94 infectious, 59 toxin, 50 challenged-and-treated and
#      15 additive. PCoA_Group additionally buckets seven levels as
#      "Other <something>", which the classification specification forbids.
#
#   2. BASIS. They were computed on merged-table.qza, which holds all 3,037
#      runs from 50 studies and therefore mixes V3-V4 with V4, V3-V5 and
#      full-length amplicons in one distance matrix. Different regions recover
#      different feature sets, so separation in that ordination is partly an
#      artefact of which region a study sequenced.
#
# This script uses the R-side table instead: V3-V4 only, MMUPHin
# batch-corrected, CLR-transformed. Distance is Euclidean on CLR, i.e.
# Aitchison distance, which is the appropriate metric for compositional data
# and needs no rarefaction.
#
# Permutations are BLOCKED WITHIN STUDY. Runs from one study share husbandry,
# genetics, diet base, sequencing run and protocol, so permuting them freely
# would test a null no one is proposing and inflate significance.
#
# betadisper is reported alongside, because PERMANOVA cannot distinguish a
# difference in group centroid from a difference in group dispersion, and
# with eight classes of very unequal size that distinction matters.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

suppressPackageStartupMessages({ library(data.table); library(vegan) })

NPERM     <- 9999
LABELS    <- file.path(PDI_ROOT, "intervention_classification_FINAL.csv")
CLR_TSV   <- file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv")
META_TSV  <- METADATA_TSV

# Studies that cannot be labelled at class resolution, with the reason each
# is dropped. Neither is a quality exclusion.
EXCLUDE <- c(
  Jansseune_2024 = "sequencing is diet-level only: Standard diet and Challenge diet each pool three treatment arms, so assigning a class would place additive birds in the control stratum",
  Massacci_2019  = "class is age-dependent (d21 faecal samples are pre-challenge, d40 post-challenge) and no sequencing-side table carries Age_d"
)

log_msg("[1/5] Loading CLR table, run metadata and final class labels")
clr <- fread(CLR_TSV)
gcol <- names(clr)[1]
X <- t(as.matrix(clr[, -1, with = FALSE]))          # samples x genera
colnames(X) <- clr[[gcol]]
log_msg("  CLR matrix: %d samples x %d genera", nrow(X), ncol(X))

md <- fread(META_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
setnames(md, 1, "Run")
md[, `:=`(Run = trimws(Run), Study_ID = trimws(Study_ID),
          Treatment_Group = trimws(Treatment_Group))]

lab <- fread(LABELS)
lab[, `:=`(Study_ID = trimws(Study_ID), Treatment_Group = trimws(Treatment_Group))]

# Class is a property of the TREATMENT ARM, not of the gut segment a sample was
# taken from, so the join key is Study_ID x Treatment_Group. Joining on Segment
# as well fails for most runs, because sequencing and biomarkers are frequently
# taken from different segments of the same bird.
arms <- unique(lab[, .(Study_ID, Treatment_Group, Class, Subclass)])
dup  <- arms[, .N, by = .(Study_ID, Treatment_Group)][N > 1]
if (nrow(dup)) {
  log_msg("  NOTE: %d arm(s) carry more than one class across segment/age:", nrow(dup))
  for (i in seq_len(nrow(dup))) log_msg("    %s / %s", dup$Study_ID[i], dup$Treatment_Group[i])
  log_msg("  These must be in EXCLUDE or the join is ambiguous.")
}

log_msg("[2/5] Joining classes onto runs")
d <- merge(md[Run %in% rownames(X), .(Run, Study_ID, Treatment_Group, Segment)],
           arms, by = c("Study_ID", "Treatment_Group"), all.x = TRUE)
log_msg("  matched %d of %d runs (%.0f%%)", sum(!is.na(d$Class)), nrow(d),
        100 * mean(!is.na(d$Class)))

for (s in names(EXCLUDE)) {
  n <- sum(d$Study_ID == s)
  if (n > 0) log_msg("  excluding %s (%d runs): %s", s, n, EXCLUDE[[s]])
}
d <- d[!is.na(Class) & !(Study_ID %in% names(EXCLUDE))]
d <- unique(d, by = "Run")
X <- X[d$Run, , drop = FALSE]
X <- X[, apply(X, 2, function(v) sd(v) > 0), drop = FALSE]
log_msg("  ordination set: %d runs, %d studies, %d classes",
        nrow(d), uniqueN(d$Study_ID), uniqueN(d$Class))

log_msg("[3/5] Aitchison ordination (PCA on CLR)")
pc <- prcomp(X, center = TRUE, scale. = FALSE)
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
d[, `:=`(PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC3 = pc$x[, 3])]
log_msg("  PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%", ve[1], ve[2], ve[3])

log_msg("[4/5] PERMANOVA, permutations blocked within study")
D   <- dist(X)                                   # Euclidean on CLR = Aitchison
ctl <- how(blocks = factor(d$Study_ID), nperm = NPERM)
set.seed(1)
aC <- adonis2(D ~ Class, data = d, permutations = ctl, by = "terms")
aS <- adonis2(D ~ Segment, data = d, permutations = ctl, by = "terms")
bd <- betadisper(D, factor(d$Class))
set.seed(1)
pbd <- permutest(bd, permutations = NPERM)

log_msg("[5/5] Writing outputs")
out <- c(
  "=== Class ordination and PERMANOVA (Aitchison, V3-V4, batch-adjusted) ===", "",
  sprintf("Runs: %d   Studies: %d   Classes: %d   Genera: %d",
          nrow(d), uniqueN(d$Study_ID), uniqueN(d$Class), ncol(X)),
  sprintf("Permutations: %d, blocked within Study_ID", NPERM), "",
  "Replaces merged/braycurtis_Intervention_Class_permanova.qzv and",
  "merged/braycurtis_PCoA_Group_permanova.qzv. See the header of this script",
  "for why those two are retired rather than relabelled.", "",
  "Excluded, neither for quality reasons:",
  unname(sprintf("  %s: %s", names(EXCLUDE), EXCLUDE)), "",
  "--- variance explained by the ordination axes ---",
  sprintf("  PC1 %.2f%%   PC2 %.2f%%   PC3 %.2f%%", ve[1], ve[2], ve[3]), "",
  "--- PERMANOVA ---",
  capture.output(print(aC)), "",
  capture.output(print(aS)), "",
  "--- betadisper (homogeneity of multivariate dispersion) ---",
  "PERMANOVA cannot separate a centroid shift from a dispersion difference.",
  "A significant betadisper means the classes differ in spread as well as, or",
  "instead of, in location, and the PERMANOVA must be read with that in mind.",
  capture.output(print(pbd)), "",
  "--- runs per class ---",
  capture.output(print(d[, .N, by = Class][order(-N)])))

dir.create(file.path(R_OUT_DIR, "verification"), recursive = TRUE, showWarnings = FALSE)
f1 <- file.path(R_OUT_DIR, "verification", "55_class_ordination.txt")
writeLines(out, f1); log_msg("  -> %s", f1)

f2 <- file.path(R_OUT_DIR, "gdi", "pcoa_strata.tsv")
fwrite(d[, .(Run, Study_ID, Treatment_Group, Segment, Class, Subclass, PC1, PC2, PC3)],
       f2, sep = "\t")
log_msg("  -> %s", f2)

f3 <- file.path(R_OUT_DIR, "gdi", "pcoa_varexp.tsv")
fwrite(data.table(axis = paste0("PC", 1:3), varexp_pct = round(ve[1:3], 3)), f3, sep = "\t")
log_msg("  -> %s", f3)
log_msg("=== DONE ===")
