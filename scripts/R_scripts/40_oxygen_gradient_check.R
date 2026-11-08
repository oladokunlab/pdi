#!/usr/bin/env Rscript
# =============================================================================
# 40_oxygen_gradient_check.R
#
# Does this corpus reproduce the anaerobe-toward-cecum gradient reported in
# Fig. 2f-g of Bindari et al. (Microbiome 2021, s40168-021-01040-x)?
#
# RUN:  Rscript 40_oxygen_gradient_check.R
#
# -----------------------------------------------------------------------------
# WHAT THIS IS AND IS NOT
# -----------------------------------------------------------------------------
# The cited paper used BugBase to predict oxygen tolerance and biofilm
# formation, then compared abundance of each category across five sample
# types (duodenum, jejunum, ileum, cecum, feces) with pairwise Wilcoxon
# rank-sum tests. Anaerobes rose sharply toward the cecum; aerobes dominated
# the small intestine.
#
# This project has no biofilm-formation trait - checked directly against the
# BactoTraits column list, it is not there. Oxygen tolerance IS available
# (oxygen_tolerance_aerobe/anaerobe/facultative/microaerophile), the same
# broad trait the cited paper used, though from a different source (BacDive-
# derived community weighted means, not BugBase's phylogeny-based
# prediction).
#
# This is NOT a test of the GDI. It makes no reference to the GHI and
# predicts nothing. It is a validity check: if this corpus's sequencing and
# taxonomy pipeline are sound, the same well-established gut physiology
# (oxygen availability falls sharply from duodenum to cecum, driving the
# aerobe-to-anaerobe shift) should be visible in it independently of
# anything built for the GDI search. A missing gradient would be a warning
# sign about the pipeline; a matching gradient is expected and mildly
# reassuring, not a discovery, since it is the guild-ratio anaerobe result
# already tested in 33_guild_ratios.R (against GHI, not segment) computed a
# different way.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
con <- file(file.path(R_OUT_DIR, "gdi", "40_oxygen_gradient_check.txt"), open = "wt"); sink(con, split = TRUE)

log_msg("[1/5] Loading BactoTraits and matching genera (same procedure as 32_bactotraits.R)")
BT <- Sys.getenv("BACTOTRAITS_CSV",
      unset = file.path(PDI_ROOT, "reference", "BACTOTRAITS_DATASET_2026-01-28_GENUSLVL.csv"))
bt <- as.data.frame(fread(BT, dec = ",", check.names = FALSE), stringsAsFactors = FALSE)
oxy_cols <- c("oxygen_tolerance_facultative","oxygen_tolerance_microaerophile",
              "oxygen_tolerance_anaerobe","oxygen_tolerance_aerobe")
for (cn in oxy_cols) {
  v <- bt[[cn]]
  if (is.character(v)) v <- suppressWarnings(as.numeric(gsub(",", ".", v)))
  bt[[cn]] <- v
}
bt$Genus <- trimws(as.character(bt$Genus))

log_msg("[2/5] Loading run-level genus proportions and segment labels")
prop <- fread(file.path(R_OUT_DIR, "genus", "genus_proportions_batch_adjusted.tsv"))
gname <- prop[[1]]; P <- as.matrix(prop[, -1, with = FALSE]); rownames(P) <- gname
gm <- gsub("^g__", "", gname); gm[gm == "Escherichia-Shigella"] <- "Escherichia"
idx <- match(gm, bt$Genus); matched <- !is.na(idx)
log_msg("  %d of %d genera matched to BactoTraits (%.0f%% by name)",
        sum(matched), length(gm), 100*sum(matched)/length(gm))

T <- as.matrix(bt[idx[matched], oxy_cols, drop = FALSE]); T[!is.finite(T)] <- 0
Pm <- P[matched, , drop = FALSE]
w <- sweep(Pm, 2, pmax(colSums(Pm, na.rm = TRUE), 1e-12), "/")
CWM <- t(T) %*% w   # 4 oxygen categories x runs
CWM <- t(CWM)       # runs x 4

log_msg("[3/5] Attaching gut segment labels")
md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2seg <- setNames(trimws(as.character(md$Segment)), rid)

common <- intersect(rownames(CWM), names(r2seg))
CWM <- CWM[common, , drop = FALSE]
seg <- r2seg[common]
seg[seg == "" | is.na(seg)] <- NA
ok <- !is.na(seg)
CWM <- CWM[ok, , drop = FALSE]; seg <- seg[ok]

log_msg("  segments present in this corpus, with sample counts:")
tb <- sort(table(seg), decreasing = TRUE)
for (s in names(tb)) log_msg("    %-20s %4d runs", s, tb[[s]])

log_msg("")
log_msg("[4/5] Anaerobe and aerobe abundance by segment (median CWM, oriented low to high O2)")
# Order segments anterior to posterior where recognizable; unrecognized
# segment labels are kept but placed at the end rather than dropped, so
# nothing is silently excluded.
canonical_order <- c("Duodenum","Jejunum","Ileum","Cecum","Colon","Feces")
present <- names(tb)
ord <- c(intersect(canonical_order, present), setdiff(present, canonical_order))

# Diagnostics first - the previous run produced a header with zero data rows
# and no error, which is not explainable from reading the code alone. Rather
# than guess a fix blind, print every intermediate quantity the table build
# depends on.
log_msg("  DIAGNOSTIC: class(CWM)=%s dim(CWM)=%s", class(CWM)[1], paste(dim(CWM), collapse=" x "))
log_msg("  DIAGNOSTIC: colnames(CWM)=%s", paste(colnames(CWM), collapse=", "))
log_msg("  DIAGNOSTIC: length(seg)=%d nrow(CWM)=%d (must match)", length(seg), nrow(CWM))
log_msg("  DIAGNOSTIC: ord=%s (length %d)", paste(ord, collapse=", "), length(ord))
log_msg("  DIAGNOSTIC: class(tb)=%s", class(tb)[1])

n_vals <- as.integer(tb[ord])
anaerobe_vals <- sapply(ord, function(s) {
  v <- CWM[seg == s, "oxygen_tolerance_anaerobe"]
  log_msg("    DIAGNOSTIC segment=%s: n matched rows=%d, any finite=%s",
          s, length(v), any(is.finite(v)))
  median(v, na.rm = TRUE)
})
aerobe_vals <- sapply(ord, function(s) median(CWM[seg == s, "oxygen_tolerance_aerobe"], na.rm = TRUE))
fac_vals <- sapply(ord, function(s) median(CWM[seg == s, "oxygen_tolerance_facultative"], na.rm = TRUE))

log_msg("  DIAGNOSTIC: length(n_vals)=%d length(anaerobe_vals)=%d", length(n_vals), length(anaerobe_vals))

log_msg("")
log_msg("  %-16s %8s %10s %10s %10s", "segment", "n", "anaerobe", "aerobe", "facultative")
for (i in seq_along(ord)) {
  log_msg("  %-16s %8d %10.3f %10.3f %10.3f",
          ord[i], n_vals[i], anaerobe_vals[i], aerobe_vals[i], fac_vals[i])
}

log_msg("")
log_msg("[5/5] Pairwise Wilcoxon rank-sum tests, same method as the cited paper")
log_msg("  (Bindari et al. 2021: 'significances of the comparisons were determined")
log_msg("  using the pairwise Wilcoxon rank-sum test')")
log_msg("")
run_pairwise <- function(trait_col, trait_name) {
  log_msg("  --- %s ---", trait_name)
  segs_with_data <- present[tb[present] >= 3]
  if (length(segs_with_data) < 2) { log_msg("    fewer than 2 segments with >=3 runs; skipped"); return(invisible(NULL)) }
  x <- CWM[seg %in% segs_with_data, trait_col]
  g <- factor(seg[seg %in% segs_with_data], levels = segs_with_data)
  pw <- pairwise.wilcox.test(x, g, p.adjust.method = "BH")
  print(pw)
}
run_pairwise("oxygen_tolerance_anaerobe", "Anaerobe")
run_pairwise("oxygen_tolerance_aerobe", "Aerobe")

log_msg("")
log_msg("=== DONE ===")
log_msg("")
log_msg("READING: if anaerobe CWM rises and aerobe CWM falls moving from the small")
log_msg("intestine (Duodenum/Jejunum/Ileum) toward the Cecum, this corpus reproduces")
log_msg("the physiological gradient the cited paper reports in its Fig. 2f-g. That is")
log_msg("a validity check on this project's sequencing and taxonomy pipeline, not a")
log_msg("new predictor - it says nothing about which genera track the GHI and is not")
log_msg("a candidate for the GDI panel.")
sink(); close(con)
