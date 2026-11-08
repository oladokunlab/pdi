#!/usr/bin/env Rscript
# =============================================================================
# 36_siga_orientation.R
#
# Two questions:
#   1. Does sIgA rise or fall under challenge? The mucosal domain scores
#      AUC 0.312 - well below chance - which says the +1 orientation is wrong.
#   2. What happens if the two immunoglobulin domains are combined?
#
# RUN:  Rscript 36_siga_orientation.R
# TIME: under a minute.
#
# -----------------------------------------------------------------------------
# WHY THE ORIENTATION IS NOW IN QUESTION
# -----------------------------------------------------------------------------
# Before the compartment fix, the mucosal domain mixed secretory IgA with serum
# IgA and scored AUC 0.588 - weak but above chance. With serum IgA moved to the
# systemic domain, the mucosal domain is sIgA alone and scores 0.312. The
# mixing had been pulling it toward the middle and hiding what sIgA does.
#
# 0.312 means challenged birds had HIGHER sIgA than controls, consistently.
# That is biologically coherent: secretory IgA is the frontline mucosal
# antibody response, upregulated when the gut is under attack. It marks an
# active response, not a quiet healthy gut.
#
# This script's own metadata already applies that reasoning to systemic Ig:
# "circulating immunoglobulin rises in RESPONSE to antigenic exposure... the
# +1 is retained ONLY so the domain is oriented consistently for PCA; it is
# not a claim that higher is healthier."
#
# The same caveat plainly applies to sIgA, and was simply not visible while the
# domain was diluted.
#
# -----------------------------------------------------------------------------
# WHAT WOULD MAKE A SIGN FLIP LEGITIMATE
# -----------------------------------------------------------------------------
# Flipping because a domain scores badly is fitting to the outcome. Flipping
# because most studies independently show challenged birds higher is a
# statement about biology that happens to improve the score.
#
# So the test is per-study direction, not the pooled AUC. If the pattern is
# consistent across studies, the sign is wrong. If it comes from two or three
# studies, it is noise and the sign stays.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

HEALTHY   <- "Control"
ADDITIVE  <- "Non_Disruptive_Additive"
CHALLENGE <- c("Pathogen_Challenge","Antibiotic_Challenge","Dietary_Feed_Challenge")

auc_of <- function(pos, neg) {
  pos <- pos[is.finite(pos)]; neg <- neg[is.finite(neg)]
  if (length(pos) < 3 || length(neg) < 3) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos)*(length(pos)+1)/2)/(length(pos)*length(neg))
}

log_msg("[1/4] Loading")
bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
bio$Study_ID <- trimws(bio$Study_ID); bio$cls <- trimws(as.character(bio$Challenge_Type))
gh <- as.data.frame(fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")))
gh$Study_ID <- trimws(gh$Study_ID); gh$cls <- trimws(as.character(gh$Challenge_Type))
dom <- grep("^GHI_dom_", names(gh), value = TRUE)

# =============================================================================
log_msg("")
log_msg("=========================================================")
log_msg("TEST 1 - sIgA direction, study by study")
log_msg("=========================================================")
log_msg("  The pooled AUC of 0.312 could come from a consistent pattern or from")
log_msg("  a few dominant studies. Only the former justifies changing the sign.")
log_msg("")
v <- suppressWarnings(as.numeric(bio$sIgA))
log_msg("  sIgA available for %d rows across %d studies",
        sum(is.finite(v)), length(unique(bio$Study_ID[is.finite(v)])))
log_msg("")
log_msg("  %-24s %12s %12s %11s  %s", "study", "control", "challenged", "difference", "direction")
hi <- 0; lo <- 0
for (s in unique(bio$Study_ID)) {
  k <- bio$Study_ID == s & is.finite(v)
  a <- v[k & bio$cls %in% HEALTHY]; b <- v[k & bio$cls %in% CHALLENGE]
  if (!length(a) || !length(b)) next
  d <- mean(b) - mean(a)
  if (d > 0) hi <- hi + 1 else lo <- lo + 1
  log_msg("  %-24s %12.3f %12.3f %+11.3f  %s", s, mean(a), mean(b), d,
          if (d > 0) "challenged HIGHER" else "challenged lower")
}
log_msg("")
log_msg("  %d studies challenged higher, %d challenged lower", hi, lo)
if (hi + lo == 0) {
  log_msg("  no study has both a control and a challenged arm with sIgA")
} else if (hi > lo) {
  log_msg("  READING: the majority show sIgA RISING under challenge, which")
  log_msg("  contradicts the +1 orientation. Secretory IgA marks an active")
  log_msg("  mucosal immune response rather than a healthy gut - the same")
  log_msg("  reasoning already documented for systemic immunoglobulin.")
} else {
  log_msg("  READING: most studies follow the assumed direction, so the low AUC")
  log_msg("  is not explained by a wrong sign.")
}

# also against additive, which should move the opposite way if the sign is wrong
log_msg("")
log_msg("  Same test against ADDITIVE arms (a beneficial supplement):")
ha <- 0; la <- 0
for (s in unique(bio$Study_ID)) {
  k <- bio$Study_ID == s & is.finite(v)
  a <- v[k & bio$cls %in% HEALTHY]; b <- v[k & bio$cls %in% ADDITIVE]
  if (!length(a) || !length(b)) next
  if (mean(b) - mean(a) > 0) ha <- ha + 1 else la <- la + 1
}
log_msg("    %d studies additive higher, %d additive lower", ha, la)
log_msg("    If sIgA marked poor health, additives should LOWER it. If additives")
log_msg("    also raise it, sIgA is tracking immune activity generally rather")
log_msg("    than direction of health - which no single sign can capture.")

# =============================================================================
log_msg("")
log_msg("=========================================================")
log_msg("TEST 2 - Effect of flipping the sIgA sign")
log_msg("=========================================================")
sc <- gh[["GHI_dom_Immunoglobulin_Mucosal"]]
log_msg("  mucosal domain as oriented : AUC %.3f",
        auc_of(-sc[gh$cls %in% CHALLENGE], -sc[gh$cls %in% HEALTHY]))
log_msg("  mucosal domain flipped     : AUC %.3f",
        auc_of(sc[gh$cls %in% CHALLENGE], sc[gh$cls %in% HEALTHY]))

combo <- function(M) { M <- as.matrix(M); n <- rowSums(!is.na(M))
  ifelse(n >= 1, rowMeans(M, na.rm = TRUE), NA_real_) }
g2 <- gh; g2[["GHI_dom_Immunoglobulin_Mucosal"]] <- -sc
for (nm in c("current","sIgA flipped","mucosal removed")) {
  d <- switch(nm, "current" = gh[, dom, drop=FALSE],
                  "sIgA flipped" = g2[, dom, drop=FALSE],
                  gh[, setdiff(dom,"GHI_dom_Immunoglobulin_Mucosal"), drop=FALSE])
  x <- combo(d)
  log_msg("  combined index %-18s AUC %.3f", nm,
          auc_of(-x[gh$cls %in% CHALLENGE], -x[gh$cls %in% HEALTHY]))
}

# =============================================================================
log_msg("")
log_msg("=========================================================")
log_msg("TEST 3 - What if the two immunoglobulin domains are combined?")
log_msg("=========================================================")
log_msg("  Combining them assumes mucosal and systemic immunoglobulin carry the")
log_msg("  same meaning. If they point in opposite directions, averaging cancels")
log_msg("  both - the failure mode that keeps IL-6 and IL-10 apart.")
log_msg("")
muc <- gh[["GHI_dom_Immunoglobulin_Mucosal"]]
sys <- gh[["GHI_dom_Immunoglobulin_Systemic"]]
log_msg("  mucosal  AUC %.3f", auc_of(-muc[gh$cls %in% CHALLENGE], -muc[gh$cls %in% HEALTHY]))
log_msg("  systemic AUC %.3f", auc_of(-sys[gh$cls %in% CHALLENGE], -sys[gh$cls %in% HEALTHY]))
ok <- is.finite(muc) & is.finite(sys)
log_msg("  correlation between them: %s (n=%d)",
        ifelse(sum(ok) > 5, sprintf("%+.3f", cor(muc[ok], sys[ok], method="spearman")), "too few"),
        sum(ok))
merged <- combo(cbind(muc, sys))
log_msg("  merged   AUC %.3f", auc_of(-merged[gh$cls %in% CHALLENGE], -merged[gh$cls %in% HEALTHY]))
log_msg("")
d7 <- setdiff(dom, c("GHI_dom_Immunoglobulin_Mucosal","GHI_dom_Immunoglobulin_Systemic"))
x7 <- combo(cbind(gh[, d7, drop=FALSE], Immunoglobulin = merged))
x8 <- combo(gh[, dom, drop=FALSE])
log_msg("  combined index, 8 domains (split)  AUC %.3f",
        auc_of(-x8[gh$cls %in% CHALLENGE], -x8[gh$cls %in% HEALTHY]))
log_msg("  combined index, 7 domains (merged) AUC %.3f",
        auc_of(-x7[gh$cls %in% CHALLENGE], -x7[gh$cls %in% HEALTHY]))
log_msg("")
log_msg("  If the two domains sit on opposite sides of 0.5, merging them averages")
log_msg("  a marker of health with a marker of immune activation. The combined")
log_msg("  AUC may barely move while the domain stops meaning anything.")
log_msg("=== DONE ===")
