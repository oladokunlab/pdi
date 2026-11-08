#!/usr/bin/env Rscript
# =============================================================================
# 29_ghi_diagnostics.R
#
# Six diagnostics on the GHI, answering one question: is the modest
# discrimination (combined AUC, computed below) caused by a noisy index, or by genuine
# biological overlap between groups?
#
# RUN:  Rscript 29_ghi_diagnostics.R
# TIME: 2-5 minutes.
#
# -----------------------------------------------------------------------------
# WHY THIS ORDER
# -----------------------------------------------------------------------------
# Test 1 decides whether any of the rest is worth acting on. If the index
# discriminates well WITHIN a single study but poorly when studies are pooled,
# the construction is sound and cross-study heterogeneity is the problem -
# which is fixable. If even the best single study discriminates poorly, the
# overlap is biological and no amount of reweighting will remove it.
#
# Everything after that is either a specific improvement to test, or a check
# on whether a proposed practical application is defensible.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

HEALTHY   <- "Unchallenged control"
ADDITIVE  <- "Non-disruptive additive"
CHALLENGE <- c("Infectious challenge","Toxin challenge","Dietary challenge","Antibiotic challenge","Environmental challenge")
set.seed(42)

auc_of <- function(pos, neg) {
  pos <- pos[is.finite(pos)]; neg <- neg[is.finite(neg)]
  if (length(pos) < 3 || length(neg) < 3) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos)*(length(pos)+1)/2) / (length(pos)*length(neg))
}

log_msg("Loading GHI scores")
gh <- as.data.frame(fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")))
gh$Study_ID <- trimws(gh$Study_ID)
gh$cls <- trimws(as.character(gh$Intervention_Class))
dom_cols <- grep("^GHI_dom_", names(gh), value = TRUE)
log_msg("  %d rows, %d domains", nrow(gh), length(dom_cols))

# subtype, for test 3
gh$subtype <- trimws(as.character(gh$Intervention_Subclass))

SC <- "GHI_Combined_raw"

# =============================================================================
log_msg("")
log_msg("=========================================================")
log_msg("TEST 1 - Within-study AUC: is the index noisy, or do the groups overlap?")
log_msg("=========================================================")
log_msg("Pooled AUC mixes two things: how well the index separates groups, and")
log_msg("how much studies differ from each other. Computing AUC inside each")
log_msg("study separately removes the second.")
log_msg("")
# --- pooled AUC ---------------------------------------------------------
# Computed from the score table with auc_of() at every reporting site, so the
# printed value and the figure always come from the same calculation.
pooled_auc <- auc_of(-gh[[SC]][gh$cls %in% CHALLENGE], -gh[[SC]][gh$cls %in% HEALTHY])
fmt_pooled <- if (is.na(pooled_auc)) "n/a" else sprintf("%.3f", pooled_auc)
log_msg("  pooled AUC computed from the score table: %s", fmt_pooled)

ws <- list()
for (s in unique(gh$Study_ID)) {
  d <- gh[gh$Study_ID == s, ]
  a <- auc_of(-d[[SC]][d$cls %in% CHALLENGE], -d[[SC]][d$cls %in% HEALTHY])
  if (!is.na(a)) ws[[s]] <- a
}
if (length(ws)) {
  v <- unlist(ws)
  log_msg("  %d studies had >=3 control and >=3 challenged samples", length(v))
  log_msg("  within-study AUC: median %.3f, range %.3f - %.3f", median(v), min(v), max(v))
  log_msg("  studies above 0.80: %d of %d", sum(v >= 0.80), length(v))
  log_msg("  POOLED AUC for comparison: %s", fmt_pooled)
  log_msg("")
  for (s in names(sort(v, decreasing = TRUE))) log_msg("    %-24s %.3f", s, v[[s]])
  log_msg("")
  if (median(v) > 0.75) {
    log_msg("  READING: the index discriminates well within a study. The pooled")
    log_msg("  figure is depressed by cross-study heterogeneity, which is fixable.")
  } else {
    log_msg("  READING: discrimination is modest even within a single study.")
    log_msg("  The overlap is biological, not an artefact of pooling. Reweighting")
    log_msg("  and reconstruction will not remove it.")
  }
} else log_msg("  no study had enough of both groups")

# =============================================================================
log_msg("")
log_msg("=========================================================")
log_msg("TEST 2 - Do control birds score the same across studies?")
log_msg("=========================================================")
log_msg("Every validation so far is RELATIVE - challenged vs control inside one")
log_msg("study. Any application comparing a score against a stored baseline needs")
log_msg("the score to mean the same thing in absolute terms across studies.")
log_msg("")
ctl <- gh[gh$cls %in% HEALTHY & is.finite(gh[[SC]]), ]
cm <- tapply(ctl[[SC]], ctl$Study_ID, mean)
cm <- cm[is.finite(cm)]
all_sd <- stats::sd(ctl[[SC]])
between_sd <- stats::sd(cm)
within_sd <- sqrt(mean(tapply(ctl[[SC]], ctl$Study_ID, function(x)
  if (length(x) > 1) stats::var(x) else NA_real_), na.rm = TRUE))
log_msg("  %d control samples across %d studies", nrow(ctl), length(cm))
log_msg("  SD of all control scores          : %.3f", all_sd)
log_msg("  SD of per-study control MEANS     : %.3f  <- between-study spread", between_sd)
log_msg("  typical SD WITHIN a study         : %.3f", within_sd)
log_msg("  ratio between/within              : %.2f", between_sd/within_sd)
log_msg("")
if (between_sd > within_sd) {
  log_msg("  READING: control birds differ MORE between studies than within them.")
  log_msg("  An absolute score is not comparable across studies, so any application")
  log_msg("  relying on a stored baseline is not supported. Same-study comparison")
  log_msg("  only.")
} else {
  log_msg("  READING: control birds score similarly across studies. Absolute")
  log_msg("  interpretation has some support, though it should be confirmed on")
  log_msg("  new data before deployment.")
}

# =============================================================================
log_msg("")
log_msg("=========================================================")
log_msg("TEST 3 - Does pooling challenge types create the overlap?")
log_msg("=========================================================")
log_msg("A mild LPS injection and a severe Salmonella infection are both labelled")
log_msg("'challenged'. If they damage the gut by very different amounts, pooling")
log_msg("them puts mild-challenge birds among the controls and depresses AUC for")
log_msg("a reason that is not the index's fault.")
log_msg("")
ch <- gh[gh$cls %in% CHALLENGE & is.finite(gh[[SC]]), ]
ctlv <- gh[[SC]][gh$cls %in% HEALTHY & is.finite(gh[[SC]])]
by_sub <- table(ch$subtype)
by_sub <- by_sub[!is.na(names(by_sub)) & nzchar(names(by_sub))]
log_msg("  challenge subtypes present:")
for (nm in names(sort(by_sub, decreasing = TRUE)))
  log_msg("    %-38s %3d samples", nm, by_sub[[nm]])
log_msg("")
log_msg("  AUC of each subtype against ALL controls:")
res3 <- c()
for (nm in names(by_sub)) {
  if (by_sub[[nm]] < 5) next
  a <- auc_of(-ch[[SC]][ch$subtype == nm], -ctlv)
  if (!is.na(a)) { res3[nm] <- a; log_msg("    %-38s %.3f (n=%d)", nm, a, by_sub[[nm]]) }
}
log_msg("")
if (length(res3) >= 2) {
  log_msg("  spread across subtypes: %.3f to %.3f", min(res3), max(res3))
  if (max(res3) - min(res3) > 0.15) {
    log_msg("  READING: challenge types differ substantially in how far they move")
    log_msg("  the index. Pooling them into one 'challenged' group understates")
    log_msg("  performance against the more severe models. Worth reporting AUC")
    log_msg("  per challenge type rather than pooled.")
  } else {
    log_msg("  READING: challenge types move the index by similar amounts. Pooling")
    log_msg("  is not what limits discrimination.")
  }
}

# =============================================================================
log_msg("")
log_msg("=========================================================")
log_msg("TEST 4 - Is the composite less noisy than any single marker?")
log_msg("=========================================================")
log_msg("If so, a trial using the GHI as its primary outcome needs fewer birds")
log_msg("than one using a single biomarker. That is a direct, quantifiable")
log_msg("saving and does not depend on the index being diagnostic.")
log_msg("")
resid_sd <- function(v, g) {
  ok <- is.finite(v)
  if (sum(ok) < 10) return(NA_real_)
  v <- v[ok]; g <- g[ok]
  sqrt(mean(tapply(v, g, function(x) if (length(x) > 1) stats::var(x) else NA_real_), na.rm = TRUE))
}
cmb <- resid_sd(gh[[SC]], gh$Study_ID)
log_msg("  within-study SD, combined GHI : %.3f", cmb)
dom_sd <- c()
for (d in dom_cols) {
  s <- resid_sd(gh[[d]], gh$Study_ID)
  if (!is.na(s)) {
    dom_sd[sub("^GHI_dom_","",d)] <- s
    log_msg("  within-study SD, %-24s %.3f%s",
            sub("^GHI_dom_","",d), s, if (s > cmb) "   (noisier)" else "")
  }
}
log_msg("")
log_msg("  A smaller SD means less noise for the same effect, so fewer birds are")
log_msg("  needed to detect a difference of a given size.")

# =============================================================================
log_msg("")
log_msg("=========================================================")
log_msg("TEST 5 - Does requiring more domains per sample improve discrimination?")
log_msg("=========================================================")
log_msg("Proration lets a sample be scored from a single domain. A one-domain")
log_msg("score and a six-domain score currently carry equal weight.")
log_msg("")
ndom <- rowSums(!is.na(gh[, dom_cols, drop = FALSE]))
cov_tab <- data.frame(k = integer(0), auc = numeric(0),
                      n_ch = integer(0), n_ct = integer(0))
for (k in 1:5) {
  sel <- ndom >= k
  a <- auc_of(-gh[[SC]][sel & gh$cls %in% CHALLENGE], -gh[[SC]][sel & gh$cls %in% HEALTHY])
  n_ch <- sum(sel & gh$cls %in% CHALLENGE); n_ct <- sum(sel & gh$cls %in% HEALTHY)
  cov_tab <- rbind(cov_tab, data.frame(k = k, auc = a, n_ch = n_ch, n_ct = n_ct))
  log_msg("  >= %d domain(s): AUC %s   (%d challenged, %d control)",
          k, ifelse(is.na(a), "n/a", sprintf("%.3f", a)), n_ch, n_ct)
}
log_msg("")
log_msg("  If AUC rises with the minimum, low-coverage samples are adding noise")
log_msg("  and a floor is worth imposing - at the cost of sample count.")

# =============================================================================
log_msg("")
log_msg("=========================================================")
log_msg("TEST 6 - Would weighting the domains beat averaging them?")
log_msg("=========================================================")
log_msg("All seven domains currently contribute equally. The per-domain AUCs")
log_msg("reported above show whether any single domain beats the combined index; if")
log_msg("one does, equal weighting is diluting the informative domains.")
log_msg("")
log_msg("Weights are fitted LEAVE-ONE-STUDY-OUT: the weight applied to a study is")
log_msg("derived from the other studies only, never from itself. Fitting weights")
log_msg("on all data and evaluating on the same data would be circular.")
log_msg("")
dom_auc <- c()
for (d in dom_cols) {
  a <- auc_of(-gh[[d]][gh$cls %in% CHALLENGE], -gh[[d]][gh$cls %in% HEALTHY])
  if (!is.na(a)) dom_auc[d] <- a
}
log_msg("  per-domain AUC (pooled, for orientation):")
for (d in names(sort(dom_auc, decreasing = TRUE)))
  log_msg("    %-26s %.3f", sub("^GHI_dom_","",d), dom_auc[d])

wscore <- rep(NA_real_, nrow(gh))
for (s in unique(gh$Study_ID)) {
  tr <- gh$Study_ID != s; te <- gh$Study_ID == s
  w <- c()
  for (d in dom_cols) {
    a <- auc_of(-gh[[d]][tr & gh$cls %in% CHALLENGE], -gh[[d]][tr & gh$cls %in% HEALTHY])
    w[d] <- if (is.na(a)) 0 else max(0, a - 0.5)   # weight = discrimination above chance
  }
  if (sum(w) == 0) next
  w <- w / sum(w)
  M <- as.matrix(gh[te, dom_cols, drop = FALSE])
  num <- rowSums(sweep(M, 2, w, "*"), na.rm = TRUE)
  den <- rowSums(sweep(!is.na(M), 2, w, "*"), na.rm = TRUE)
  wscore[te] <- ifelse(den > 0, num/den, NA_real_)
}
a_w <- auc_of(-wscore[gh$cls %in% CHALLENGE], -wscore[gh$cls %in% HEALTHY])
log_msg("")
log_msg("  equal-weight combined AUC (current) : %s", fmt_pooled)
log_msg("  LOSO-weighted combined AUC          : %s",
        ifelse(is.na(a_w), "n/a", sprintf("%.3f", a_w)))
if (!is.na(a_w)) {
  if (a_w > 0.70) log_msg("  READING: weighting gives a material improvement and is worth adopting.")
  else if (!is.na(pooled_auc) && a_w > pooled_auc + 0.02) log_msg("  READING: a modest improvement. Worth reporting as a sensitivity analysis.")
  else log_msg("  READING: weighting does not help. Equal weighting is not the constraint.")
}

log_msg("")
log_msg("Writing summary")
out <- c("GHI DIAGNOSTICS", strrep("=", 74), "",
  sprintf("Within-study AUC: median %.3f (range %.3f-%.3f, %d studies)",
          if(length(ws)) median(unlist(ws)) else NA, if(length(ws)) min(unlist(ws)) else NA,
          if(length(ws)) max(unlist(ws)) else NA, length(ws)),
  sprintf("Pooled AUC: %s", fmt_pooled), "",
  sprintf("Control-score spread: between-study SD %.3f vs within-study SD %.3f (ratio %.2f)",
          between_sd, within_sd, between_sd/within_sd), "",
  sprintf("LOSO-weighted AUC: %s (equal-weight: %s)", 
          ifelse(is.na(a_w),"n/a",sprintf("%.3f", a_w)), fmt_pooled))

# Tests 3, 4 and 5 were computed above and printed to the log, but were not
# written to this file, so three of the six diagnostics had no citable source.
# They are appended here. Nothing is recomputed; these are the same values.
if (exists("res3") && length(res3)) {
  out <- c(out, "", "AUC by challenge subtype (each subtype against all controls):",
           sprintf("  %-38s %.3f", names(res3), as.numeric(res3)),
           sprintf("  spread across subtypes: %.3f to %.3f", min(res3), max(res3)))
}
if (exists("cmb") && is.finite(cmb)) {
  out <- c(out, "", "Within-study SD, composite against each single domain:",
           sprintf("  %-26s %.3f", "combined GHI", cmb))
  if (exists("dom_sd") && length(dom_sd)) {
    out <- c(out, sprintf("  %-26s %.3f%s", names(dom_sd), as.numeric(dom_sd),
                          ifelse(as.numeric(dom_sd) > cmb, "   (noisier)", "")))
  }
}
if (exists("cov_tab") && nrow(cov_tab)) {
  out <- c(out, "", "AUC by minimum domain coverage:",
           sprintf("  >= %d domain(s): AUC %s   (%d challenged, %d control)",
                   cov_tab$k, ifelse(is.na(cov_tab$auc), "n/a", sprintf("%.3f", cov_tab$auc)),
                   cov_tab$n_ch, cov_tab$n_ct))
}
writeLines(out, file.path(R_OUT_DIR, "ghi", "ghi_diagnostics.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "ghi", "ghi_diagnostics.txt"))
log_msg("=== DONE ===")
