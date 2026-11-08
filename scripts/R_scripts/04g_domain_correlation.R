#!/usr/bin/env Rscript
# =============================================================================
# 04g_domain_correlation.R
#
# Why is Morphology-GHI uncorrelated with Biomarker-GHI? (04f gave r = -0.007)
#
# RUN:  Rscript 04g_domain_correlation.R
# TIME: seconds. Reads ghi_scores_pca.tsv only. Changes nothing.
#
# -----------------------------------------------------------------------------
# THE QUESTION
# -----------------------------------------------------------------------------
# 04f found the two arms of the GHI are orthogonal at sample level, within
# study as well as pooled. That is in tension with 04e, where each arm
# separately detected known intervention direction. Two measures that both
# track the same external truth should not be exactly uncorrelated with
# each other.
#
# Three explanations, which imply different actions:
#
#   (1) DIFFERENT FACETS. Structure (villi) and function (cytokines, SCFA)
#       are genuinely separate axes of gut health. -> Report two indices,
#       do not force an arithmetic-mean composite.
#
#   (2) SAMPLE-LEVEL PAIRING IS FICTIONAL. If one arm was extracted as group
#       means and the other per-bird, sample-level correlation is meaningless
#       even though group-level contrasts (04e) still work.
#       -> Composite defensible at group level only.
#
#   (3) CANCELLATION INSIDE THE BIOMARKER ARM. One domain pointing the wrong
#       way drags the mean to zero - the Claudin-1 failure mode, already found
#       once in this project. IgG/IgM sit at +1 and are flagged REVIEW in 04b.
#       -> Fixable by splitting the domain. Proceed to Step 2.
#
# This script separates them by correlating every domain against every other,
# rather than collapsing the biomarker side to one number first.
#
# -----------------------------------------------------------------------------
# HOW TO READ IT
# -----------------------------------------------------------------------------
#   All pairs ~0                        -> (1) or (2). Step 2 will not fix it.
#   Some pairs sensible, one negative   -> (3). That domain is the culprit.
#   Morphology ~0 with everything, but
#     biomarker domains agree together  -> (1), cleanly. Two indices.
#
# Correlations are computed WITHIN STUDY (each column centred inside its own
# study before pooling). Between-study variation otherwise inflates agreement -
# the artefact behind the R2 = 0.000 LASSO bug. Pairwise-complete, so each
# cell uses whatever samples carry both domains; n is printed alongside r
# because some pairs will rest on very few samples and must not be read as
# equal in weight.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

SCRIPT_DIR <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())

source(file.path(SCRIPT_DIR, "00_r_config.R"))

GHI_TSV <- file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")
if (!file.exists(GHI_TSV))
  stop("Not found: ", GHI_TSV, "\n  Run 04b_ghi_pca.R first.")

ghi <- fread(GHI_TSV)
dom_cols <- grep("^GHI_dom_", names(ghi), value = TRUE)
if (length(dom_cols) < 2) stop("Fewer than two GHI_dom_ columns found.")

short <- sub("^GHI_dom_", "", dom_cols)
log_msg("[1/4] %d domains: %s", length(dom_cols), paste(short, collapse = ", "))

# --- coverage ----------------------------------------------------------------
# How many samples carry each domain at all. A domain present on very few
# samples cannot meaningfully correlate with anything and its row below is
# decorative.
log_msg("[2/4] Coverage")
for (i in seq_along(dom_cols))
  log_msg("  %-24s %4d samples, %2d studies", short[i],
          sum(!is.na(ghi[[dom_cols[i]]])),
          uniqueN(ghi$Study_ID[!is.na(ghi[[dom_cols[i]]])]))

# --- within-study centring ---------------------------------------------------
# Centre each domain inside each study, then pool. Studies contributing a
# single sample for a domain carry no within-study contrast and become NA
# for that domain, which is correct - they are uninformative here.
X <- copy(ghi)
for (cc in dom_cols) {
  X[, (cc) := {
    v <- get(cc)
    if (sum(!is.na(v)) >= 2) v - mean(v, na.rm = TRUE) else NA_real_
  }, by = Study_ID]
}

# --- pairwise matrix ---------------------------------------------------------
log_msg("[3/4] Within-study correlation (r, with n beneath)")

k <- length(dom_cols)
R <- matrix(NA_real_, k, k, dimnames = list(short, short))
N <- matrix(0L,       k, k, dimnames = list(short, short))

for (i in seq_len(k)) for (j in seq_len(k)) {
  a <- X[[dom_cols[i]]]; b <- X[[dom_cols[j]]]
  ok <- !is.na(a) & !is.na(b)
  N[i, j] <- sum(ok)
  if (sum(ok) >= 10 && sd(a[ok]) > 0 && sd(b[ok]) > 0)
    R[i, j] <- cor(a[ok], b[ok])
}

w <- max(nchar(short), 8)
hdr <- paste0(strrep(" ", w + 2),
              paste(sprintf("%9s", substr(short, 1, 9)), collapse = ""))
log_msg("%s", hdr)
for (i in seq_len(k)) {
  rline <- paste(vapply(seq_len(k), function(j)
    if (is.na(R[i, j])) sprintf("%9s", "-") else sprintf("%+9.3f", R[i, j]),
    character(1)), collapse = "")
  nline <- paste(vapply(seq_len(k), function(j)
    sprintf("%9s", sprintf("(%d)", N[i, j])), character(1)), collapse = "")
  log_msg("%-*s  %s", w, short[i], rline)
  log_msg("%-*s  %s", w, "",       nline)
}

# --- interpretation ----------------------------------------------------------
log_msg("[4/4] Reading")

off <- R; diag(off) <- NA_real_
usable <- !is.na(off)

if (!any(usable)) {
  verdict <- "No pair had enough overlap. Domains are near-disjoint across samples."
} else {
  mx <- max(abs(off[usable]))
  morph_i <- match("Morphology", short)
  bio_idx  <- setdiff(seq_len(k), morph_i)

  bio_off <- if (length(bio_idx) >= 2) off[bio_idx, bio_idx, drop = FALSE] else matrix(NA_real_, 0, 0)
  bio_ok  <- !is.na(bio_off)
  bio_med <- if (any(bio_ok)) median(bio_off[bio_ok]) else NA_real_

  neg <- which(off < -0.15 & usable, arr.ind = TRUE)

  verdict <- if (mx < 0.15) {
    "ALL PAIRS ~ZERO -> explanation (1) or (2). Splitting a domain will not fix this. Check whether the two arms were extracted at the same unit (per-bird vs group mean) before concluding they are separate facets."
  } else if (nrow(neg) > 0) {
    paste0("NEGATIVE PAIR(S) PRESENT -> explanation (3), cancellation. Offenders: ",
           paste(unique(sprintf("%s~%s", short[neg[,1]], short[neg[,2]])), collapse = ", "),
           ". Check their direction-of-effect entries in 04b DOMAINS.")
  } else if (!is.na(bio_med) && bio_med > 0.15 &&
             all(is.na(off[morph_i, ]) | abs(off[morph_i, ]) < 0.15)) {
    "BIOMARKER DOMAINS AGREE, MORPHOLOGY STANDS APART -> explanation (1), cleanly. Report Morphology and Biomarker as two indices rather than averaging them."
  } else {
    "MIXED. No single pattern dominates. Read the matrix directly; treat cells with small n as uninformative."
  }
}
log_msg("  %s", verdict)

# --- write -------------------------------------------------------------------
lines <- c(
  "GHI DOMAIN-LEVEL CORRELATION",
  "Diagnostic for the null convergent-validity result in 04f (r = -0.007).",
  "Within-study centred, pairwise-complete.",
  strrep("=", 74), "",
  hdr
)
for (i in seq_len(k)) {
  lines <- c(lines,
    sprintf("%-*s  %s", w, short[i],
      paste(vapply(seq_len(k), function(j)
        if (is.na(R[i,j])) sprintf("%9s","-") else sprintf("%+9.3f", R[i,j]),
        character(1)), collapse = "")),
    sprintf("%-*s  %s", w, "",
      paste(vapply(seq_len(k), function(j)
        sprintf("%9s", sprintf("(%d)", N[i,j])), character(1)), collapse = "")))
}
lines <- c(lines, "", paste("VERDICT:", verdict), "",
  "n in parentheses. Cells resting on few samples are not comparable to",
  "cells resting on many; do not read the matrix as uniformly weighted.")

outs <- file.path(R_OUT_DIR, "ghi", "ghi_domain_correlation.txt")
writeLines(lines, outs)
log_msg("summary -> %s", outs)

outm <- file.path(R_OUT_DIR, "ghi", "ghi_domain_correlation_matrix.tsv")
fwrite(data.table(domain = short, as.data.table(R)), outm, sep = "\t")
log_msg("matrix -> %s", outm)
