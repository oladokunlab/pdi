#!/usr/bin/env Rscript
# =============================================================================
# score_new_ghi.R
#
# Computes the Gut Health Index for a study outside the derivation corpus.
#
# Usage:  Rscript score_new_ghi.R <biomarkers.tsv> [out.tsv]
#
# Input: a tab-separated table, one row per treatment group or per sample.
#   Study_ID          required; markers are standardised within each value
#   Treatment_Group   optional; carried through to the output
#   marker columns    any subset of the panel below, numeric
#
# Panel and direction of effect (read at run time from ghi_index_definition.tsv;
# this list is the definition shipped with the deposit):
#   Morphology                VH_um (+)   CD_um (-)
#   Barrier                   ZO1 (+)     Ocln (+)    Cldn1 (+)   MUC2 (+)
#   Inflammation              IL6 (-)     TNFa (-)    IL1b (-)
#   Anti_Inflammatory         IL10 (+)
#   SCFA                      Ace (+)     Prop (+)    Buty (+)
#   Performance               BWG_g (+)   FCR (-)
#   Immunoglobulin_Systemic   IgG (+)     IgM (+)     IgA_serum (+)
#
# A negative direction indicates that higher raw values correspond to poorer
# gut health; those markers are inverted before use. Absent markers are
# accommodated: each domain is scored from the markers present, and the index
# is the mean of the domains scored. Where a domain was fitted by PCA and a row
# lacks one of its markers, the row is scored from the markers it carries and
# placed on the PC1 scale using the calibration constants in the definition
# file. Absent domains are expected, as the biomarkers measured are determined
# by study design.
#
# NOT SCORED, DELIBERATELY. Secretory IgA (sIgA) is excluded: it rises both in
# a healthy, well-regulated mucosa and in a bird mounting a defence, so neither
# sign is correct for a domain scored on its own. IFN-gamma is excluded on
# estimation grounds: adding it to Inflammation cuts complete cases for the
# principal-component fit from 85 to 29 while leaving the validated contrasts
# unchanged. Both were tested; see the amendments record.
#
# -----------------------------------------------------------------------------
# INTERPRETATION
# -----------------------------------------------------------------------------
# Markers are z-scored within study, so scores are relative to the baseline of
# the data supplied rather than to an external reference population. Higher
# values indicate better gut health. Zero corresponds to the mean of the
# supplied data, and SCORES FROM SEPARATE RUNS ARE NOT COMPARABLE.
#
# Two rows per study are the minimum for the within-study standardisation. The
# script exits if fewer are supplied. Scoring one treatment group alone is not
# possible and not meaningful: the index measures a difference from a baseline
# that must be present in the same input.
#
# -----------------------------------------------------------------------------
# VALIDATION (derivation corpus, 50 studies)
# -----------------------------------------------------------------------------
# Paired within-study differences on GHI_Combined_raw, pooled across studies:
#
#   beneficial additive against control   +0.818 SD, k = 36, p = 2.3e-07
#   control against challenged            +0.824 SD, k = 24, p = 0.018
#
# The index separates treatment from control in both directions. Its domains
# are only loosely correlated with one another, so it summarises several
# weakly coupled axes rather than one latent construct; read the domain scores
# alongside the composite, not instead of it.
#
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  cat("Usage: Rscript score_new_ghi.R <biomarkers.tsv> [out.tsv]\n"); quit(status = 1)
}
inp <- args[1]
outp <- if (length(args) > 1) args[2] else sub("\\.tsv$", "_ghi.tsv", inp)
if (!file.exists(inp)) stop("Not found: ", inp)
suppressPackageStartupMessages(library(data.table))

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
SD <- .locate()

cand <- c(file.path(SD, "ghi_index_definition.tsv"),
          file.path(SD, "..", "r_analysis", "ghi", "ghi_index_definition.tsv"),
          file.path(SD, "..", "..", "r_analysis", "ghi", "ghi_index_definition.tsv"))
dfn <- cand[file.exists(cand)][1]
if (is.na(dfn))
  stop("ghi_index_definition.tsv not found. Run 04b_ghi_pca.R, or place the\n",
       "  file beside this script.")
defn <- fread(dfn)
cat(sprintf("Definition: %d markers across %d domains (%s)\n",
            nrow(defn), uniqueN(defn$domain), basename(dfn)))

d <- fread(inp)
if (!"Study_ID" %in% names(d)) stop("Input needs a Study_ID column.")
present <- intersect(defn$marker, names(d))
if (!length(present)) stop("None of the panel markers are present in the input.")
cat(sprintf("Markers found: %d of %d (%s)\n", length(present), uniqueN(defn$marker),
            paste(present, collapse = ", ")))
absent <- setdiff(unique(defn$marker), present)
if (length(absent)) cat("Absent from input:", paste(absent, collapse = ", "), "\n")

small <- d[, .N, by = Study_ID][N < 2]
if (nrow(small))
  stop("Fewer than 2 rows for: ", paste(small$Study_ID, collapse = ", "),
       "\n  Within-study standardisation requires a minimum of 2 rows per study.")

# ---- z-score within study, then orient --------------------------------------
Z <- list()
for (m in present) {
  v <- suppressWarnings(as.numeric(d[[m]]))
  z <- rep(NA_real_, length(v))
  for (s in unique(d$Study_ID)) {
    i <- d$Study_ID == s
    mu <- mean(v[i], na.rm = TRUE); sdv <- stats::sd(v[i], na.rm = TRUE)
    if (is.finite(sdv) && sdv > 0) z[i] <- (v[i] - mu) / sdv
  }
  dir <- defn$direction[match(m, defn$marker)]
  Z[[m]] <- z * dir
}

# ---- domain scores -----------------------------------------------------------
# PCA loadings are applied where the domain was fitted by PCA and every one of
# its markers is present. Otherwise the domain is the mean of whichever
# oriented markers are available, which is how sparse domains were handled in
# the original construction.
dom_scores <- list()
for (dom in unique(defn$domain)) {
  sub <- defn[domain == dom]
  mk <- intersect(sub$marker, present)
  if (!length(mk)) next
  M <- do.call(cbind, Z[mk]); if (is.null(dim(M))) M <- matrix(M, ncol = 1)
  colnames(M) <- mk

  use_pca <- sub$method[1] == "PCA" && all(!is.na(sub$pc1_loading)) &&
             length(mk) == nrow(sub)
  n_av <- rowSums(!is.na(M))
  if (use_pca) {
    w <- sub$pc1_loading[match(mk, sub$marker)]
    sc <- as.numeric(M %*% w)
    cc <- stats::complete.cases(M)
    sc[!cc] <- NA_real_
    how <- "PC1"
    # Rows missing a marker of this domain cannot take a PC1 score. They are
    # scored from the markers they do carry and placed on the PC1 scale using
    # the calibration constants estimated during construction.
    off <- sub$fill_offset[1]; scl <- sub$fill_scale[1]
    if (!is.na(off) && !is.na(scl)) {
      need <- !cc & n_av >= 1
      if (any(need)) {
        sc[need] <- off + scl * rowMeans(M, na.rm = TRUE)[need]
        how <- sprintf("PC1+%d", sum(need))
      }
    }
  } else {
    sc <- ifelse(n_av >= 1, rowMeans(M, na.rm = TRUE), NA_real_)
    how <- "mean"
  }
  sc[!is.finite(sc)] <- NA_real_
  dom_scores[[dom]] <- sc
  note <- if (sum(!is.na(sc)) == 0 && how == "PC1")
            "  (no row carries every marker of this domain)" else ""
  cat(sprintf("  %-24s %-8s %3d scored%s\n", dom, how, sum(!is.na(sc)), note))
}
if (!length(dom_scores)) stop("No domain could be scored.")

# ---- combine -----------------------------------------------------------------
DM <- do.call(cbind, dom_scores)
if (is.null(dim(DM))) DM <- matrix(DM, ncol = 1)
nd <- rowSums(!is.na(DM))
combined <- ifelse(nd >= 1, rowMeans(DM, na.rm = TRUE), NA_real_)

morph <- if ("Morphology" %in% names(dom_scores)) dom_scores[["Morphology"]] else NA_real_
biod <- setdiff(names(dom_scores), "Morphology")
biom <- if (length(biod)) {
  B <- do.call(cbind, dom_scores[biod]); if (is.null(dim(B))) B <- matrix(B, ncol = 1)
  ifelse(rowSums(!is.na(B)) >= 1, rowMeans(B, na.rm = TRUE), NA_real_)
} else rep(NA_real_, nrow(d))

# Assembled in one call rather than by successive column assignment, which
# triggers a shallow-copy warning from data.table.
id_cols <- intersect(c("Study_ID", "Sample_ID", "Treatment_Group", "Segment"), names(d))
res <- as.data.table(c(
  setNames(lapply(id_cols, function(cc) d[[cc]]), id_cols),
  list(GHI_Combined = combined, GHI_Morphology = morph,
       GHI_Biomarker = biom, n_domains = nd),
  setNames(dom_scores, paste0("dom_", names(dom_scores)))))

cat("\nGut Health Index. Higher values indicate better gut health.\n")
lab_cols <- intersect(c("Study_ID", "Treatment_Group", "Segment"), names(res))
for (i in seq_len(nrow(res)))
  cat(sprintf("  %-46s %s  (%d domain%s)\n",
              paste(vapply(lab_cols, function(cc) as.character(res[[cc]][i]),
                           character(1)), collapse = " / "),
              if (is.na(res$GHI_Combined[i])) "     n/a" else sprintf("%+8.3f", res$GHI_Combined[i]),
              res$n_domains[i], if (res$n_domains[i] == 1) "" else "s"))
n_one <- sum(res$n_domains == 1, na.rm = TRUE)
if (n_one > 0)
  cat(sprintf("\n%d row(s) scored on a single domain. A single-domain score\n", n_one),
      "reflects that one biological aspect only.\n", sep = "")

cat("\nScores are standardised within study and are relative to the data\n")
cat("supplied. They are not comparable across separate runs.\n")

fwrite(res, outp, sep = "\t")
cat(sprintf("\nwritten -> %s\n", outp))
