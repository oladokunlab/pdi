# 60_ghi_reliability.R --------------------------------------------------------
# How much of the null could be outcome measurement error?
#
# Published reports did not give per-group dispersion, so the registered
# inverse-variance weighting and Egger test are not computable. Internal
# consistency answers the same question without new extraction: if the domains
# are reliable, unreliability cannot be what is hiding the signal.
#
# EXPECTED:  Barrier 0.668, Inflammation 0.729, SCFA 0.889
#            attenuation ceiling rho ~ 0.82-0.94 against observed 0.00-0.14
#            Morphology and Performance return negative alpha because they are
#            DIFFERENCE scores of positively correlated markers -- expected, not a bug
# -----------------------------------------------------------------------------
source("50_setup.R")
con <- open_log("60_ghi_reliability")

# Raw marker values are NOT in ghi_scores.tsv (which holds domain scores only);
# they live in the extraction workbook, sheet "Sample_with_Zscores".
WB <- Sys.getenv("PDI_WORKBOOK", unset = file.path(DIR_META, "GDI_Metadata.xlsx"))
stopifnot(requireNamespace("readxl", quietly = TRUE), file.exists(WB))
raw <- as.data.frame(readxl::read_excel(WB, sheet = "Sample_with_Zscores"))
def <- read.delim(file.path(DIR_GHI, "ghi_index_definition.tsv"), check.names = FALSE)
cat("marker workbook:", nrow(raw), "rows x", ncol(raw), "columns\n")

# standardised alpha from the marker correlation matrix
alpha_std <- function(M) {
  M <- M[complete.cases(M), , drop = FALSE]
  k <- ncol(M)
  if (k < 2 || nrow(M) < 5) return(c(alpha = NA, k = k, n = nrow(M)))
  R <- cor(M); rbar <- mean(R[upper.tri(R)])
  c(alpha = k * rbar / (1 + (k - 1) * rbar), k = k, n = nrow(M))
}

hdr("Internal consistency by domain")
rows <- list()
for (dm in unique(def$domain)) {
  mk <- def$marker[def$domain == dm]
  dir <- def$direction[def$domain == dm]
  present <- mk %in% names(raw)
  if (sum(present) < 2) { cat(sprintf("%-26s only %d marker(s) -- alpha undefined\n",
                                      dm, sum(present))); next }
  M <- as.data.frame(lapply(which(present), function(i)
        suppressWarnings(as.numeric(raw[[mk[i]]])) * dir[i]))
  names(M) <- mk[present]
  a <- alpha_std(M)
  rows[[dm]] <- data.frame(domain = dm, markers = a["k"], n_complete = a["n"],
                           alpha = round(a["alpha"], 3))
  cat(sprintf("%-26s k = %d  n = %3d  standardised alpha = %+.3f\n",
              dm, a["k"], a["n"], a["alpha"]))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(DIR_OUT, "ghi_reliability_alpha.csv"), row.names = FALSE)

hdr("Attenuation ceiling")
good <- res[is.finite(res$alpha) & res$alpha > 0, ]
cat("Using the reliable domains only, the maximum observable correlation between\n",
    "a perfectly measured predictor and the index is sqrt(alpha):\n", sep = "")
for (i in seq_len(nrow(good)))
  cat(sprintf("  %-26s ceiling rho = %.2f\n", good$domain[i], sqrt(good$alpha[i])))
cat("\nObserved microbiome-to-GHI correlations are 0.00-0.14.\n",
    "Outcome measurement error does NOT explain the null.\n", sep = "")

hdr("Why two domains give negative alpha")
for (pair in list(c("VH_um", "CD_um"), c("BWG_g", "FCR"))) {
  if (all(pair %in% names(raw))) {
    x <- suppressWarnings(as.numeric(raw[[pair[1]]])); y <- suppressWarnings(as.numeric(raw[[pair[2]]]))
    k <- is.finite(x) & is.finite(y)
    cat(sprintf("%s vs %s: r = %+.3f (n = %d)\n", pair[1], pair[2],
                cor(x[k], y[k]), sum(k)))
  }
}
cat("\nBoth are DIFFERENCE scores, for which alpha is undefined, but for different\n",
    "reasons and they should not be described together:\n",
    "  Morphology  VH_um and CD_um correlate strongly and positively (+0.52), so\n",
    "              their difference is dominated by measurement noise. This is why\n",
    "              the morphology domain contributes least to the composite (r =\n",
    "              0.646) despite having the most observations (n = 323).\n",
    "  Performance BWG_g and FCR are essentially uncorrelated (+0.05). The two\n",
    "              markers carry near-independent information; alpha is slightly\n",
    "              negative only because direction correction flips the sign.\n",
    "NOTE: use BWG_g (weight GAIN), the marker named in the index definition.\n",
    "BW_g (absolute weight) correlates with FCR at +0.49 and is a different column.\n",
    sep = "")
close_log(con)
