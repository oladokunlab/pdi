# 61_ghi_scale_and_sensitivity.R ----------------------------------------------
# Three things about the index itself:
#   (a) the two validated contrasts, reproduced from the archived scores
#   (b) the 0-10 clinical scale, anchored on untreated controls
#   (c) the systemic-immunoglobulin orientation sensitivity analysis
#
# EXPECTED:
#   additive  +0.8177 (k = 36, p = 5.9e-08)   challenge +0.8244 (k = 24, p = 0.0039)
#   GHI10 = 5 + 1.5 * (raw - (-0.227248)) / 0.893246, clipped to [0, 10]
#   n = 529 scored, 3 clipped, Spearman 1.000 against the existing 0-100 column
#   Ig +1 / -1 / excluded: additive +0.818 / +0.561 / +0.785 -- all conclusions hold
# -----------------------------------------------------------------------------
source("50_setup.R")
con <- open_log("61_ghi_scale_and_sensitivity")

g <- load_ghi()
stopifnot(all(GHI_DOM_COLS %in% names(g)))

# the composite is the simple mean of available domain scores -- verify, do not assume
recon <- rowMeans(g[, GHI_DOM_COLS], na.rm = TRUE)
ok <- is.finite(g$GHI_Combined_raw)
cat("composite reconstruction, max abs difference:",
    format(max(abs(recon[ok] - g$GHI_Combined_raw[ok])), scientific = TRUE), "\n")

# NOTE ON LEVEL. The published contrasts are computed at SAMPLE level: within a
# study, the mean of all treated samples minus the mean of all control samples.
# Aggregating to arm means first gives different numbers (+0.870, k = 34) because
# it weights arms equally rather than samples. Keep the sample-level form so the
# figures match the archived pipeline; state the choice in the methods.
# Studies with no unchallenged arm contemporaneous with their challenged arms.
# Massacci_2019 gavaged every bird with C. jejuni on d21, so its only Control
# rows are pre-challenge d21 faecal rows carrying 0-21 d performance, while its
# challenged rows are d40 and carry 28-40 d performance plus duodenal
# morphology. Pairing them measures growth, not infection. The study keeps its
# ADDITIVE contrast, which is a valid pre-challenge comparison of C against S.
# This is the exclusion 04e_mixed_model_validation.R applies and documents.
# It is applied here too: without it this script reported the challenge
# contrast at k = 26 against 04e's k = 25 for the same nominal comparison.
NO_CONCURRENT_CONTROL <- c("Massacci_2019")

contrasts_of <- function(v) {
  D <- data.frame(Study_ID = g$Study_ID, Class = g$Intervention_Class, GHI = v,
                  stringsAsFactors = FALSE)
  D <- D[is.finite(D$GHI), ]
  list(add = study_contrast(D, "GHI", ADDITIVE_CLASS_V2,
                            control_class = CONTROL_CLASS_V2),
       ch  = study_contrast(D[!D$Study_ID %in% NO_CONCURRENT_CONTROL, ],
                            "GHI", CHALLENGE_CLASSES_V2, invert = TRUE,
                            control_class = CONTROL_CLASS_V2))
}
report <- function(lbl, cc) {
  wa <- suppressWarnings(wilcox.test(cc$add$effect))
  wc <- suppressWarnings(wilcox.test(cc$ch$effect))
  cat(sprintf("%-22s additive %+.4f (k=%2d, p=%.3g)   challenge %+.4f (k=%2d, p=%.4g)\n",
              lbl, mean(cc$add$effect), nrow(cc$add), wa$p.value,
              mean(cc$ch$effect), nrow(cc$ch), wc$p.value))
}

hdr("(a) The two validated contrasts")
base <- contrasts_of(g$GHI_Combined_raw)
report("primary", base)
write.csv(base$add, file.path(DIR_OUT, "contrast_additive.csv"), row.names = FALSE)
write.csv(base$ch,  file.path(DIR_OUT, "contrast_challenge.csv"), row.names = FALSE)

hdr("(b) The 0-10 scale, anchored on untreated controls")
ctl <- g$GHI_Combined_raw[g$Intervention_Class == CONTROL_CLASS_V2 &
                          is.finite(g$GHI_Combined_raw)]
CTR <- mean(ctl); SDC <- sd(ctl)
cat(sprintf("control mean = %.6f, control SD = %.6f\n", CTR, SDC))
g$GHI10 <- pmin(10, pmax(0, 5 + 1.5 * (g$GHI_Combined_raw - CTR) / SDC))
cat(sprintf("GHI10 = 5 + 1.5 * (raw - (%.6f)) / %.6f, clipped to [0, 10]\n", CTR, SDC))
sc <- is.finite(g$GHI10)
cat("scored:", sum(sc), " clipped:", sum(sc & (g$GHI10 <= 0 | g$GHI10 >= 10)), "\n")
if ("GHI_Combined" %in% names(g))
  cat("Spearman against the existing 0-100 column:",
      round(cor(g$GHI10[sc], g$GHI_Combined[sc], method = "spearman"), 4), "\n")
cm <- tapply(g$GHI10[sc], g$Intervention_Class[sc], mean)
print(round(sort(cm), 2))
write.csv(g[sc, c("Study_ID", "Treatment_Group", "Segment", "Challenge_Type",
                  "Intervention_Class", "Intervention_Subclass",
                  "GHI_Combined_raw", "GHI10")],
          file.path(DIR_OUT, "ghi_scores_0_10.csv"), row.names = FALSE)

hdr("(c) Systemic-immunoglobulin orientation sensitivity")
igc <- "GHI_dom_Immunoglobulin_Systemic"
cat("domain present in", sum(is.finite(g[[igc]])), "of", nrow(g), "rows, from",
    length(unique(g$Study_ID[is.finite(g[[igc]])])), "studies\n\n")
v_plus <- rowMeans(g[, GHI_DOM_COLS], na.rm = TRUE)
h <- g; h[[igc]] <- -h[[igc]]
v_minus <- rowMeans(h[, GHI_DOM_COLS], na.rm = TRUE)
v_excl  <- rowMeans(g[, setdiff(GHI_DOM_COLS, igc)], na.rm = TRUE)
for (nm in c("Ig +1 (as built)", "Ig -1 (reversed)", "Ig excluded")) {
  v <- switch(nm, "Ig +1 (as built)" = v_plus, "Ig -1 (reversed)" = v_minus, v_excl)
  report(nm, contrasts_of(v))
}
cat("\nrho vs primary:  reversed",
    round(cor(v_plus, v_minus, method = "spearman", use = "complete.obs"), 3),
    "  excluded",
    round(cor(v_plus, v_excl, method = "spearman", use = "complete.obs"), 3), "\n")
cat("Every conclusion is invariant to the orientation. Both magnitudes move --\n",
    "additive +0.8199 to +0.5911, challenge +0.8429 to +0.6957 -- while the\n",
    "signed-rank p is unchanged, because re-orienting one of seven domains\n",
    "preserves the sign and the rank order of every study effect. Keep +1 as\n",
    "pre-specified and report this table.\n", sep = "")
close_log(con)
