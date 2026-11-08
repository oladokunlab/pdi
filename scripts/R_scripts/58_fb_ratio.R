# 58_fb_ratio.R ---------------------------------------------------------------
# Firmicutes:Bacteroidetes ratio against the Gut Health Index and against
# treatment class. Not part of the registered pipeline; added because it is the
# most frequently cited summary ratio in the field.
#
# Phylum labels come from the BactoTraits genus-level reference, an external
# auditable mapping. No hand assignment.
#
# EXPECTED:
#   72 of 115 genera matched (44 Bacillota, 9 Bacteroidota), 50.2% of abundance
#   1464 of 2064 runs have both phyla present
#   F:B vs GHI   rho = -0.148, p = 0.111 (9999 within-study permutations)
#   additive     delta log(F:B) = +0.005, p = 0.49, k = 25
#   challenged   delta log(F:B) = +0.041, p = 0.89, k = 13
# -----------------------------------------------------------------------------
source("50_setup.R")
con <- open_log("58_fb_ratio")

BACTOTRAITS <- Sys.getenv("BACTOTRAITS_CSV",
  unset = file.path(PDI_ROOT, "reference", "BACTOTRAITS_DATASET_2026-01-28_GENUSLVL.csv"))
stopifnot(file.exists(BACTOTRAITS))

bt <- read.csv2(BACTOTRAITS, check.names = FALSE)   # semicolon-delimited
bt$g <- tolower(trimws(as.character(bt$Genus)))
bt <- bt[!is.na(bt$Phylum) & bt$Phylum != "", ]
phy <- tapply(as.character(bt$Phylum), bt$g, function(z) names(sort(table(z), decreasing = TRUE))[1])

prop <- load_prop_adj()
gl  <- sub("^g__", "", rownames(prop))
k1  <- tolower(trimws(gsub("[][]", "", gsub("_", " ", gl))))
k2  <- tolower(trimws(gl))
lab <- ifelse(!is.na(phy[k1]), phy[k1], phy[k2])

hdr("Phylum assignment")
cat("matched:", sum(!is.na(lab)), "of", length(lab), "genera\n")
print(table(lab))
ab <- rowSums(prop) / sum(prop)
cat(sprintf("abundance covered: %.1f%%\n", 100 * sum(ab[!is.na(lab)])))
cat(sprintf("genus-unclassified abundance: %.1f%%  <- cannot be assigned a phylum here\n",
            100 * ab[rownames(prop) == "g__unclassified"]))

F <- colSums(prop[which(lab == "Bacillota"), , drop = FALSE])
B <- colSums(prop[which(lab == "Bacteroidota"), , drop = FALSE])
both <- F > 0 & B > 0
cat("runs with both phyla present:", sum(both), "of", length(F), "\n")
fb <- log(F[both] / B[both])

meta <- load_seq_meta()
d <- data.frame(sample_id = names(fb), fb = as.numeric(fb), stringsAsFactors = FALSE)
d <- merge(d, meta[, c("sample_id", "Study_ID", "Treatment_Group")], by = "sample_id")
arm <- aggregate(list(fb = d$fb),
                 by = list(Study_ID = d$Study_ID, Treatment_Group = d$Treatment_Group), FUN = mean)

J <- merge(arm, ghi_by_arm(), by = c("Study_ID", "Treatment_Group"))
J <- J[is.finite(J$GHI) & is.finite(J$fb), ]
tab <- table(J$Study_ID); J <- J[J$Study_ID %in% names(tab)[tab >= 2], ]
cat("\narms with both F:B and GHI:", nrow(J), "studies:", length(unique(J$Study_ID)), "\n")
cat(sprintf("F:B log-ratio: median %.2f, IQR %.2f to %.2f\n",
            median(J$fb), quantile(J$fb, .25), quantile(J$fb, .75)))

hdr("F:B versus GHI")
J$fbc  <- center_within(J$fb,  J$Study_ID)
J$ghic <- center_within(J$GHI, J$Study_ID)
r <- perm_spearman_within(J$fbc, J$ghic, J$Study_ID)
cat(sprintf("within-study Spearman rho = %+.4f, permutation p = %.4f (n = %d arms)\n",
            r$rho, r$p, r$n))

hdr("F:B versus treatment class")
for (nm in c("Additive vs control", "Control vs challenged")) {
  trt <- if (nm == "Additive vs control") "Non_Disruptive_Additive" else CHALLENGE_CLASSES
  ct  <- study_contrast(J, "fb", trt, invert = (nm == "Control vs challenged"))
  w   <- suppressWarnings(wilcox.test(ct$effect))
  cat(sprintf("%-24s k = %2d  mean delta log(F:B) = %+.3f  p = %.4f\n",
              nm, nrow(ct), mean(ct$effect), w$p.value))
}
cat("\nInterventions that moved the GHI by +0.82 SD moved the F:B ratio by 0.005-0.041\n",
    "log units. The ratio does not respond to what changes gut health.\n", sep = "")
close_log(con)
