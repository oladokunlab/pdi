# 63_rob_and_publication_bias.R -----------------------------------------------
# Registered sensitivity analyses that had not been run, plus the substitute
# for the registered Egger test.
#
# EXPECTED:
#   two SYRCLE domains have zero variance (ROB_5 Low x50, ROB_6 Unclear x50)
#   7 High ratings across 7 distinct studies, none with more than one
#   excluding high risk : additive +0.849 (k=31)  challenge +0.805 (k=21)
#   Illumina only       : additive +0.775 (k=33)  challenge +0.869 (k=23)
#   publication bias: additive clean; challenge asymmetry from Lucke_2018 alone
#                     (p = 0.017 with, p = 0.58 without)
#   The n >= 6 per group sensitivity is NOT executable: 'Sample size' is free text.
# -----------------------------------------------------------------------------
source("50_setup.R")
con <- open_log("63_rob_and_publication_bias")

WB <- Sys.getenv("PDI_WORKBOOK", unset = file.path(DIR_META, "GDI_Metadata.xlsx"))
have_xl <- requireNamespace("readxl", quietly = TRUE) && file.exists(WB)

if (have_xl) {
  sheets <- readxl::excel_sheets(WB)
  rb <- grep("risk|rob", sheets, ignore.case = TRUE, value = TRUE)
  hdr("Risk-of-bias domain distribution")
  if (length(rb)) {
    R <- as.data.frame(readxl::read_excel(WB, sheet = rb[1]))
    dom <- grep("^ROB", names(R), value = TRUE)
    tabs <- sapply(dom, function(d) table(factor(R[[d]], levels = c("Low", "Unclear", "High"))))
    print(t(tabs))
    dead <- dom[apply(tabs, 2, function(z) sum(z > 0) == 1)]
    cat("\nzero-variance domains (assess nothing):", paste(dead, collapse = ", "), "\n")
    high <- apply(R[, dom] == "High", 1, sum, na.rm = TRUE)
    idcol <- grep("study", names(R), ignore.case = TRUE, value = TRUE)[1]
    cat("studies with >=1 High:", sum(high > 0), " High ratings total:", sum(high),
        " max per study:", max(high), "\n")
    HIGH_STUDIES <- R[[idcol]][high > 0]
    cat(paste(" ", HIGH_STUDIES, collapse = "\n"), "\n")
    write.csv(data.frame(study = R[[idcol]], n_high = high),
              file.path(DIR_OUT, "rob_severity.csv"), row.names = FALSE)
  }
} else {
  cat("NOTE: readxl or the workbook is unavailable; supply PDI_WORKBOOK to run\n",
      "the risk-of-bias section. Set HIGH_STUDIES manually to continue.\n", sep = "")
  HIGH_STUDIES <- c("Gao_2024", "Han_2023", "Miska_2024", "Zhang_2018",
                    "Zhang_2023", "Zhao_2022", "Zhao_2023b")
}

g <- load_ghi()
# Arm means, keyed on Challenge_Type as well as Treatment_Group. The earlier
# form keyed on (Study_ID, Treatment_Group) alone and took the FIRST class for
# each group. Where one group label spans arms of different class -- as in
# Massacci_2019, whose group "C" holds a pre-challenge d21 Control row and
# three post-challenge d40 Pathogen_Challenge rows -- that averaged rows of
# different class into one arm and gave it a single label, erasing the control
# arm and dropping the study from both contrasts. Keying on the class keeps the
# arms distinct; the challenge contrast then excludes that study explicitly,
# through NO_CONCURRENT_CONTROL, for the documented reason.
A <- aggregate(list(GHI = g$GHI_Combined_raw),
               by = list(Study_ID = g$Study_ID,
                         Treatment_Group = g$Treatment_Group,
                         Class = g$Intervention_Class),
               FUN = function(z) mean(z, na.rm = TRUE))

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

# Cohort identifiers that are not their own publication. Ma_2024 reports two
# independent cohorts -- ~160 chicks on the FMT trial and ~160 hens on the ALA
# trial, one BioProject (PRJNA1086343) -- and the hen cohort carries the
# biomarker Study_ID "Ma_2024_Layers", which has no row of its own in the Study
# or Risk_of_Bias sheets. Study-level exclusions must be resolved through the
# parent ID or that cohort silently survives every subset below regardless of
# the parent study's platform and risk-of-bias rating. k therefore counts
# independent experiments, not publications: additive k = 36 arises from 35
# publications and challenge k = 25 from 24.
parent_study <- function(x) sub("_Layers$", "", x)

# The archived pipeline tested the mean study effect with a one-sample t-test.
# Both that and the rank-based Wilcoxon are reported; the manuscript quotes the
# t-test p so that the figures match the pipeline outputs already on disk.
both_contrasts <- function(A, label) {
  ad <- study_contrast(A, "GHI", ADDITIVE_CLASS_V2,
                       control_class = CONTROL_CLASS_V2)
  ch <- study_contrast(A[!A$Study_ID %in% NO_CONCURRENT_CONTROL, ],
                       "GHI", CHALLENGE_CLASSES_V2, invert = TRUE,
                       control_class = CONTROL_CLASS_V2)
  pt <- function(e) t.test(e)$p.value
  pw <- function(e) suppressWarnings(wilcox.test(e))$p.value
  cat(sprintf("%-30s additive %+.3f (k=%2d, t p=%.3g; W p=%.3g)   challenge %+.3f (k=%2d, t p=%.3g; W p=%.3g)\n",
              label, mean(ad$effect), nrow(ad), pt(ad$effect), pw(ad$effect),
              mean(ch$effect), nrow(ch), pt(ch$effect), pw(ch$effect)))
  data.frame(analysis = label,
             additive = mean(ad$effect), k_add = nrow(ad),
             add_p_t = pt(ad$effect), add_p_wilcox = pw(ad$effect),
             challenge = mean(ch$effect), k_ch = nrow(ch),
             ch_p_t = pt(ch$effect), ch_p_wilcox = pw(ch$effect))
}

cat("Unit of analysis: k counts independent experiments, not publications.\n",
    "Ma_2024 contributes two cohorts (chicks/FMT as Ma_2024, hens/ALA as\n",
    "Ma_2024_Layers) from one paper and one BioProject, so additive k = 36\n",
    "arises from 35 publications and challenge k = 25 from 24. The challenge\n",
    "contrast excludes ", paste(NO_CONCURRENT_CONTROL, collapse = ", "),
    " for want of a contemporaneous control arm.\n\n", sep = "")

hdr("Registered sensitivity analyses")
res <- list()
res[[1]] <- both_contrasts(A, "All studies (primary)")
res[[2]] <- both_contrasts(A[!parent_study(A$Study_ID) %in% HIGH_STUDIES, ], "Excluding high risk of bias")
# IMPORTANT. The Illumina restriction is defined from the STUDY-level record
# (workbook sheet "Study", column "Seq"), not from the per-run Platform field in
# the sequencing metadata. The two disagree: run-level metadata reports ILLUMINA
# for runs belonging to studies whose own platform was never confirmed, and
# excludes two studies that the study-level record accepts. Using the run-level
# field silently substitutes a different subset of the same size (k = 33 / 23)
# and changes the answer to +0.849 / +0.826. Three studies are excluded here:
# Ma_2026 (PacBio), Tian_2026 (not provided), Wang_2024 (platform unspecified).
if (have_xl) {
  S <- as.data.frame(readxl::read_excel(WB, sheet = "Study"))
  non_ill <- S$Study_ID[!grepl("ILLUMINA", toupper(S$Seq))]
} else {
  non_ill <- c("Ma_2026", "Tian_2026", "Wang_2024")
}
cat("excluded from the Illumina-only subset:", paste(non_ill, collapse = ", "), "\n")
res[[3]] <- both_contrasts(A[!parent_study(A$Study_ID) %in% non_ill, ], "Illumina platforms only")
res[[4]] <- both_contrasts(A[A$Study_ID != "Lucke_2018", ], "Excluding Lucke_2018")
out <- do.call(rbind, res)
write.csv(out, file.path(DIR_OUT, "rob_sensitivity_analyses.csv"), row.names = FALSE)

hdr("Publication bias -- precision proxy")
cat("The registered Egger test needs per-arm standard errors, which published\n",
    "reports do not provide. Substitute: regress the study effect on 1/sqrt(n_arms),\n",
    "a precision proxy, and test rank correlation (Begg).\n\n", sep = "")
for (nm in c("additive", "challenge")) {
  ct <- if (nm == "additive") study_contrast(A, "GHI", ADDITIVE_CLASS_V2,
                                            control_class = CONTROL_CLASS_V2)
        else study_contrast(A[!A$Study_ID %in% NO_CONCURRENT_CONTROL, ],
                            "GHI", CHALLENGE_CLASSES_V2, invert = TRUE,
                            control_class = CONTROL_CLASS_V2)
  prec <- 1 / sqrt(ct$n_arms)   # contrast-specific, not the study's total arm count
  # Deposit the exact contrast object this test is computed on. Script 61
  # writes contrast_additive.csv and contrast_challenge.csv from a different
  # construction with a different n_arms, so re-running this spec on those
  # files does not reproduce the slopes below. These two files do.
  write.csv(data.frame(study = ct$study, effect = ct$effect,
                       n_arms = ct$n_arms, precision_proxy = prec),
            file.path(DIR_OUT, sprintf("funnel_contrast_%s.csv", nm)),
            row.names = FALSE)
  fit <- lm(ct$effect ~ prec)
  tau <- suppressWarnings(cor.test(ct$effect, prec, method = "kendall"))
  cat(sprintf("%-10s k=%2d  Egger-style slope %+.3f (p=%.3f)   Begg tau %+.3f (p=%.3f)\n",
              nm, nrow(ct), coef(fit)[2], summary(fit)$coefficients[2, 4],
              tau$estimate, tau$p.value))
  if (nm == "challenge") {
    k <- ct$study != "Lucke_2018"
    f2 <- lm(ct$effect[k] ~ prec[k])
    cat(sprintf("%-10s   without Lucke_2018: slope %+.3f (p=%.3f)  <- the asymmetry was one study\n",
                "", coef(f2)[2], summary(f2)$coefficients[2, 4]))
  }
}

hdr("Not executable, and why")
cat("The registered n >= 6 birds per group sensitivity analysis cannot be run:\n",
    "the extraction recorded 'Sample size' as free text such as\n",
    "'32 pooled samples (from 64 birds)'. Declare this rather than omit it.\n", sep = "")
close_log(con)
