#!/usr/bin/env Rscript
# =============================================================================
# 42_displacement_coupling_metaregression.R
#
# Reframes the panel search's null result as a moderator question: not
# "which genera predict health" but "under what conditions does taxonomic
# displacement from control track GHI displacement from control." Built
# following the reframing recommended independently by external literature
# review (Perplexity, September 2026) and cross-checked against this
# project's own 50-paper corpus: three independent sources converged on
# this as the more precedented, more answerable next question, rather than
# a tenth panel-search method.
#
# RUN:  Rscript 42_displacement_coupling_metaregression.R
# TIME: under a minute. Reads existing outputs only; writes nothing back
# into the registered pipeline's files.
#
# -----------------------------------------------------------------------------
# THE DESIGN
# -----------------------------------------------------------------------------
# For every arm in every study with a control arm, two displacements are
# computed relative to that study's own control centroid:
#
#   Taxonomic displacement:  Aitchison (Euclidean-on-CLR) distance from the
#                            arm's genus profile to its study's control mean
#   GHI displacement:        the arm's GHI (Combined and per-domain) minus
#                            its study's control mean
#
# The primary model asks whether taxonomic displacement predicts GHI
# displacement, with study as a random intercept:
#
#   GHI_displacement ~ taxonomic_displacement + (1 | Study_ID)
#
# Four candidate moderators are then tested, one at a time via an
# interaction term, matching the four primary hypotheses this reframing
# was built to test:
#
#   1. Challenge_Type  - is coupling stronger for challenge arms than additive?
#   2. Segment         - does cecal/ileal tissue couple better than others?
#   3. GHI_n_domains    - does coupling strengthen as GHI measurement
#                         precision increases (fewer missing domains)?
#   4. Bird_Type        - does coupling differ between broilers and layers?
#
# Age_d and Region_16S were considered and dropped: Age_d is missing for a
# large fraction of arms once collapsed (biomarker-side field, incomplete
# join coverage) and Region_16S has almost no variation left after the
# V3-V4 restriction (only Wang_2024 and five hard-excluded studies differ,
# and the five are excluded from this corpus entirely) - neither would give
# a moderator test any real power.
#
# -----------------------------------------------------------------------------
# WHAT THIS DOES NOT CLAIM
# -----------------------------------------------------------------------------
# This is not a fitted predictor and produces no panel. It cannot be used
# to score a new sample. It answers a structural question about this
# corpus - where does the weak, real signal documented in Chapter 35
# concentrate - not a diagnostic question about any individual bird.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

CONTROL <- "Control"
set.seed(42)

log_msg("[1/7] Loading genus CLR table and collapsing to arms")
clr <- as.data.frame(fread(file.path(R_OUT_DIR, "genus", "genus_clr_batch_adjusted.tsv")))
gname_col <- names(clr)[1]
X <- as.matrix(clr[, -1, drop = FALSE]); rownames(X) <- clr[[gname_col]]
runs <- colnames(X)

md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2study <- setNames(trimws(as.character(md$Study_ID)), rid)
r2arm   <- setNames(trimws(as.character(md$Treatment_Group)), rid)
r2seg   <- setNames(trimws(as.character(md$Segment)), rid)
r2bird  <- setNames(trimws(as.character(md$Bird_Type)), rid)

common <- intersect(runs, names(r2study))
arm_key <- paste(r2study[common], r2arm[common], sep = "||")
uniq_arms <- unique(arm_key)
A <- t(vapply(uniq_arms, function(a) rowMeans(X[, common[arm_key == a], drop = FALSE], na.rm = TRUE),
              numeric(nrow(X))))
colnames(A) <- rownames(X)
arm_study <- sub("\\|\\|.*$", "", uniq_arms)
log_msg("  %d arms across %d studies (genus space)", nrow(A), length(unique(arm_study)))

log_msg("[2/7] Loading GHI scores (Combined + 7 domains + coverage)")
ghi <- as.data.frame(fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")))
dom_cols <- grep("^GHI_dom_", names(ghi), value = TRUE)
ghi_arm <- aggregate(ghi[, c("GHI_Combined_raw", dom_cols, "GHI_n_domains")],
                      by = list(Study_ID = trimws(ghi$Study_ID), Treatment_Group = trimws(ghi$Treatment_Group)),
                      FUN = function(x) mean(x, na.rm = TRUE))
ghi_arm$key <- paste(ghi_arm$Study_ID, ghi_arm$Treatment_Group, sep = "||")
ghi_map <- setNames(seq_len(nrow(ghi_arm)), ghi_arm$key)

# challenge type, segment, bird type per arm - take the first non-NA per key
cls_map <- setNames(rep(NA_character_, length(uniq_arms)), uniq_arms)
seg_map <- setNames(rep(NA_character_, length(uniq_arms)), uniq_arms)
bird_map <- setNames(rep(NA_character_, length(uniq_arms)), uniq_arms)
ghi_cls <- setNames(trimws(as.character(ghi$Challenge_Type)),
                     paste(trimws(ghi$Study_ID), trimws(ghi$Treatment_Group), sep = "||"))
for (a in uniq_arms) {
  if (a %in% names(ghi_cls)) cls_map[a] <- ghi_cls[a]
  rr <- common[arm_key == a]
  if (length(rr)) {
    seg_map[a]  <- names(sort(table(r2seg[rr]), decreasing = TRUE))[1]
    bird_map[a] <- names(sort(table(r2bird[rr]), decreasing = TRUE))[1]
  }
}

log_msg("[3/7] Computing per-arm displacement from each study's own control centroid")
studies_with_control <- unique(arm_study[cls_map[uniq_arms] == CONTROL & !is.na(cls_map[uniq_arms])])
log_msg("  %d studies have a Control arm", length(studies_with_control))

keep_idx <- arm_study %in% studies_with_control
uniq_arms2 <- uniq_arms[keep_idx]; arm_study2 <- arm_study[keep_idx]

# Two displacement measures, deliberately: unsigned Aitchison distance (does
# raw magnitude of change matter, regardless of direction) and the signed
# control-to-challenge axis (does displacement IN THE DIRECTION already
# shown informative in Chapter 35 track GHI). Chapter 35 proved direction
# carries the real signal; reporting only the unsigned distance and calling
# it null would be testing the weaker of the two measures and mislabelling
# the result. Both are computed and reported.
CHALLENGE <- c("Pathogen_Challenge", "Antibiotic_Challenge", "Dietary_Feed_Challenge")
donors <- unique(arm_study[cls_map[uniq_arms] %in% CHALLENGE])
donors <- intersect(donors, studies_with_control)
log_msg("  %d studies can donate a signed control-to-challenge direction", length(donors))

pooled_dir <- function(exclude_study) {
  ds <- setdiff(donors, exclude_study)
  vecs <- list()
  for (s in ds) {
    ctrl <- uniq_arms[arm_study == s & cls_map[uniq_arms] == CONTROL & !is.na(cls_map[uniq_arms])]
    chal <- uniq_arms[arm_study == s & cls_map[uniq_arms] %in% CHALLENGE & !is.na(cls_map[uniq_arms])]
    if (!length(ctrl) || !length(chal)) next
    d <- colMeans(A[chal, , drop = FALSE]) - colMeans(A[ctrl, , drop = FALSE])
    n <- sqrt(sum(d^2)); if (n > 0) vecs[[length(vecs) + 1]] <- d / n
  }
  if (!length(vecs)) return(NULL)
  v <- Reduce(`+`, vecs) / length(vecs)
  v / sqrt(sum(v^2))
}

tax_disp <- numeric(length(uniq_arms2))
signed_axis <- rep(NA_real_, length(uniq_arms2))
ghi_disp <- numeric(length(uniq_arms2))
dom_disp <- matrix(NA_real_, nrow = length(uniq_arms2), ncol = length(dom_cols),
                    dimnames = list(NULL, dom_cols))

for (i in seq_along(uniq_arms2)) {
  s <- arm_study2[i]
  ctrl_arms <- uniq_arms[arm_study == s & cls_map[uniq_arms] == CONTROL & !is.na(cls_map[uniq_arms])]
  if (!length(ctrl_arms)) { tax_disp[i] <- NA; ghi_disp[i] <- NA; next }
  ctrl_centroid <- colMeans(A[ctrl_arms, , drop = FALSE])
  centred <- A[uniq_arms2[i], ] - ctrl_centroid
  tax_disp[i] <- sqrt(sum(centred^2))   # Aitchison = Euclidean on CLR

  v <- pooled_dir(s)   # leave-one-study-out: this study never donates its own direction
  if (!is.null(v)) signed_axis[i] <- sum(centred * v)

  g_row <- ghi_map[uniq_arms2[i]]
  g_ctrl_rows <- ghi_map[ctrl_arms]
  g_ctrl_rows <- g_ctrl_rows[!is.na(g_ctrl_rows)]
  if (is.na(g_row) || !length(g_ctrl_rows)) { ghi_disp[i] <- NA; next }
  ghi_disp[i] <- ghi_arm$GHI_Combined_raw[g_row] - mean(ghi_arm$GHI_Combined_raw[g_ctrl_rows], na.rm = TRUE)
  for (dc in dom_cols) {
    dom_disp[i, dc] <- ghi_arm[[dc]][g_row] - mean(ghi_arm[[dc]][g_ctrl_rows], na.rm = TRUE)
  }
}

df <- data.frame(
  arm = uniq_arms2, Study_ID = arm_study2,
  tax_disp = tax_disp, signed_axis = signed_axis, ghi_disp = ghi_disp,
  Challenge_Type = cls_map[uniq_arms2], Segment = seg_map[uniq_arms2], Bird_Type = bird_map[uniq_arms2],
  GHI_n_domains = ghi_arm$GHI_n_domains[ghi_map[uniq_arms2]],
  dom_disp, check.names = FALSE, stringsAsFactors = FALSE
)
df <- df[!is.na(df$tax_disp) & !is.na(df$ghi_disp) & df$Challenge_Type != CONTROL, ]
df <- df[!is.na(df$Challenge_Type), ]
log_msg("  %d treated arms with both displacements computed, across %d studies",
        nrow(df), length(unique(df$Study_ID)))
log_msg("  of those, %d have a signed-axis value (excludes studies with no LOSO-eligible direction)",
        sum(!is.na(df$signed_axis)))

log_msg("")
log_msg("[4/7] PRIMARY MODELS: does displacement predict GHI displacement overall?")
log_msg("  --- unsigned Aitchison distance (magnitude only) ---")
m0 <- lmer(ghi_disp ~ tax_disp + (1 | Study_ID), data = df)
s0 <- summary(m0)
log_msg("  coefficient: %.4f, SE %.4f, p = %.4f",
        s0$coefficients["tax_disp", "Estimate"], s0$coefficients["tax_disp", "Std. Error"],
        s0$coefficients["tax_disp", "Pr(>|t|)"])

log_msg("  --- signed control-to-challenge axis (direction, proven informative in Chapter 35) ---")
df_sa <- df[!is.na(df$signed_axis), ]
m0b <- tryCatch(lmer(ghi_disp ~ signed_axis + (1 | Study_ID), data = df_sa), error = function(e) NULL)
if (!is.null(m0b)) {
  s0b <- summary(m0b)
  log_msg("  coefficient: %.4f, SE %.4f, p = %.4f (n=%d)",
          s0b$coefficients["signed_axis", "Estimate"], s0b$coefficients["signed_axis", "Std. Error"],
          s0b$coefficients["signed_axis", "Pr(>|t|)"], nrow(df_sa))
} else {
  log_msg("  FAILED to fit")
}

log_msg("")
log_msg("[5/7] MODERATOR TESTS on signed axis - interaction, one moderator at a time")

test_moderator <- function(modname) {
  d2 <- df_sa[!is.na(df_sa[[modname]]), ]
  if (length(unique(d2[[modname]])) < 2 || nrow(d2) < 20) {
    log_msg("  %-16s SKIPPED - insufficient variation or n (%d rows)", modname, nrow(d2))
    return(invisible(NULL))
  }
  form_int <- as.formula(sprintf("ghi_disp ~ signed_axis * %s + (1 | Study_ID)", modname))
  m <- tryCatch(lmer(form_int, data = d2), error = function(e) NULL)
  if (is.null(m)) { log_msg("  %-16s FAILED to fit", modname); return(invisible(NULL)) }
  s <- summary(m)
  int_rows <- grep(":", rownames(s$coefficients))
  log_msg("  --- %s (n=%d) ---", modname, nrow(d2))
  for (r in int_rows) {
    log_msg("    interaction %-40s coef=%+.4f  p=%.4f",
            rownames(s$coefficients)[r], s$coefficients[r, "Estimate"], s$coefficients[r, "Pr(>|t|)"])
  }
}

test_moderator("Challenge_Type")
test_moderator("Segment")
test_moderator("Bird_Type")

log_msg("  --- GHI_n_domains (continuous) ---")
m_cov <- tryCatch(lmer(ghi_disp ~ signed_axis * GHI_n_domains + (1 | Study_ID), data = df_sa),
                   error = function(e) NULL)
if (!is.null(m_cov)) {
  s <- summary(m_cov)
  r <- grep(":", rownames(s$coefficients))
  log_msg("    interaction signed_axis:GHI_n_domains  coef=%+.4f  p=%.4f",
          s$coefficients[r, "Estimate"], s$coefficients[r, "Pr(>|t|)"])
}

log_msg("")
log_msg("[6/7] DOMAIN-SPECIFIC COUPLING (signed axis): which domains track displacement best?")
for (dc in dom_cols) {
  d3 <- df_sa[!is.na(df_sa[[dc]]), ]
  if (nrow(d3) < 15) { log_msg("  %-32s SKIPPED (n=%d)", dc, nrow(d3)); next }
  form <- as.formula(sprintf("`%s` ~ signed_axis + (1 | Study_ID)", dc))
  m <- tryCatch(lmer(form, data = d3), error = function(e) NULL)
  if (is.null(m)) { log_msg("  %-32s FAILED to fit", dc); next }
  s <- summary(m)
  log_msg("  %-32s n=%-4d coef=%+.4f  p=%.4f", dc, nrow(d3),
          s$coefficients["signed_axis", "Estimate"], s$coefficients["signed_axis", "Pr(>|t|)"])
}

log_msg("")
log_msg("[7/7] LEAVE-ONE-STUDY-OUT: does the signed-axis primary model generalise?")
studies_u <- unique(df_sa$Study_ID)
pred <- rep(NA_real_, nrow(df_sa))
for (s in studies_u) {
  train <- df_sa[df_sa$Study_ID != s, ]; test_idx <- which(df_sa$Study_ID == s)
  if (length(unique(train$Study_ID)) < 3 || length(test_idx) == 0) next
  m <- tryCatch(lmer(ghi_disp ~ signed_axis + (1 | Study_ID), data = train), error = function(e) NULL)
  if (is.null(m)) next
  fe <- fixef(m)
  pred[test_idx] <- fe["(Intercept)"] + fe["signed_axis"] * df_sa$signed_axis[test_idx]
}
ok <- !is.na(pred)
rho <- suppressWarnings(cor(df_sa$ghi_disp[ok], pred[ok], method = "spearman"))
r2  <- max(0, 1 - sum((df_sa$ghi_disp[ok] - pred[ok])^2) / sum((df_sa$ghi_disp[ok] - mean(df_sa$ghi_disp[ok]))^2))
log_msg("  LOSO rho = %.3f, R2 = %.3f, n = %d arms, %d studies", rho, r2, sum(ok), length(studies_u))

out <- df
fwrite(out, file.path(R_OUT_DIR, "gdi", "displacement_coupling_data.tsv"), sep = "\t")
log_msg("")
log_msg("=== DONE ===")
log_msg("Arm-level displacement data -> %s", file.path(R_OUT_DIR, "gdi", "displacement_coupling_data.tsv"))
log_msg("")
log_msg("READING: this tests WHEN taxonomic displacement tracks GHI displacement, not WHICH")
log_msg("genera predict health. A significant interaction means the coupling strength differs")
log_msg("by that moderator - e.g. stronger for challenge arms than additive, or in one gut")
log_msg("segment over another. A null primary-model result (Step 4) with no significant")
log_msg("moderators would mean the weak Chapter 35 signal is diffuse across conditions rather")
log_msg("than concentrated in an identifiable subset - itself an informative, reportable finding.")
