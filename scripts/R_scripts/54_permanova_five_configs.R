# 54_permanova_five_configs.R -------------------------------------------------
# The anchor result: community composition versus the Gut Health Index, under
# five feature/sample configurations, with bird age as a positive control.
#
# EXPECTED:
#   A  84 dense, published        141 arms / 41 studies  GHI R2 0.0025  p ~ 0.36
#   B  115 restored, same samples 141 / 41               GHI R2 ~0.002
#   C  115 restored, all samples  141 / 41               GHI R2 ~0.002
#   D  A, organ-matched            58 / 19               GHI R2 ~0.005
#   E  C, organ-matched            58 / 19               GHI R2 ~0.005
#   Age R2 0.08-0.14 in every configuration, p ~ 0.001 -- the positive control.
#   (Age is tested with FREE permutation; the GHI with blocked -- see below.)
# GHI p-values move with the seed; report to two decimals, never three.
# If any GHI p prints as 1.0000 the permutations are not blocked -- see below.
# -----------------------------------------------------------------------------
source("50_setup.R")
suppressMessages(library(vegan))
suppressMessages(library(zCompositions))
con <- open_log("54_permanova_five_configs")

meta <- load_seq_meta()
ghiA <- ghi_by_arm()

clr_from_counts <- function(cnt, z.delete) {
  p <- cmultRepl(t(cnt), label = 0, method = "CZM", z.delete = z.delete, suppress.print = TRUE)
  t(apply(log(p), 1, function(r) r - mean(r)))          # samples x features
}

build_config <- function(label, feature_mat, restrict_matched) {
  mm <- meta
  if (restrict_matched) mm <- mm[mm$Segment_matched == "YES", ]
  keep <- intersect(rownames(feature_mat), mm$sample_id)
  X <- feature_mat[keep, , drop = FALSE]
  mk <- mm[match(keep, mm$sample_id), ]
  arm <- paste(mk$Study_ID, mk$Treatment_Group, sep = "||")
  Xa <- t(sapply(split(seq_len(nrow(X)), arm), function(ii) colMeans(X[ii, , drop = FALSE])))
  ids <- do.call(rbind, strsplit(rownames(Xa), "\\|\\|", fixed = FALSE))
  d <- data.frame(Study_ID = ids[, 1], Treatment_Group = ids[, 2], stringsAsFactors = FALSE)
  d <- merge(d, ghiA, by = c("Study_ID", "Treatment_Group"), sort = FALSE)
  Xa <- Xa[paste(d$Study_ID, d$Treatment_Group, sep = "||"), , drop = FALSE]
  ok <- is.finite(d$GHI)
  d <- d[ok, , drop = FALSE]; Xa <- Xa[ok, , drop = FALSE]
  # Drop single-arm studies. A study contributing one arm carries no
  # within-study contrast and, once blocked permutation is used, contributes a
  # block that cannot be permuted. This is what takes 142 arms / 42 studies to
  # the 141 / 41 of the published analysis.
  tb <- table(d$Study_ID); multi <- names(tb)[tb >= 2]
  k <- d$Study_ID %in% multi
  list(label = label, X = Xa[k, , drop = FALSE], d = d[k, , drop = FALSE])
}

# Restrict a configuration to a reference set of arms, so that configurations
# differing only in FEATURES are compared on identical samples.
restrict_to <- function(cfg, keys) {
  kk <- paste(cfg$d$Study_ID, cfg$d$Treatment_Group, sep = "||") %in% keys
  cfg$X <- cfg$X[kk, , drop = FALSE]; cfg$d <- cfg$d[kk, , drop = FALSE]
  cfg
}

# PERMUTATION SCHEME. Study identity explains ~43% of community variation even
# after batch correction, so a free permutation is not a valid null: it asks
# whether the GHI beats a variable with no study structure at all, which it
# always loses, and the p-value pins at 1.0000. Permutations are therefore
# restricted WITHIN study (how(blocks = Study_ID)) and the predictor is entered
# raw -- the blocking, not a centring of the predictor, is what removes the
# between-study signal. This reproduces the published R2 = 0.0025, p ~ 0.36.
run_permanova <- function(cfg, nperm = N_PERM) {
  d <- cfg$d; X <- cfg$X
  ageok <- is.finite(d$Age_d)
  D <- vegdist(X, method = "euclidean")           # Aitchison distance on CLR
  ctrl <- how(nperm = nperm, blocks = factor(d$Study_ID))
  set.seed(SEED)
  a1 <- adonis2(D ~ GHI, data = d, permutations = ctrl, by = "terms")
  r_age <- NA; p_age <- NA
  # AGE USES FREE PERMUTATION, DELIBERATELY. Bird age at sampling is constant
  # within study in 40 of the 41 studies, so permuting it within study is a
  # no-op and forces p = 1 by construction. Age is a BETWEEN-study variable and
  # its null is the free permutation; the GHI varies WITHIN study and its null
  # is the blocked one. Using one scheme for both is wrong in one direction or
  # the other. State this in the methods -- a reviewer will ask why the two
  # terms in the same table use different permutation schemes.
  if (sum(ageok) > 10) {
    set.seed(SEED)
    a2 <- adonis2(as.dist(as.matrix(D)[ageok, ageok]) ~ Age_d, data = d[ageok, ],
                  permutations = nperm, by = "terms")
    r_age <- a2$R2[1]; p_age <- a2$`Pr(>F)`[1]
  }
  data.frame(config = cfg$label, arms = nrow(d), studies = length(unique(d$Study_ID)),
             GHI_R2 = round(a1$R2[1], 5), GHI_p = a1$`Pr(>F)`[1],
             Age_R2 = round(r_age, 5), Age_p = p_age, stringsAsFactors = FALSE)
}

hdr("Building feature matrices")
cnt <- load_counts()
clr_dense <- t(load_clr_published())                       # 1787 x 84
cat("published CLR:", nrow(clr_dense), "samples x", ncol(clr_dense), "features\n")
clr_full  <- clr_from_counts(cnt, z.delete = FALSE)        # 2064 x 115
cat("restored CLR :", nrow(clr_full), "samples x", ncol(clr_full), "features\n")
clr_fullB <- clr_full[rownames(clr_full) %in% rownames(clr_dense), , drop = FALSE]

cfgA <- build_config("A 84 dense (published)",       clr_dense, FALSE)
cfgD <- build_config("D A, organ-matched",           clr_dense, TRUE)
keyA <- paste(cfgA$d$Study_ID, cfgA$d$Treatment_Group, sep = "||")
keyD <- paste(cfgD$d$Study_ID, cfgD$d$Treatment_Group, sep = "||")
configs <- list(
  cfgA,
  restrict_to(build_config("B 115 restored, same samples", clr_fullB, FALSE), keyA),
  restrict_to(build_config("C 115 restored, all samples",  clr_full,  FALSE), keyA),
  cfgD,
  restrict_to(build_config("E C, organ-matched",           clr_full,  TRUE),  keyD)
)

hdr(paste0("PERMANOVA, ", N_PERM, " permutations"))
res <- do.call(rbind, lapply(configs, run_permanova))
print(res, row.names = FALSE)
write.csv(res, file.path(DIR_OUT, "permanova_five_configs.csv"), row.names = FALSE)

cat("\nThe age column is the positive control: the same test on the same data\n",
    "detects a large, significant effect for a trivial covariate.\n", sep = "")
close_log(con)
