#!/usr/bin/env Rscript
# =============================================================================
# 32_bactotraits.R
#
# Trait-based analysis: describe the community by what its bacteria ARE
# (measured phenotypes) rather than by which genera are present.
#
# RUN:  Rscript 32_bactotraits.R
# TIME: 10-20 minutes.
#
# -----------------------------------------------------------------------------
# WHY THIS IS NOT A REPEAT OF THE PICRUSt2 ATTEMPT
# -----------------------------------------------------------------------------
# PICRUSt2 PREDICTS gene content by placing each ASV on a reference phylogeny.
# Nothing is measured, and 31_picrust_scfa_check.R showed those predictions did
# not track measured SCFAs in this corpus.
#
# BactoTraits is different in kind. It records phenotypes MEASURED in culture -
# oxygen tolerance, Gram stain, spore formation, cell shape - taken from
# BacDive. Nothing is inferred from a tree. The PICRUSt2 result therefore says
# nothing about whether this will work.
#
# It also targets the specific failure diagnosed earlier in this project.
# Run-level classifiers were memorising study identity, and part of why that is
# possible is that different primers and regions recover different GENERA. A
# trait like "anaerobic spore-former" is more stable across those technical
# differences than a genus name is: two studies may detect different organisms
# filling the same role, and at trait level those converge.
#
# -----------------------------------------------------------------------------
# THE COVERAGE PROBLEM, AND WHY IT IS REPORTED BY ABUNDANCE
# -----------------------------------------------------------------------------
# 45 of 84 genera match BactoTraits. The 39 that do not are almost entirely
# SILVA placeholders for uncultured organisms - CHKCI002, CAG-352, UCG-005,
# Clostridia_vadinBB60_group and similar. They have no traits precisely because
# nobody has cultured them.
#
# Counting genera understates or overstates the problem depending on how
# abundant the missing ones are. What matters is what FRACTION OF THE COMMUNITY
# is covered, so that is computed and reported per sample before anything else.
# If the uncultured genera are abundant, the trait profile describes a minority
# of the bacteria actually present and every downstream result inherits that.
#
# Several unmatched genera were among the strongest hits in the genus-level
# work (CHKCI002 had q = 1.9e-10 in MaAsLin2; CAG-352 led the challenge
# contrast). This analysis necessarily discards them.
#
# -----------------------------------------------------------------------------
# WHAT IS COMPUTED
# -----------------------------------------------------------------------------
#   CWM  - community weighted mean. For each trait modality, the
#          abundance-weighted proportion of the community carrying it.
#          "What is this community like on average."
#
#   Rao  - Rao's quadratic entropy. The abundance-weighted mean trait distance
#          between pairs of taxa. "How varied is this community."
#          Two communities can share an identical CWM and differ completely in
#          Rao; diversity often responds to stress before composition does,
#          which is why the source paper reports both.
#
# Traits irrelevant to a chicken gut are excluded: halophily (salt tolerance),
# colony colour, pigmentation, phototrophy. Retained are oxygen tolerance,
# Gram stain, spore formation, morphology, motility, pH and temperature optima,
# GC content, genome size, 16S copy number, antibiotic resistance and biosafety
# level.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' required.")
suppressPackageStartupMessages(library(glmnet))

BT <- Sys.getenv("BACTOTRAITS_CSV",
      unset = file.path(PDI_ROOT, "reference", "BACTOTRAITS_DATASET_2026-01-28_GENUSLVL.csv"))
if (!file.exists(BT)) stop("Not found: ", BT)
N_BOOT <- 1000; SUBFRAC <- 0.50; STABILITY_THR <- 0.70; N_PERM <- 5000
set.seed(42)
HEALTHY <- "Control"; CHALLENGE <- c("Pathogen_Challenge","Antibiotic_Challenge","Dietary_Feed_Challenge")

log_msg("[1/7] Loading BactoTraits (genus level, 2026 release)")
# decimals are comma-separated in this file; read as character then convert
bt <- as.data.frame(fread(BT, dec = ",", check.names = FALSE), stringsAsFactors = FALSE)
log_msg("  %d genera x %d columns", nrow(bt), ncol(bt))

TAXCOLS <- c("Total_BacDive_ids","Phylum","Class","Order","Family","Genus")
trait_all <- setdiff(names(bt), TAXCOLS)
# drop trait groups with no plausible meaning in a chicken caecum
DROP <- "halophily|colony_color|pigmentation|phototroph|biosafety"
trait_cols <- trait_all[!grepl(DROP, trait_all, ignore.case = TRUE)]
log_msg("  %d trait modalities retained (%d dropped as irrelevant to a gut)",
        length(trait_cols), length(trait_all) - length(trait_cols))

for (cn in trait_cols) {
  v <- bt[[cn]]
  if (is.character(v)) v <- suppressWarnings(as.numeric(gsub(",", ".", v)))
  bt[[cn]] <- v
}
bt$Genus <- trimws(as.character(bt$Genus))
bt$n_strains <- suppressWarnings(as.numeric(bt$Total_BacDive_ids))

log_msg("[2/7] Matching to this project's genera")
prop <- fread(file.path(R_OUT_DIR, "genus", "genus_proportions_batch_adjusted.tsv"))
gname <- prop[[1]]
P <- as.matrix(prop[, -1, with = FALSE]); rownames(P) <- gname
gclean <- gsub("^g__", "", gname)
# Escherichia-Shigella is a SILVA compound label; BactoTraits holds Escherichia
gmatch <- gclean
gmatch[gmatch == "Escherichia-Shigella"] <- "Escherichia"

idx <- match(gmatch, bt$Genus)
matched <- !is.na(idx)
log_msg("  %d of %d genera matched (%.0f%%)", sum(matched), length(gclean), 100*mean(matched))

# ---- coverage BY ABUNDANCE, the figure that actually matters ---------------
colsum_all <- colSums(P, na.rm = TRUE)
colsum_cov <- colSums(P[matched, , drop = FALSE], na.rm = TRUE)
cov_frac <- ifelse(colsum_all > 0, colsum_cov/colsum_all, NA_real_)
log_msg("")
log_msg("  ABUNDANCE-WEIGHTED COVERAGE per sample:")
log_msg("    median %.1f%%, mean %.1f%%, range %.1f%% - %.1f%%",
        100*median(cov_frac, na.rm=TRUE), 100*mean(cov_frac, na.rm=TRUE),
        100*min(cov_frac, na.rm=TRUE), 100*max(cov_frac, na.rm=TRUE))
log_msg("    samples below 50%% covered: %d of %d",
        sum(cov_frac < 0.5, na.rm=TRUE), length(cov_frac))
if (median(cov_frac, na.rm=TRUE) < 0.5) {
  log_msg("    *** The traited genera are a MINORITY of the community. Every")
  log_msg("    *** result below describes less than half the bacteria present.")
} else {
  log_msg("    The traited genera account for the majority of the community.")
}

# strain support: a genus backed by one cultured strain gives a shakier profile
ns <- bt$n_strains[idx[matched]]
log_msg("")
log_msg("  strains per matched genus: median %.0f, %d genera backed by a single strain",
        median(ns, na.rm=TRUE), sum(ns == 1, na.rm=TRUE))

T <- as.matrix(bt[idx[matched], trait_cols, drop = FALSE])
rownames(T) <- gname[matched]
Pm <- P[matched, , drop = FALSE]
keep_t <- apply(T, 2, function(v) sum(is.finite(v)) > 0 && stats::sd(v, na.rm=TRUE) > 0)
T <- T[, keep_t, drop = FALSE]
T[!is.finite(T)] <- 0
log_msg("  %d traits x %d genera usable", ncol(T), nrow(T))

log_msg("[3/7] Computing community weighted means (CWM)")
w <- sweep(Pm, 2, pmax(colSums(Pm, na.rm=TRUE), 1e-12), "/")   # renormalise to covered taxa
CWM <- t(t(T) %*% w)     # samples x traits
log_msg("  %d samples x %d trait CWMs", nrow(CWM), ncol(CWM))

# ---------------------------------------------------------------------------
# CWMs ARE COMPOSITIONAL AND THE MATRIX IS RANK DEFICIENT
#
# BactoTraits is fuzzy-coded: within each trait, the modalities sum to 1 for
# every genus, and that property survives abundance weighting. The CWM matrix
# therefore contains exact linear dependencies - 79 pairs at |r| = 1.000, with
# rank 61 of 73 columns in this corpus.
#
# An earlier attempt tried to fix this by parsing trait names and dropping one
# modality per group. That failed: names like "cell_length_<=_0.9" and
# "cell_length_>_0.9_to_1.3" do not share a stem under simple suffix
# stripping, so only 9 of the ~20 groups were detected and the dependencies
# remained. glmnet returned max stability 0% and R2 = NA - the solver could
# not fit at all, which is not the same as finding no signal.
#
# The reliable fix does not depend on naming. QR decomposition with column
# pivoting identifies exactly which columns are linearly dependent on the
# others, whatever they are called, and those are dropped. What remains is a
# full-rank basis spanning the same space: no information is lost, because a
# dropped column was by definition reconstructible from the ones kept.
#
# Near-constant columns are removed first. glmnet standardises internally by
# dividing by the standard deviation, so a column with SD near zero produces
# enormous standardised values and destabilises the fit.
# ---------------------------------------------------------------------------
sdv <- apply(CWM, 2, stats::sd, na.rm = TRUE)
n_const <- sum(!is.finite(sdv) | sdv < 1e-4)
if (n_const) log_msg("  dropping %d near-constant trait(s) (SD < 1e-4)", n_const)
CWM <- CWM[, is.finite(sdv) & sdv >= 1e-4, drop = FALSE]

qrd <- qr(CWM, tol = 1e-7)
if (qrd$rank < ncol(CWM)) {
  keep_idx <- sort(qrd$pivot[seq_len(qrd$rank)])
  dropped <- setdiff(colnames(CWM), colnames(CWM)[keep_idx])
  log_msg("  matrix rank %d of %d columns - dropping %d linearly dependent trait(s)",
          qrd$rank, ncol(CWM), length(dropped))
  for (d in head(dropped, 8)) log_msg("      %s", d)
  if (length(dropped) > 8) log_msg("      ... and %d more", length(dropped) - 8)
  CWM <- CWM[, keep_idx, drop = FALSE]
} else log_msg("  matrix is already full rank")
log_msg("  %d traits retained for modelling (full rank, no exact dependencies)", ncol(CWM))

log_msg("[4/7] Computing Rao's quadratic entropy (trait diversity)")
Tz <- scale(T)
Tz[!is.finite(Tz)] <- 0
D <- as.matrix(dist(Tz))
D <- D / max(D)
rao <- apply(w, 2, function(p) as.numeric(t(p) %*% D %*% p))
log_msg("  Rao: median %.4f, range %.4f - %.4f", median(rao), min(rao), max(rao))

log_msg("[5/7] Collapsing to arms and attaching the GHI")
md <- fread(METADATA_TSV, sep="\t", header=TRUE, na.strings=c("","NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2s <- setNames(trimws(as.character(md$Study_ID)), rid)
r2t <- setNames(trimws(as.character(md$Treatment_Group)), rid)
runs <- intersect(rownames(CWM), names(r2s))
CWM <- CWM[runs, , drop=FALSE]; rao <- rao[runs]; covr <- cov_frac[runs]
ak <- paste(r2s[runs], r2t[runs], sep="||")
A   <- t(vapply(unique(ak), function(a) colMeans(CWM[ak==a, , drop=FALSE]), numeric(ncol(CWM))))
raoA <- vapply(unique(ak), function(a) mean(rao[ak==a]), numeric(1))
covA <- vapply(unique(ak), function(a) mean(covr[ak==a]), numeric(1))
astudy <- sub("\\|\\|.*$", "", rownames(A))

gh <- as.data.frame(fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")))
gh$Study_ID <- trimws(gh$Study_ID)
gk <- paste(gh$Study_ID, trimws(gh$Treatment_Group), sep="||")
dom <- grep("^GHI_dom_", names(gh), value=TRUE)
gh$ndom <- rowSums(!is.na(gh[, dom, drop=FALSE]))
armGHI <- function(col, mind = 1) {
  s <- gh$ndom >= mind & is.finite(gh[[col]])
  unname(tapply(gh[[col]][s], gk[s], mean)[rownames(A)])
}
y_all  <- armGHI("GHI_Combined_raw", 1)
y_hi   <- armGHI("GHI_Combined_raw", 3)
y_bio  <- armGHI("GHI_Biomarker_raw", 3)
cls <- unname(tapply(trimws(as.character(gh$Challenge_Type)), gk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing=TRUE))[1] })[rownames(A)])
log_msg("  %d arms across %d studies; GHI available for %d",
        nrow(A), length(unique(astudy)), sum(is.finite(y_all)))

log_msg("[6/7] Rao's Q and coverage against the GHI")
wc <- function(x, y, g) {
  ok <- is.finite(x) & is.finite(y); x<-x[ok]; y<-y[ok]; g<-g[ok]
  if (length(x) < 10) return(c(NA,NA,0,0))
  for (s in unique(g)) { i <- g==s
    if (sum(i)>=2) { x[i]<-x[i]-mean(x[i]); y[i]<-y[i]-mean(y[i]) } else { x[i]<-NA; y[i]<-NA } }
  o <- is.finite(x)&is.finite(y)
  if (sum(o)<10) return(c(NA,NA,sum(o),length(unique(g[o]))))
  ct <- suppressWarnings(cor.test(x[o],y[o],method="spearman",exact=FALSE))
  c(unname(ct$estimate), ct$p.value, sum(o), length(unique(g[o])))
}
for (nm in c("GHI all samples","GHI >=3 domains","Biomarker >=3 domains")) {
  y <- switch(nm, "GHI all samples"=y_all, "GHI >=3 domains"=y_hi, y_bio)
  r <- wc(raoA, y, astudy)
  log_msg("  Rao vs %-24s rho=%+.3f  p=%s  (%d arms, %d studies)",
          nm, r[1], ifelse(is.na(r[2]),"n/a",sprintf("%.4f",r[2])), r[3], r[4])
}
rc <- wc(covA, y_all, astudy)
log_msg("  coverage vs GHI (sanity)     rho=%+.3f  p=%s", rc[1],
        ifelse(is.na(rc[2]),"n/a",sprintf("%.4f",rc[2])))
log_msg("    (a strong correlation here would mean the coverage fraction itself")
log_msg("     tracks health, which would confound the trait results)")

log_msg("[7/7] LASSO on trait CWMs, same rules as the genus analysis")
run_lasso <- function(y, label) {
  ok <- is.finite(y)
  Xr <- A[ok,,drop=FALSE]; yy <- y[ok]; ss <- astudy[ok]
  tb <- table(ss); keep <- ss %in% names(tb)[tb>=2]
  Xr <- Xr[keep,,drop=FALSE]; yy <- yy[keep]; ss <- ss[keep]
  if (nrow(Xr) < 30) { log_msg("  %s: too few arms (%d)", label, nrow(Xr)); return(NULL) }
  Xc <- Xr; yc <- yy
  for (s in unique(ss)) { i <- ss==s
    Xc[i,] <- sweep(Xr[i,,drop=FALSE],2,colMeans(Xr[i,,drop=FALSE]),"-"); yc[i] <- yy[i]-mean(yy[i]) }
  Xc <- Xc[, apply(Xc,2,function(v) stats::sd(v)>0), drop=FALSE]
  studies <- unique(ss); ntake <- max(3, floor(length(studies)*SUBFRAC))
  cnt <- setNames(numeric(ncol(Xc)), colnames(Xc)); nok <- 0
  for (b in seq_len(N_BOOT)) {
    tk <- sample(studies, ntake); ii <- which(ss %in% tk)
    if (length(ii)<10 || length(unique(yc[ii]))<5) next
    f <- tryCatch(cv.glmnet(Xc[ii,,drop=FALSE], yc[ii], alpha=1, nfolds=min(5,length(ii)),
                            standardize=TRUE, intercept=FALSE), error=function(e) NULL)
    if (is.null(f)) next
    co <- as.matrix(coef(f, s=f$lambda.1se)); co <- co[rownames(co)!="(Intercept)",,drop=FALSE]
    sel <- rownames(co)[co[,1]!=0]; cnt[sel] <- cnt[sel]+1; nok <- nok+1
  }
  freq <- if (nok>0) cnt/nok else cnt
  r2 <- tryCatch({ f <- cv.glmnet(Xc,yc,alpha=1,standardize=TRUE,intercept=FALSE)
                   p <- as.numeric(predict(f,newx=Xc,s=f$lambda.1se))
                   if (sd(p)==0) 0 else max(0,1-sum((yc-p)^2)/sum((yc-mean(yc))^2)) },
                 error=function(e) NA_real_)
  log_msg("")
  log_msg("  --- %s: %d arms, %d studies, %d traits ---", label, nrow(Xc), length(studies), ncol(Xc))
  log_msg("      max stability %.0f%%, traits >=70%%: %d, R2 = %.3f",
          100*max(freq), sum(freq>=STABILITY_THR), r2)
  top <- sort(freq, decreasing=TRUE)[1:min(8,length(freq))]
  for (i in seq_along(top)) log_msg("        %-44s %3.0f%%", substr(names(top)[i],1,44), 100*top[i])
  list(label=label, n=nrow(Xc), k=length(studies), max=max(freq),
       n70=sum(freq>=STABILITY_THR), r2=r2, top=top)
}
L1 <- run_lasso(y_all, "GHI, all samples")
L2 <- run_lasso(y_hi,  "GHI, >=3 domains")
L3 <- run_lasso(y_bio, "Biomarker, >=3 domains")

log_msg("")
log_msg("=========================================================")
log_msg("  genus-level comparison: max stability 31%%, 0 at 70%%, R2 0.000")
log_msg("=========================================================")

out <- c("BACTOTRAITS - TRAIT-BASED ANALYSIS", strrep("=",74), "",
  sprintf("Genera matched: %d of %d (%.0f%%)", sum(matched), length(gclean), 100*mean(matched)),
  sprintf("Abundance-weighted coverage: median %.1f%%, range %.1f-%.1f%%",
          100*median(cov_frac,na.rm=TRUE), 100*min(cov_frac,na.rm=TRUE), 100*max(cov_frac,na.rm=TRUE)),
  sprintf("Traits used: %d modalities", ncol(T)), "",
  "Unmatched genera are almost entirely SILVA placeholders for uncultured",
  "organisms, several of which were among the strongest hits in the",
  "genus-level analysis. This analysis necessarily discards them.", "")
for (L in Filter(Negate(is.null), list(L1,L2,L3)))
  out <- c(out, sprintf("%-26s %3d arms %3d studies  max stab %3.0f%%  >=70%%: %d  R2 %.3f",
           L$label, L$n, L$k, 100*L$max, L$n70, L$r2))
out <- c(out, "", "Genus-level comparison: max stability 31%, 0 at 70%, R2 0.000")
writeLines(out, file.path(R_OUT_DIR, "gdi", "bactotraits.txt"))
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "bactotraits.txt"))
log_msg("=== DONE ===")
