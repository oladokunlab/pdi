#!/usr/bin/env Rscript
# =============================================================================
# 33_guild_ratios.R
#
# Do ratios between functional guilds change with treatment?
#
# RUN:  Rscript 33_guild_ratios.R
# TIME: 3-6 minutes.
#
# -----------------------------------------------------------------------------
# WHY A RATIO MIGHT SUCCEED WHERE THE MODELS DID NOT
# -----------------------------------------------------------------------------
# Every earlier attempt handed a model a high-dimensional space - 84 genera,
# 8,892 KEGG orthologs, 79 traits - and asked it to find structure. Each had a
# penalty parameter and a feature-selection step, and the run-level work showed
# what those latch onto: study identity, not dysbiosis.
#
# A ratio between two pre-specified guilds has nothing to fit. No feature
# selection, no penalty, no tuning. It is a single number defined before the
# data is examined, which removes the mechanism that produced the memorisation.
#
# Ratios are also naturally robust to compositional artefacts: a batch effect
# that inflates a whole sample cancels between numerator and denominator.
#
# -----------------------------------------------------------------------------
# TWO INDEPENDENT WAYS OF DEFINING THE GUILDS
# -----------------------------------------------------------------------------
# LITERATURE - genera assigned from published poultry and mammalian gut work.
#   Defensible and interpretable, but it is a judgment call about which genus
#   belongs where, and genera are not functionally uniform (some Clostridium
#   species produce butyrate, others do not).
#
# BACTOTRAITS - guilds defined from measured phenotypes: anaerobe vs
#   facultative, spore-former vs not, Gram-positive vs negative. Objective and
#   free of the author's judgment, but covers only ~48% of the community by
#   abundance because uncultured genera have no trait record.
#
# Both are computed. Agreement between two independently constructed
# definitions would be far stronger evidence than either alone; disagreement
# would indicate the result depends on how the guilds were drawn, which is
# itself worth knowing.
#
# -----------------------------------------------------------------------------
# HOW THE RATIO IS TESTED
# -----------------------------------------------------------------------------
# Log-ratios, since a raw ratio is bounded below by zero and badly skewed.
# Tested by paired within-study differences pooled across studies, the same
# design used for the GHI validation - so a ratio never has to be comparable
# in absolute terms between studies, only in direction within one.
#
# Every ratio is also correlated against the GHI, and a leave-one-study-out
# check reports whether the direction holds when each study is excluded in
# turn. A ratio that only works because of one study is not a finding.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
if (!requireNamespace("metafor", quietly = TRUE)) stop("Package 'metafor' required.")

HEALTHY   <- "Control"
ADDITIVE  <- "Non_Disruptive_Additive"
CHALLENGE <- c("Pathogen_Challenge","Antibiotic_Challenge","Dietary_Feed_Challenge")
BT <- Sys.getenv("BACTOTRAITS_CSV",
      unset = file.path(PDI_ROOT, "reference", "BACTOTRAITS_DATASET_2026-01-28_GENUSLVL.csv"))
set.seed(42)

# =============================================================================
# LITERATURE-DEFINED GUILDS
# Assigned from published poultry and mammalian gut work, before examining any
# result in this project. Genera absent from the corpus are simply unmatched.
# =============================================================================
GUILDS_LIT <- list(
  SCFA_butyrate = c("Faecalibacterium","Butyricicoccus","Roseburia","Anaerostipes",
                    "Blautia","Coprococcus","Subdoligranulum","Eubacterium",
                    "Butyrivibrio","Anaerobutyricum","Agathobacter","Oscillibacter",
                    "Pseudobutyrivibrio","Shuttleworthia"),
  Lactic_acid   = c("Lactobacillus","Limosilactobacillus","Ligilactobacillus",
                    "Lactiplantibacillus","Lacticaseibacillus","Bifidobacterium",
                    "Enterococcus","Streptococcus","Weissella","Pediococcus"),
  Opportunist   = c("Escherichia-Shigella","Escherichia","Salmonella","Campylobacter",
                    "Helicobacter","Clostridium_sensu_stricto_1","Proteus",
                    "Klebsiella","Staphylococcus","Enterobacter","Citrobacter",
                    "Gallibacterium","Avibacterium"),
  Mucin_degrader= c("Akkermansia","Ruminococcus_torques_group","Bacteroides",
                    "Barnesiella","Alistipes","Parabacteroides","Prevotella")
)

log_msg("[1/6] Loading abundances")
prop <- fread(file.path(R_OUT_DIR, "genus", "genus_proportions_batch_adjusted.tsv"))
gname <- prop[[1]]
P <- as.matrix(prop[, -1, with = FALSE]); rownames(P) <- gname
gclean <- gsub("^g__", "", gname)
log_msg("  %d genera x %d runs", nrow(P), ncol(P))

log_msg("[2/6] Assigning literature guilds")
lit_idx <- lapply(GUILDS_LIT, function(v) which(gclean %in% v))
for (nm in names(lit_idx)) {
  hit <- gclean[lit_idx[[nm]]]
  log_msg("  %-16s %2d genera: %s", nm, length(hit),
          substr(paste(hit, collapse=", "), 1, 62))
}
unassigned <- setdiff(gclean, unlist(GUILDS_LIT))
log_msg("  unassigned: %d of %d genera (%.0f%% of the corpus by count)",
        length(unassigned), length(gclean), 100*length(unassigned)/length(gclean))

log_msg("[3/6] Assigning BactoTraits guilds")
trait_idx <- list()
if (file.exists(BT)) {
  bt <- as.data.frame(fread(BT, dec = ",", check.names = FALSE), stringsAsFactors = FALSE)
  bt$Genus <- trimws(as.character(bt$Genus))
  gm <- gclean; gm[gm == "Escherichia-Shigella"] <- "Escherichia"
  mi <- match(gm, bt$Genus)
  numv <- function(cn) { v <- bt[[cn]][mi]
    if (is.character(v)) v <- suppressWarnings(as.numeric(gsub(",", ".", v))); v }
  anaer <- numv("oxygen_tolerance_anaerobe")
  facul <- numv("oxygen_tolerance_facultative")
  aerob <- numv("oxygen_tolerance_aerobe")
  spore <- numv("spore_formation_yes")
  gramp <- numv("gram_stain_positive")
  # a genus joins a guild when the trait clearly dominates (>0.6 of strains)
  trait_idx <- list(
    Strict_anaerobe = which(!is.na(anaer) & anaer > 0.6),
    Facultative     = which(!is.na(facul) & facul > 0.6),
    Aerobe          = which(!is.na(aerob) & aerob > 0.6),
    Spore_former    = which(!is.na(spore) & spore > 0.6),
    Gram_positive   = which(!is.na(gramp) & gramp > 0.6)
  )
  for (nm in names(trait_idx))
    log_msg("  %-16s %2d genera", nm, length(trait_idx[[nm]]))
  log_msg("  (a genus joins a guild when >60%% of its strains carry the trait)")
} else stop(sprintf(paste("BactoTraits not found at %s.",
  "Three of the seven guild ratios are trait-derived and cannot be computed",
  "without it; degrading silently to literature guilds produced a partial",
  "table indistinguishable from a complete one. Set BACTOTRAITS_CSV or place",
  "the file under reference/."), BT))

# =============================================================================
log_msg("[4/6] Computing guild abundances and log-ratios")
guild_abund <- function(idx_list) {
  t(vapply(idx_list, function(i)
    if (!length(i)) rep(NA_real_, ncol(P)) else colSums(P[i, , drop = FALSE], na.rm = TRUE),
    numeric(ncol(P))))
}
GA_lit <- guild_abund(lit_idx)
GA_tr  <- if (length(trait_idx)) guild_abund(trait_idx) else NULL

PSEUDO <- 1e-6
lr <- function(GA, a, b) {
  if (is.null(GA) || !a %in% rownames(GA) || !b %in% rownames(GA)) return(NULL)
  x <- log((GA[a, ] + PSEUDO) / (GA[b, ] + PSEUDO))
  x[is.finite(x)] <- x[is.finite(x)]; x
}
RATIOS <- list()
RATIOS[["LIT: butyrate / opportunist"]]      <- lr(GA_lit, "SCFA_butyrate", "Opportunist")
RATIOS[["LIT: butyrate / lactic acid"]]      <- lr(GA_lit, "SCFA_butyrate", "Lactic_acid")
RATIOS[["LIT: lactic acid / opportunist"]]   <- lr(GA_lit, "Lactic_acid", "Opportunist")
RATIOS[["LIT: mucin degr / butyrate"]]       <- lr(GA_lit, "Mucin_degrader", "SCFA_butyrate")
if (!is.null(GA_tr)) {
  RATIOS[["TRAIT: anaerobe / facultative"]]  <- lr(GA_tr, "Strict_anaerobe", "Facultative")
  RATIOS[["TRAIT: anaerobe / aerobe"]]       <- lr(GA_tr, "Strict_anaerobe", "Aerobe")
  RATIOS[["TRAIT: spore / non-spore"]]       <- {
    sp <- GA_tr["Spore_former", ]
    tot <- colSums(P, na.rm = TRUE)
    log((sp + PSEUDO) / (pmax(tot - sp, 0) + PSEUDO)) }
}
RATIOS <- Filter(Negate(is.null), RATIOS)
log_msg("  %d ratios computed", length(RATIOS))

log_msg("[5/6] Collapsing to arms")
md <- fread(METADATA_TSV, sep="\t", header=TRUE, na.strings=c("","NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2s <- setNames(trimws(as.character(md$Study_ID)), rid)
r2t <- setNames(trimws(as.character(md$Treatment_Group)), rid)
runs <- intersect(colnames(P), names(r2s))
ak <- paste(r2s[runs], r2t[runs], sep="||")
ua <- unique(ak)
armv <- function(v) { v <- v[runs]; vapply(ua, function(a) mean(v[ak==a], na.rm=TRUE), numeric(1)) }
RA <- lapply(RATIOS, armv)
astudy <- sub("\\|\\|.*$", "", ua)

gh <- as.data.frame(fread(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv")))
gh$Study_ID <- trimws(gh$Study_ID)
gk <- paste(gh$Study_ID, trimws(gh$Treatment_Group), sep="||")
dom <- grep("^GHI_dom_", names(gh), value=TRUE)
gh$ndom <- rowSums(!is.na(gh[, dom, drop=FALSE]))
cls <- unname(tapply(trimws(as.character(gh$Challenge_Type)), gk, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) NA_character_ else names(sort(table(z), decreasing=TRUE))[1] })[ua])
ghi3 <- { s <- gh$ndom >= 3 & is.finite(gh$GHI_Combined_raw)
          unname(tapply(gh$GHI_Combined_raw[s], gk[s], mean)[ua]) }
log_msg("  %d arms across %d studies", length(ua), length(unique(astudy)))

# =============================================================================
log_msg("[6/6] Testing each ratio")

# -----------------------------------------------------------------------------
# POOLING - deliberately identical to 04e_mixed_model_validation.R
#
# An earlier version of this script weighted each study by var(c(a,b))/n. That
# is circular: with one arm per group, var() of two values IS the squared
# difference, so a study showing a LARGE effect was assigned a LARGE variance
# and therefore LESS weight. It biased every pooled estimate toward zero and
# produced the "sampling variance extremely large" warnings from metafor.
#
# This corpus cannot support proper inverse-variance weighting - that needs
# each group's SD and sample size, and published papers report group means
# only. So studies are weighted by the NUMBER OF GROUP MEANS behind each
# difference, which is coarser but not circular, and is labelled as such
# rather than presented as a true inverse-variance analysis.
#
# Heterogeneity is estimated by DerSimonian-Laird. Random rather than fixed
# effects, because studies differ in bird age, breed, challenge organism and
# husbandry; assuming one common true effect across that range would give
# artificially narrow intervals.
# -----------------------------------------------------------------------------
contrast_test <- function(v, treat_set) {
  d <- list()
  for (s in unique(astudy)) {
    i <- astudy == s
    a <- v[i & cls %in% HEALTHY]; b <- v[i & cls %in% treat_set]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (!length(a) || !length(b)) next
    d[[length(d) + 1]] <- list(diff = mean(b) - mean(a), n = length(a) + length(b))
  }
  if (length(d) < 4) return(NULL)
  y  <- vapply(d, function(z) z$diff, numeric(1))
  w0 <- vapply(d, function(z) z$n,    numeric(1))
  k  <- length(y)
  mu_f <- sum(w0 * y) / sum(w0)
  Q    <- sum(w0 * (y - mu_f)^2)
  C    <- sum(w0) - sum(w0^2) / sum(w0)
  tau2 <- max(0, (Q - (k - 1)) / C)
  wi   <- 1 / (1 / w0 + tau2)
  mu   <- sum(wi * y) / sum(wi)
  se   <- sqrt(1 / sum(wi))
  p    <- 2 * stats::pnorm(-abs(mu / se))
  I2   <- if (Q > k - 1) 100 * (Q - (k - 1)) / Q else 0
  list(k = k, est = mu, p = p, i2 = I2)
}
wcor <- function(x, y, g) {
  ok <- is.finite(x)&is.finite(y); x<-x[ok]; y<-y[ok]; g<-g[ok]
  if (length(x)<10) return(c(NA,NA,0))
  for (s in unique(g)) { i<-g==s
    if (sum(i)>=2) { x[i]<-x[i]-mean(x[i]); y[i]<-y[i]-mean(y[i]) } else { x[i]<-NA; y[i]<-NA } }
  o <- is.finite(x)&is.finite(y)
  if (sum(o)<10) return(c(NA,NA,sum(o)))
  ct <- suppressWarnings(cor.test(x[o],y[o],method="spearman",exact=FALSE))
  c(unname(ct$estimate), ct$p.value, sum(o))
}

log_msg("")
log_msg("  CHALLENGE vs control (a negative estimate means the ratio FALLS under challenge)")
log_msg("  %-30s %4s %10s %9s %6s", "ratio", "k", "estimate", "p", "I2")
res <- list()
for (nm in names(RA)) {
  r <- contrast_test(RA[[nm]], CHALLENGE)
  if (is.null(r)) { log_msg("  %-30s  too few studies", nm); next }
  res[[paste(nm,"CH")]] <- r
  log_msg("  %-30s %4d %+10.3f %9.4f %5.0f%%%s", substr(nm,1,30), r$k, r$est, r$p, r$i2,
          if (r$p <= 0.05) "  *" else "")
}
log_msg("")
log_msg("  ADDITIVE vs control")
log_msg("  %-30s %4s %10s %9s %6s", "ratio", "k", "estimate", "p", "I2")
for (nm in names(RA)) {
  r <- contrast_test(RA[[nm]], ADDITIVE)
  if (is.null(r)) next
  res[[paste(nm,"AD")]] <- r
  log_msg("  %-30s %4d %+10.3f %9.4f %5.0f%%%s", substr(nm,1,30), r$k, r$est, r$p, r$i2,
          if (r$p <= 0.05) "  *" else "")
}
log_msg("")
log_msg("  Correlation with the GHI (>=3 domains, within-study)")
log_msg("  %-30s %10s %9s %6s", "ratio", "rho", "p", "arms")
for (nm in names(RA)) {
  r <- wcor(RA[[nm]], ghi3, astudy)
  if (is.na(r[1])) next
  log_msg("  %-30s %+10.3f %9.4f %6d%s", substr(nm,1,30), r[1], r[2], r[3],
          if (!is.na(r[2]) && r[2] <= 0.05) "  *" else "")
}

# leave-one-study-out on any ratio significant for challenge
# Leave-one-study-out on EVERY significant ratio, both contrasts. A result
# that survives only while one particular study is present is not a finding,
# and with k as low as 14 that is a live possibility.
log_msg("")
log_msg("  Leave-one-study-out check on every significant ratio")
any_sig <- FALSE
for (tag in c("CH","AD")) {
  tset <- if (tag == "CH") CHALLENGE else ADDITIVE
  cname <- if (tag == "CH") "challenge" else "additive"
  for (nm in names(RA)) {
    r <- res[[paste(nm, tag)]]
    if (is.null(r) || is.na(r$p) || r$p > 0.05) next
    any_sig <- TRUE
    ps <- c()
    for (s in unique(astudy)) {
      v2 <- RA[[nm]]; v2[astudy == s] <- NA
      sub <- contrast_test(v2, tset)
      if (!is.null(sub)) ps <- c(ps, sub$p)
    }
    log_msg("  %-28s (%s) p %.4f-%.4f across %d folds; <=0.05 in %d of %d",
            substr(nm,1,28), cname, min(ps), max(ps), length(ps),
            sum(ps <= 0.05), length(ps))
    if (sum(ps <= 0.05) < length(ps))
      log_msg("      NOT robust: significance depends on which studies are included.")
  }
}
if (!any_sig) log_msg("  (no ratio reached p <= 0.05 on either contrast)")

out <- c("FUNCTIONAL GUILD RATIOS", strrep("=",74), "",
  "Guilds defined two independent ways: from published poultry literature, and",
  "from measured BactoTraits phenotypes. Log-ratios tested by paired",
  "within-study differences pooled by random-effects meta-analysis.", "",
  sprintf("Arms: %d across %d studies", length(ua), length(unique(astudy))), "",
  sprintf("  %-30s %-10s %4s %10s %9s", "ratio", "contrast", "k", "estimate", "p"))
for (nm in names(res)) { r <- res[[nm]]
  out <- c(out, sprintf("  %-42s %4d %+10.3f %9.4f", nm, r$k, r$est, r$p)) }
writeLines(out, file.path(R_OUT_DIR, "gdi", "guild_ratios.txt"))
log_msg("")
log_msg("summary -> %s", file.path(R_OUT_DIR, "gdi", "guild_ratios.txt"))
log_msg("=== DONE ===")
