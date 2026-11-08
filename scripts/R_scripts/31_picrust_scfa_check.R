#!/usr/bin/env Rscript
# =============================================================================
# 31_picrust_scfa_check.R
#
# Do PICRUSt2's predicted SCFA-producing enzymes track the SCFAs actually
# measured in the same birds?
#
# RUN:  Rscript 31_picrust_scfa_check.R
# TIME: under a minute.
#
# -----------------------------------------------------------------------------
# WHY THIS RUNS BEFORE ANY MODELLING
# -----------------------------------------------------------------------------
# PICRUSt2 output is an INFERENCE. Each ASV is placed on a reference phylogeny
# and its gene content predicted from nearby genomes. Nothing is measured.
#
# This corpus offers a rare chance to check that inference directly: the GHI's
# SCFA domain holds measured acetate, propionate and butyrate, and PICRUSt2
# independently predicts the enzymes that produce those exact molecules. If
# prediction and measurement agree, the functional data carries real
# information here. If they do not, every downstream analysis would be
# modelling noise.
#
# The NSTI from the pipeline run (mean 0.057, max 0.185) says the ASVs sit
# close to reference genomes - encouraging, but that is a statement about
# phylogenetic distance, not about biological meaning in a chicken gut. This
# test addresses the latter.
#
# -----------------------------------------------------------------------------
# WHY KEGG ORTHOLOGS RATHER THAN METACYC PATHWAYS
# -----------------------------------------------------------------------------
# The MetaCyc pathway output contains 539 entries and almost nothing
# SCFA-specific: no butyrate pathway at all, and the only acetate matches were
# an aromatic-degradation pathway and methanogenesis, neither of which is
# bacterial acetate fermentation. Matching on pathway NAMES is unreliable.
#
# KEGG orthologs are individual enzymes and unambiguous. K00929 is butyrate
# kinase; it does one thing. Note the identifiers carry a "ko:" prefix in the
# PICRUSt2 output, which must be included when matching.
#
# -----------------------------------------------------------------------------
# HOW THE COMPARISON IS MADE
# -----------------------------------------------------------------------------
# Predictions are per RUN; measured SCFAs are per ARM, since published papers
# report group means. Runs are collapsed to arm means first, as everywhere
# else in this project.
#
# Correlation is computed WITHIN study and pooled, because absolute enzyme
# abundance and absolute SCFA concentration sit on study-specific scales -
# different labs, assays and units. Only the within-study relationship means
# anything.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))
con <- file(file.path(R_OUT_DIR, "gdi", "31_picrust_scfa_check.txt"), open = "wt"); sink(con, split = TRUE)

KO_FILE <- Sys.getenv("PICRUST_KO_TSV", unset = "")
if (!nzchar(KO_FILE)) {
  # The deposit ships this table gzipped; prefer whichever exists.
  .cand <- file.path(PDI_ROOT, "picrust2", "out", "KO_metagenome_out",
                     c("pred_metagenome_unstrat.tsv", "pred_metagenome_unstrat.tsv.gz"))
  .hit <- .cand[file.exists(.cand)]
  KO_FILE <- if (length(.hit)) .hit[1] else .cand[1]
}
if (!file.exists(KO_FILE)) stop("Not found: ", KO_FILE)

# enzyme -> which measured SCFA it should track
ENZ <- list(
  Butyrate = c("ko:K00929" = "butyrate kinase (buk)",
               "ko:K01034" = "butyryl-CoA:acetate CoA-transferase alpha",
               "ko:K01035" = "butyryl-CoA:acetate CoA-transferase beta",
               "ko:K00634" = "phosphotransbutyrylase (ptb)",
               "ko:K00248" = "butyryl-CoA dehydrogenase"),
  Propionate = c("ko:K01026" = "propionate CoA-transferase",
                 "ko:K00932" = "propionate kinase"),
  Acetate = c("ko:K00925" = "acetate kinase (ackA)",
              "ko:K13788" = "phosphotransacetylase (pta)")
)

log_msg("[1/4] Loading predicted KEGG ortholog abundances")
k <- fread(KO_FILE)
kid <- k[[1]]
K <- as.matrix(k[, -1, with = FALSE]); rownames(K) <- kid
log_msg("  %d KOs x %d runs", nrow(K), ncol(K))

present <- unlist(lapply(ENZ, names))
found <- present[present %in% kid]
log_msg("  %d of %d target enzymes present in the output", length(found), length(present))
if (!length(found)) stop("None of the target enzymes are present - check the ID format.")

# CLR: predicted gene abundances are compositional for the same reason taxon
# abundances are, so untransformed values would give artefactual correlations.
K[K <= 0] <- min(K[K > 0], na.rm = TRUE) / 2
gm <- exp(colMeans(log(K)))
Kc <- log(sweep(K, 2, gm, "/"))

log_msg("[2/4] Collapsing runs to arms")
md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
md <- md[!startsWith(as.character(md[[1]]), "#")]
rid <- trimws(as.character(md[[1]]))
r2s <- setNames(trimws(as.character(md$Study_ID)), rid)
r2t <- setNames(trimws(as.character(md$Treatment_Group)), rid)
runs <- intersect(colnames(Kc), names(r2s))
Kc <- Kc[, runs, drop = FALSE]
ak <- paste(r2s[runs], r2t[runs], sep = "||")
A <- t(vapply(unique(ak), function(a) rowMeans(Kc[, ak == a, drop = FALSE]), numeric(nrow(Kc))))
arm_study <- sub("\\|\\|.*$", "", rownames(A))
log_msg("  %d arms across %d studies", nrow(A), length(unique(arm_study)))

log_msg("[3/4] Attaching measured SCFAs")
bio <- as.data.frame(readxl::read_excel(METADATA_XLSX, sheet = BIOMARKER_SHEET),
                     stringsAsFactors = FALSE)
bk <- paste(trimws(bio$Study_ID), trimws(bio$Treatment_Group), sep = "||")
getm <- function(col) {
  if (!col %in% names(bio)) return(setNames(rep(NA_real_, nrow(A)), rownames(A)))
  unname(tapply(suppressWarnings(as.numeric(bio[[col]])), bk, mean, na.rm = TRUE)[rownames(A)])
}
meas <- list(Butyrate = getm("Buty"), Propionate = getm("Prop"), Acetate = getm("Ace"))
for (nm in names(meas))
  log_msg("  measured %-11s available for %d arms", nm, sum(is.finite(meas[[nm]])))

log_msg("[4/4] Within-study correlation, predicted enzyme vs measured molecule")
wcor <- function(x, y, g) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; g <- g[ok]
  if (length(x) < 8) return(c(NA_real_, 0L, 0L, NA_real_))
  for (s in unique(g)) { i <- g == s
    if (sum(i) >= 2) { x[i] <- x[i]-mean(x[i]); y[i] <- y[i]-mean(y[i]) } else { x[i] <- NA; y[i] <- NA } }
  o <- is.finite(x) & is.finite(y)
  if (sum(o) < 8) return(c(NA_real_, sum(o), length(unique(g[o])), NA_real_))
  ct <- suppressWarnings(cor.test(x[o], y[o], method = "spearman", exact = FALSE))
  c(unname(ct$estimate), sum(o), length(unique(g[o])), ct$p.value)
}

log_msg("")
log_msg("  %-11s %-44s %8s %8s %6s %8s", "SCFA", "enzyme", "rho", "p", "arms", "studies")
rows <- list()
for (nm in names(ENZ)) {
  for (id in names(ENZ[[nm]])) {
    if (!id %in% colnames(A)) next
    r <- wcor(A[, id], meas[[nm]], arm_study)
    if (is.na(r[1])) next
    rows[[paste(nm,id)]] <- c(list(scfa=nm, enz=ENZ[[nm]][[id]]), as.list(r))
    log_msg("  %-11s %-44s %+8.3f %8.4f %6d %8d%s",
            nm, substr(ENZ[[nm]][[id]], 1, 44), r[1], r[4], r[2], r[3],
            if (!is.na(r[4]) && r[4] <= 0.05) "  *" else "")
  }
}

pos <- vapply(rows, function(z) !is.na(z[[3]]) && z[[3]] > 0, logical(1))
sig <- vapply(rows, function(z) !is.na(z[[6]]) && z[[6]] <= 0.05, logical(1))
log_msg("")
log_msg("  %d of %d enzyme/molecule pairs correlate positively; %d reach p <= 0.05",
        sum(pos), length(rows), sum(sig))
log_msg("")
if (sum(sig) > 0) {
  log_msg("  READING: predicted synthesis capacity tracks the measured molecule.")
  log_msg("  The functional predictions carry real information in this corpus, so")
  log_msg("  downstream pathway analysis rests on something meaningful.")
} else if (sum(pos) > length(rows)/2) {
  log_msg("  READING: correlations lean positive but none reaches significance.")
  log_msg("  With this few arms that is weak support rather than a negative -")
  log_msg("  suggestive, not sufficient to rely on.")
} else {
  log_msg("  READING: no consistent relationship between predicted enzymes and")
  log_msg("  measured SCFAs. The predictions do not appear to capture actual")
  log_msg("  metabolic output here, and downstream functional analysis would be")
  log_msg("  modelling something unvalidated.")
}
log_msg("")
log_msg("  Measured SCFAs exist for roughly 50 arms across 14 studies, so these")
log_msg("  are indicative. A null here is weak evidence; a clear positive is")
log_msg("  stronger, since it would have to survive that limited power.")
log_msg("=== DONE ===")
sink(); close(con)
