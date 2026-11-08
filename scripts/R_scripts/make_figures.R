#!/usr/bin/env Rscript
# =============================================================================
# make_figures.R
#
# Generates every deck figure that can be produced from files already on disk.
#
# RUN:  Rscript make_figures.R
# OUT:  <repo>/figures/    one PNG per slide, named for the slide it fills
# TIME: seconds. Reads only; writes nothing outside figures/.
#
# Produces:
#   slide11_batch_correction.png   PCoA before and after batch correction
#   slide16_domain_heatmap.png     7 x 7 domain correlation matrix
#   slide18_maaslin_volcano.png    coefficient against -log10(q)
#   slide22_loso_scatter.png       GDI against GHI on held-out studies
#   slide22_loso_forest.png        per-study Spearman, sorted
#   slide19_stability.png          stability selection frequencies
#   slide23_fold_consistency.png   taxon appearance across LOSO folds
#
# Anything a source file is missing is skipped with a message rather than
# stopping the run, so a partial set still comes out.
#
# Base graphics only - no ggplot2, no pheatmap, nothing to install.
# Palette matches the deck: deep teal, seafoam, amber on white.
# =============================================================================

.locate <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
}
source(file.path(.locate(), "00_r_config.R"))

FIG <- file.path(PDI_ROOT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

INK <- "#0E333D"; TEAL <- "#2E7D8F"; SEA <- "#8FBFB4"
AMBER <- "#D4A24C"; MUTE <- "#63797F"; RULE <- "#D6E0DC"
OXIDE <- "#A34428"

png_open <- function(name, w = 2000, h = 1400, res = 200) {
  grDevices::png(file.path(FIG, name), width = w, height = h, res = res,
                 type = "cairo", bg = "white")
  par(family = "sans", mar = c(4.2, 4.6, 3.0, 1.4), col.axis = "#2A3E45",
      col.lab = INK, col.main = INK, cex.axis = 0.85, cex.lab = 0.95,
      cex.main = 1.15, font.main = 2, bty = "n", las = 1)
}
ok <- function(f) { if (file.exists(f)) TRUE else { message("  SKIP - not found: ", basename(f)); FALSE } }
done <- function(n) message("  wrote ", n)

log_msg("Writing figures to %s", FIG)

# =============================================================================
# SLIDE 11 - BATCH CORRECTION, BEFORE AND AFTER
#
# Principal coordinates on Aitchison distance (Euclidean on CLR), computed
# before and after MMUPHin. Coloured by study. If the correction worked, the
# study clusters loosen between the panels.
#
# Uses the R-side tables, which are already restricted to V3-V4 - unlike the
# QIIME artefacts in merged/, which still hold all 50 studies.
# =============================================================================
message("[slide 11] batch correction")
f_adj <- file.path(R_OUT_DIR, "genus", "genus_proportions_batch_adjusted.tsv")
f_raw <- file.path(R_OUT_DIR, "genus", "genus_counts.tsv")
if (ok(f_adj) && ok(f_raw)) {
  adj <- fread(f_adj); raw <- fread(f_raw)
  gcol_a <- names(adj)[1]; gcol_r <- names(raw)[1]
  A <- as.matrix(adj[, -1, with = FALSE]); rownames(A) <- adj[[gcol_a]]
  R0 <- as.matrix(raw[, -1, with = FALSE]); rownames(R0) <- raw[[gcol_r]]
  common <- intersect(colnames(A), colnames(R0))
  A <- A[, common, drop = FALSE]; R0 <- R0[, common, drop = FALSE]

  # study label per run, from the QIIME metadata
  md <- fread(METADATA_TSV, sep = "\t", header = TRUE, na.strings = c("", "NA"))
  md <- md[!startsWith(as.character(md[[1]]), "#")]
  s_of <- setNames(trimws(as.character(md$Study_ID)), trimws(as.character(md[[1]])))
  studies <- s_of[common]

  clr <- function(M) {
    P <- sweep(M, 2, colSums(M), "/")
    P[P <= 0] <- min(P[P > 0], na.rm = TRUE) / 2
    L <- log(P); sweep(L, 2, colMeans(L), "-")
  }
  pc <- function(M) {
    X <- t(clr(M)); X <- X[, apply(X, 2, function(v) sd(v) > 0), drop = FALSE]
    p <- prcomp(X, center = TRUE, scale. = FALSE)
    list(x = p$x[, 1], y = p$x[, 2],
         v = round(100 * (p$sdev^2 / sum(p$sdev^2))[1:2], 1))
  }
  set.seed(1)
  pb <- pc(R0); pa <- pc(A)
  us <- sort(unique(studies))
  pal <- grDevices::hcl.colors(length(us), "Zissou 1")
  cols <- grDevices::adjustcolor(pal[match(studies, us)], alpha.f = 0.65)

  png_open("slide11_batch_correction.png", w = 2400, h = 1200, res = 200)
  par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3.0, 1.0))
  plot(pb$x, pb$y, pch = 19, cex = 0.5, col = cols,
       xlab = sprintf("PC1  (%.1f%%)", pb$v[1]), ylab = sprintf("PC2  (%.1f%%)", pb$v[2]),
       main = "Before correction")
  mtext("samples cluster by study", side = 3, line = 0.1, cex = 0.78, col = MUTE)
  plot(pa$x, pa$y, pch = 19, cex = 0.5, col = cols,
       xlab = sprintf("PC1  (%.1f%%)", pa$v[1]), ylab = sprintf("PC2  (%.1f%%)", pa$v[2]),
       main = "After MMUPHin")
  mtext(sprintf("%d studies, %d runs, coloured by study", length(us), length(common)),
        side = 3, line = 0.1, cex = 0.78, col = MUTE)
  dev.off(); done("slide11_batch_correction.png")
}

# =============================================================================
# SLIDE 16 - DOMAIN CORRELATION HEATMAP
#
# Within-study correlations between the seven GHI domains. The point of the
# figure is that the domains form three groups rather than one construct, so
# the cells are annotated with r and the diagonal is left blank.
# =============================================================================
message("[slide 16] domain heatmap")
f_dom <- file.path(R_OUT_DIR, "ghi", "ghi_domain_correlation_matrix.tsv")
if (ok(f_dom)) {
  d <- fread(f_dom)
  nm <- d[[1]]
  M <- as.matrix(d[, -1, with = FALSE]); rownames(M) <- nm; colnames(M) <- nm
  diag(M) <- NA
  short <- gsub("Immunoglobulin_", "Ig ", nm)

  png_open("slide16_domain_heatmap.png", w = 1700, h = 1500, res = 200)
  par(mar = c(8.2, 8.2, 2.4, 1.0))
  n <- nrow(M)
  ramp <- colorRampPalette(c(OXIDE, "#FFFFFF", TEAL))(101)
  idx <- function(v) if (is.na(v)) NA else ramp[round((max(-1, min(1, v)) + 1) / 2 * 100) + 1]
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = "Within-study correlation between GHI domains", asp = 1)
  for (i in 1:n) for (j in 1:n) {
    v <- M[i, j]
    rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5,
         col = if (is.na(v)) "#F2F5F3" else idx(v), border = "white", lwd = 2)
    if (!is.na(v))
      text(j, n - i + 1, sprintf("%+.2f", v), cex = 0.72,
           col = if (abs(v) > 0.55) "white" else INK, font = if (abs(v) > 0.55) 2 else 1)
  }
  axis(1, at = 1:n, labels = short, las = 2, tick = FALSE, cex.axis = 0.78, line = -0.4)
  axis(2, at = n:1, labels = short, las = 1, tick = FALSE, cex.axis = 0.78, line = -0.4)
  mtext("teal = positive, red = negative, blank = diagonal", side = 3, line = 0.4,
        cex = 0.74, col = MUTE)
  dev.off(); done("slide16_domain_heatmap.png")
}

# =============================================================================
# SLIDE 18 - MAASLIN2 VOLCANO
#
# Every genus tested, coefficient against significance. The dashed line is the
# FDR 0.25 screening threshold; genera above it went forward to the LASSO.
# =============================================================================
message("[slide 18] volcano")
f_ma <- file.path(R_OUT_DIR, "gdi", "maaslin2_all_results.tsv")
if (ok(f_ma)) {
  m <- fread(f_ma)
  cf <- intersect(c("coef", "coefficient"), names(m))[1]
  qv <- intersect(c("qval", "qvalue", "q.value"), names(m))[1]
  ft <- intersect(c("feature", "genus", "Feature"), names(m))[1]
  if (is.na(cf) || is.na(qv)) {
    message("  SKIP - expected coef/qval columns not present")
  } else {
    m <- m[is.finite(get(cf)) & is.finite(get(qv))]
    y <- -log10(pmax(m[[qv]], 1e-12)); x <- m[[cf]]
    sig <- m[[qv]] <= 0.25
    png_open("slide18_maaslin_volcano.png", w = 2000, h = 1400, res = 200)
    plot(x, y, pch = 19, cex = ifelse(sig, 1.15, 0.75),
         col = ifelse(sig, grDevices::adjustcolor(TEAL, .85), grDevices::adjustcolor(MUTE, .4)),
         xlab = "Association with the GHI  (coefficient)",
         ylab = expression(-log[10] * "(q)"),
         main = "Differential abundance against the Gut Health Index")
    abline(h = -log10(0.25), lty = 2, col = AMBER, lwd = 2)
    abline(v = 0, col = RULE, lwd = 1)
    text(par("usr")[2], -log10(0.25), "  FDR 0.25", adj = c(1.02, -0.6),
         col = AMBER, cex = 0.78, font = 2)
    if (!is.na(ft)) {
      top <- m[order(get(qv))][1:min(10, .N)]
      ty <- -log10(pmax(top[[qv]], 1e-12))
      text(top[[cf]], ty, gsub("^g__", "", top[[ft]]), pos = 4, cex = 0.66,
           col = INK, offset = 0.35, xpd = NA)
    }
    mtext(sprintf("%d of %d genera passed the screening threshold", sum(sig), nrow(m)),
          side = 3, line = 0.25, cex = 0.78, col = MUTE)
    dev.off(); done("slide18_maaslin_volcano.png")
  }
}

# =============================================================================
# SLIDE 22 - LOSO SCATTER AND FOREST
#
# Both panels come from held-out predictions only. A working index gives a
# NEGATIVE slope: the GDI rises with dysbiosis, the GHI with health.
# =============================================================================
message("[slide 22] LOSO")
f_lp <- file.path(R_OUT_DIR, "gdi", "loso_predictions.tsv")
if (ok(f_lp)) {
  lp <- fread(f_lp)
  if (all(c("gdi_centered", "ghi_centered", "study") %in% names(lp))) {
    d <- lp[is.finite(gdi_centered) & is.finite(ghi_centered)]
    us <- sort(unique(d$study))
    pal <- grDevices::hcl.colors(length(us), "Zissou 1")
    png_open("slide22_loso_scatter.png", w = 2000, h = 1500, res = 200)
    plot(d$gdi_centered, d$ghi_centered, pch = 19, cex = 1.05,
         col = grDevices::adjustcolor(pal[match(d$study, us)], .8),
         xlab = "GDI  (higher = more dysbiotic)",
         ylab = "GHI, within-study centred  (higher = healthier)",
         main = "Leave-one-study-out predictions")
    abline(lm(ghi_centered ~ gdi_centered, data = d), col = INK, lwd = 2.4)
    abline(h = 0, v = 0, col = RULE, lwd = 1)
    rho <- suppressWarnings(cor(d$gdi_centered, d$ghi_centered, method = "spearman"))
    legend("topright", bty = "n", text.col = INK, cex = 0.9,
           legend = c(sprintf("Spearman rho = %+.3f", rho),
                      sprintf("%d arms, %d held-out studies", nrow(d), length(us))))
    mtext("each point is one treatment arm, coloured by study; every point scored by a model that never saw its study",
          side = 3, line = 0.25, cex = 0.72, col = MUTE)
    dev.off(); done("slide22_loso_scatter.png")
  } else message("  SKIP - expected columns not present")
}

f_ps <- file.path(R_OUT_DIR, "gdi", "loso_per_study.tsv")
if (ok(f_ps)) {
  ps <- fread(f_ps)
  ps <- ps[is.finite(rho)][order(rho)]
  if (nrow(ps)) {
    png_open("slide22_loso_forest.png", w = 1700, h = 1900, res = 200)
    par(mar = c(4.4, 9.6, 3.0, 1.6))
    n <- nrow(ps)
    plot(NA, xlim = c(-1.05, 1.05), ylim = c(0.5, n + 0.5), yaxt = "n",
         xlab = "Spearman rho  (negative = working direction)", ylab = "",
         main = "Per-study performance, held out")
    abline(v = 0, col = INK, lwd = 1.6)
    for (v in c(-0.5, 0.5)) abline(v = v, col = RULE, lty = 3)
    segments(0, 1:n, ps$rho, 1:n, col = RULE, lwd = 1.6)
    points(ps$rho, 1:n, pch = 19, cex = 1.15,
           col = ifelse(ps$rho < 0, TEAL, AMBER))
    axis(2, at = 1:n, labels = ps$study, las = 1, tick = FALSE, cex.axis = 0.72, line = -0.5)
    legend("bottomright", bty = "n", cex = 0.8, pch = 19, col = c(TEAL, AMBER),
           text.col = INK, legend = c("working direction", "against direction"))
    mtext(sprintf("%d of %d studies negative; per-study values rest on 2-7 arms and are shown for pattern, not inference",
                  sum(ps$rho < 0), n), side = 3, line = 0.25, cex = 0.7, col = MUTE)
    dev.off(); done("slide22_loso_forest.png")
  }
}

# =============================================================================
# SLIDE 19 - STABILITY SELECTION
# The registered 70% line is drawn so the shortfall is visible rather than
# asserted.
# =============================================================================
message("[slide 19] stability")
f_st <- file.path(R_OUT_DIR, "gdi", "gdi_stability.tsv")
if (ok(f_st)) {
  st <- fread(f_st)[order(-stability)][1:min(12, .N)]
  st <- st[stability > 0]
  png_open("slide19_stability.png", w = 2000, h = 1400, res = 200)
  par(mar = c(4.4, 11.5, 3.0, 2.2))
  b <- barplot(rev(st$stability * 100), horiz = TRUE, col = rev(ifelse(st$stability >= 0.20, TEAL, SEA)),
               border = NA, xlim = c(0, 100), xlab = "Selection frequency across 1,000 resamples (%)",
               names.arg = rev(gsub("^g__", "", st$genus)), las = 1, cex.names = 0.78,
               main = "Stability selection")
  abline(v = 70, col = AMBER, lwd = 2.4, lty = 2)
  text(70, max(b) + 0.6, "registered threshold", col = AMBER, cex = 0.74, font = 2, adj = c(0.5, -0.3), xpd = NA)
  text(rev(st$stability * 100), b, sprintf(" %.0f%%", rev(st$stability * 100)),
       pos = 4, cex = 0.72, col = INK, xpd = NA)
  mtext("no taxon reached 70%; the panel requires 20% stability AND a non-zero coefficient in the full-data fit",
        side = 3, line = 0.25, cex = 0.72, col = MUTE)
  dev.off(); done("slide19_stability.png")
}

# =============================================================================
# SLIDE 23 - FOLD CONSISTENCY
# Recomputed from the fold panels if present; otherwise skipped, since the
# summary text file does not carry the per-taxon counts in parseable form.
# =============================================================================
message("[slide 23] fold consistency")
f_def <- file.path(R_OUT_DIR, "gdi", "gdi_index_definition.tsv")
f_sum <- file.path(R_OUT_DIR, "gdi", "loso_summary.txt")
if (ok(f_def) && ok(f_sum)) {
  panel <- unique(fread(f_def)$genus)
  tx <- readLines(f_sum, warn = FALSE)
  hits <- grep("^\\s+g__\\S+\\s+\\d+/\\s*\\d+", tx, value = TRUE)
  if (!length(hits)) {
    message("  SKIP - fold counts not present in loso_summary.txt; take them from the 08 console output")
  } else {
    nmv <- sub("^\\s+(g__\\S+).*", "\\1", hits)
    num <- as.numeric(sub(".*?(\\d+)\\s*/\\s*(\\d+).*", "\\1", hits))
    den <- as.numeric(sub(".*?(\\d+)\\s*/\\s*(\\d+).*", "\\2", hits))
    o <- order(num)
    png_open("slide23_fold_consistency.png", w = 2000, h = 1300, res = 200)
    par(mar = c(4.4, 11.5, 3.0, 2.2))
    b <- barplot(num[o], horiz = TRUE, col = ifelse(num[o] / den[o] >= 0.7, TEAL, SEA),
                 border = NA, xlim = c(0, max(den) * 1.15),
                 xlab = sprintf("Folds selected (of %d)", max(den)),
                 names.arg = gsub("^g__", "", nmv[o]), las = 1, cex.names = 0.8,
                 main = "Panel consistency across leave-one-study-out folds")
    text(num[o], b, sprintf(" %d", num[o]), pos = 4, cex = 0.74, col = INK, xpd = NA)
    dev.off(); done("slide23_fold_consistency.png")
  }
}

log_msg("")
log_msg("Done. Figures in %s", FIG)
log_msg("Not generated here: the PRISMA flow (build from the Rayyan export) and")
log_msg("the diversity panels (05_diversity.sh has not produced output).")
