# 65_multiplicity.R -----------------------------------------------------------
# Benjamini-Hochberg correction applied WITHIN each family of tests, so that
# every reported p-value is accompanied by its adjusted value.
#
# Families are defined by what was tested together, not by what was
# significant. Reading the inputs from the files the other scripts write means
# the correction cannot drift out of step with the analyses.
#
# EXPECTED: 5 families, 36 tests, 12 nominally significant -> 9 survive q <= 0.05.
#   survives: six of eight axis-vs-coverage, both axis-by-target, caecal butyrate (q=0.038)
#   fails:    acetate 0.118, anti-inflammatory coupling 0.110,
#             lactic/opportunist guild 0.136, anaerobe/facultative 0.290
# If fewer than 5 families are listed, a parse failed -- check the family names
# printed below against the five expected.
# -----------------------------------------------------------------------------
source("50_setup.R")
con <- open_log("65_multiplicity")

fam <- list()

# --- family 1: functional guild ratios (14 tests), parsed from the archive ---
gf <- file.path(DIR_GDI, "guild_ratios.txt")
if (file.exists(gf)) {
  L <- readLines(gf, warn = FALSE)
  L <- L[grepl("^\\s+(LIT|TRAIT):", L)]
  parts <- regmatches(L, regexpr("[-+][0-9.]+\\s+[0-9.]+\\s*$", L))
  lab <- trimws(sub("\\s+[-+][0-9.]+\\s+[0-9.]+\\s*$", "", L))
  num <- do.call(rbind, lapply(strsplit(trimws(parts), "\\s+"), as.numeric))
  fam[["Guild ratios"]] <- data.frame(test = lab, p = num[, 2])
}

# --- family 2: displacement domain couplings (7 tests) ----------------------
df <- file.path(DIR_OUT, "displacement_coupling.csv")
if (file.exists(df)) {
  d <- read.csv(df)
  fam[["Displacement couplings"]] <- data.frame(test = d$domain, p = d$p)
}

# --- family 3: caecal SCFA (3 tests) ----------------------------------------
sf <- file.path(DIR_OUT, "scfa_guild.csv")
if (file.exists(sf)) {
  d <- read.csv(sf)
  fam[["Caecal SCFA guild"]] <- data.frame(test = d$scfa, p = d$p)
}

# --- family 4: unsupervised axis vs GHI at increasing coverage (8 tests) ----
# Columns are: target, ">=", N, rho, perm p, arms, studies. The p is field 5.
af <- file.path(DIR_GDI, "axis_vs_coverage.txt")
if (file.exists(af)) {
  L <- grep("^\\s+(Combined|Biomarker|Morphology|Performance)\\s+>=", readLines(af, warn = FALSE), value = TRUE)
  f <- strsplit(trimws(L), "\\s+")
  pv <- suppressWarnings(as.numeric(sapply(f, `[`, 5)))
  lb <- sapply(f, function(z) paste(z[1], z[2], z[3]))
  if (any(is.finite(pv)))
    fam[["Axis vs coverage"]] <- data.frame(test = lb[is.finite(pv)], p = pv[is.finite(pv)])
}

# --- family 5: unsupervised axis by target (4 tests) ------------------------
# From the UNSUPERVISED block of target_comparison_domains.txt:
# columns are target, rho, perm p, n -- the p is field 3.
tf <- file.path(DIR_GDI, "target_comparison_domains.txt")
if (file.exists(tf)) {
  L <- readLines(tf, warn = FALSE)
  start <- grep("^UNSUPERVISED", L)
  if (length(start)) {
    blk <- L[seq(start[1], length(L))]
    hit <- grep("^\\s+(Morphology|Biomarker|Performance|Combined)\\s+[-+]", blk, value = TRUE)
    f <- strsplit(trimws(hit), "\\s+")
    pv <- suppressWarnings(as.numeric(sapply(f, `[`, 3)))
    lb <- sapply(f, `[`, 1)
    if (any(is.finite(pv)))
      fam[["Unsupervised axis by target"]] <- data.frame(test = lb[is.finite(pv)], p = pv[is.finite(pv)])
  }
}

cat("families found:", length(fam), "--", paste(names(fam), collapse = "; "), "\n")
if (length(fam) < 5) warning("Expected 5 families; a parser may have failed.")
if (length(fam) == 0) stop("No input families found. Run scripts 57 and 59 first.")

hdr("Benjamini-Hochberg within each family")
out <- list()
for (nm in names(fam)) {
  d <- fam[[nm]]
  d <- d[is.finite(d$p), ]
  d$family <- nm
  d$BH_q <- p.adjust(d$p, method = "BH")
  d$survives <- d$BH_q <= 0.05
  d <- d[order(d$p), c("family", "test", "p", "BH_q", "survives")]
  cat("\n--", nm, "(", nrow(d), "tests )\n")
  print(d, row.names = FALSE, digits = 3)
  out[[nm]] <- d
}
res <- do.call(rbind, out)
write.csv(res, file.path(DIR_OUT, "multiplicity_corrected_pvalues.csv"), row.names = FALSE)

hdr("Summary")
cat("nominally significant (p <= 0.05):", sum(res$p <= 0.05), "\n")
cat("surviving BH within family (q <= 0.05):", sum(res$survives), "\n")
cat("\nReport both values throughout. Findings that do not survive are labelled\n",
    "exploratory; they are not deleted.\n", sep = "")
close_log(con)
