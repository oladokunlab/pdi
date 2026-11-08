#!/usr/bin/env Rscript
# 04f_ghi_domain_table.R
# Manuscript Table 3: response of each of the seven GHI domains and the
# combined index. Pooling as in 04e: per-study paired differences, then a
# DerSimonian-Laird random-effects model weighted by the number of group
# means. Massacci_2019 has no concurrent control and is left out of the
# challenge contrast. Reads ghi_scores_pca.tsv written by 04b.
.locate_and_source_config <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  d <- if (length(f) > 0) dirname(normalizePath(sub("^--file=", "", f))) else normalizePath(".")
  source(file.path(d, "00_r_config.R"))
}
.locate_and_source_config()

g <- utils::read.delim(file.path(R_OUT_DIR, "ghi", "ghi_scores_pca.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
HEALTHY   <- "Unchallenged control"
ADDITIVE  <- "Non-disruptive additive"
CHALLENGE <- c("Dietary challenge", "Infectious challenge", "Toxin challenge",
               "Antibiotic challenge", "Environmental challenge")
SKIP_CH   <- "Massacci_2019"
DOMS <- c("Morphology" = "GHI_dom_Morphology",
          "Barrier" = "GHI_dom_Barrier",
          "Inflammation" = "GHI_dom_Inflammation",
          "Anti-inflammatory" = "GHI_dom_Anti_Inflammatory",
          "SCFA" = "GHI_dom_SCFA",
          "Performance" = "GHI_dom_Performance",
          "Systemic immunoglobulin" = "GHI_dom_Immunoglobulin_Systemic",
          "Combined index" = "GHI_Combined_raw")

diffs <- function(x, better, worse, skip = character(0)) {
  out <- data.frame(study = character(0), diff = numeric(0), n = numeric(0))
  for (s in unique(g$Study_ID)) {
    if (s %in% skip) next
    k <- g$Study_ID == s & !is.na(x)
    a <- x[k & g$Intervention_Class %in% better]
    b <- x[k & g$Intervention_Class %in% worse]
    if (length(a) && length(b))
      out <- rbind(out, data.frame(study = s, diff = mean(a) - mean(b), n = length(a) + length(b)))
  }
  out
}
pool <- function(d) {
  y <- d$diff; w0 <- d$n; k <- length(y)
  mu_f <- sum(w0 * y) / sum(w0); Q <- sum(w0 * (y - mu_f)^2)
  C <- sum(w0) - sum(w0^2) / sum(w0); tau2 <- max(0, (Q - (k - 1)) / C)
  wi <- 1 / (1 / w0 + tau2); mu <- sum(wi * y) / sum(wi); se <- sqrt(1 / sum(wi))
  q <- stats::qt(0.975, max(1, k - 2))
  list(k = k, est = mu, lo = mu - 1.96 * se, hi = mu + 1.96 * se,
       I2 = if (Q > k - 1) 100 * (Q - (k - 1)) / Q else 0,
       pi_lo = mu - q * sqrt(tau2 + se^2), pi_hi = mu + q * sqrt(tau2 + se^2))
}
auc <- function(x) {
  m <- g$Intervention_Class %in% c(HEALTHY, CHALLENGE) & !is.na(x)
  y <- g$Intervention_Class[m] %in% HEALTHY; r <- rank(x[m])
  n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

rows <- list()
for (nm in names(DOMS)) {
  x <- g[[DOMS[[nm]]]]
  a <- pool(diffs(x, ADDITIVE, HEALTHY))
  dc <- diffs(x, HEALTHY, CHALLENGE, SKIP_CH); c <- pool(dc)
  s <- sum(dc$diff > 0); n <- nrow(dc)
  p <- stats::binom.test(s, n, 0.5, alternative = "greater")$p.value
  rows[[nm]] <- data.frame(domain = nm, groups = sum(!is.na(x)), above = s, studies = n, P = p,
    add_est = a$est, add_lo = a$lo, add_hi = a$hi, add_k = a$k, add_I2 = a$I2, add_pi_lo = a$pi_lo, add_pi_hi = a$pi_hi,
    ch_est = c$est, ch_lo = c$lo, ch_hi = c$hi, ch_k = c$k, ch_I2 = c$I2, ch_pi_lo = c$pi_lo, ch_pi_hi = c$pi_hi,
    AUC = auc(x))
}
tab <- do.call(rbind, rows)
utils::write.table(tab, file.path(R_OUT_DIR, "ghi", "ghi_domain_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
lines <- c("=== Table 3: GHI domains and combined index ===", "",
  sprintf("%-24s %6s %9s %7s %27s %3s %4s %27s %3s %4s %6s", "domain", "groups", "ctrl>chal", "P",
          "additive (95% CI)", "k", "I2", "challenge (95% CI)", "k", "I2", "AUC"))
for (i in seq_len(nrow(tab))) with(tab[i, ], {
  lines <<- c(lines, sprintf("%-24s %6d %4d of %-2d %7.3f %+8.3f (%+.3f to %+.3f) %3d %4.0f %+8.3f (%+.3f to %+.3f) %3d %4.0f %6.3f",
    domain, groups, above, studies, P, add_est, add_lo, add_hi, add_k, add_I2, ch_est, ch_lo, ch_hi, ch_k, ch_I2, AUC))
})
writeLines(lines, file.path(R_OUT_DIR, "ghi", "ghi_domain_table.txt"))
writeLines(lines)
