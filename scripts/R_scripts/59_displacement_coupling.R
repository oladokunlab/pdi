# 59_displacement_coupling.R --------------------------------------------------
# Does a change in the community track a change in gut health?
# Replaces script 42, which computed these models but wrote no output file.
#
# Two predictors, and the distinction is the result:
#   tax_disp    unsigned distance from the arm's own control centroid
#   signed_axis position along the control-to-challenge direction
#
# EXPECTED (87 arms, 37 studies):
#   unsigned  +0.040  p = 0.504     signed  -0.186  p = 0.0399
#   Anti_Inflammatory -0.402 p = 0.0157 (BH q = 0.110 across the 7 domains)
# -----------------------------------------------------------------------------
source("50_setup.R")
suppressMessages(library(lmerTest))
con <- open_log("59_displacement_coupling")

df <- read.delim(file.path(DIR_GDI, "displacement_coupling_data.tsv"), check.names = FALSE)
cat("arms:", nrow(df), " studies:", length(unique(df$Study_ID)), "\n")

pfrom <- function(m) { co <- summary(m)$coefficients; c(est = co[2, 1], p = co[2, ncol(co)]) }

hdr("Composite GHI displacement")
m1 <- lmer(ghi_disp ~ tax_disp + (1 | Study_ID), data = df)
m2 <- lmer(ghi_disp ~ signed_axis + (1 | Study_ID), data = df[!is.na(df$signed_axis), ])
a <- pfrom(m1); b <- pfrom(m2)
cat(sprintf("unsigned displacement : coef %+.4f  p = %.4f\n", a["est"], a["p"]))
cat(sprintf("signed challenge axis : coef %+.4f  p = %.4f\n", b["est"], b["p"]))
cat("\nHow far the community moved says nothing; which direction it moved does.\n")

hdr("Per domain, signed axis as predictor")
rows <- list()
for (dc in GHI_DOM_COLS) {
  d3 <- df[!is.na(df[[dc]]) & !is.na(df$signed_axis), ]
  if (nrow(d3) < 12 || length(unique(d3$Study_ID)) < 3) next
  m <- tryCatch(lmer(as.formula(sprintf("`%s` ~ signed_axis + (1 | Study_ID)", dc)), data = d3),
                error = function(e) NULL)
  if (is.null(m)) next
  e <- pfrom(m)
  rows[[dc]] <- data.frame(domain = sub("^GHI_dom_", "", dc), n_arms = nrow(d3),
                           n_studies = length(unique(d3$Study_ID)),
                           coef = unname(e["est"]), p = unname(e["p"]))
}
res <- do.call(rbind, rows)
res <- res[order(res$p), ]
res$BH_q <- p.adjust(res$p, method = "BH")
print(res, row.names = FALSE, digits = 3)
write.csv(res, file.path(DIR_OUT, "displacement_coupling.csv"), row.names = FALSE)
cat(sprintf("\n%d of %d domains share the negative sign.\n", sum(res$coef < 0), nrow(res)))
close_log(con)
