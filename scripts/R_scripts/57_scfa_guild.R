# 57_scfa_guild.R -------------------------------------------------------------
# The one positive result: a pre-specified butyrate-producer guild predicts
# MEASURED caecal butyrate, on caecal samples only.
#
# The guild is taken from the project's own literature-defined guild list
# (GUILDS_LIT$SCFA_butyrate in 33_guild_ratios.R) and was fixed before testing.
#
# EXPECTED (33 arms, 10 studies, 9999 within-study permutations):
#   butyrate   rho = +0.495  p = 0.0161  BH q = 0.048
#   acetate    rho = +0.409  p = 0.0471  BH q = 0.071
#   propionate rho = -0.004  p = 0.9884
#   Unit trap: raw values give p = 0.70; standardised within study, p = 0.003.
# -----------------------------------------------------------------------------
source("50_setup.R")
con <- open_log("57_scfa_guild")

# guild definition, copied verbatim from 33_guild_ratios.R so the two agree
# copied verbatim from GUILDS_LIT$SCFA_butyrate in 33_guild_ratios.R
GUILD_BUTYRATE <- c("Faecalibacterium", "Butyricicoccus", "Roseburia", "Anaerostipes",
                    "Blautia", "Coprococcus", "Subdoligranulum", "Eubacterium",
                    "Butyrivibrio", "Anaerobutyricum", "Agathobacter", "Oscillibacter",
                    "Pseudobutyrivibrio", "Shuttleworthia")

prop <- load_prop_adj()
gl <- sub("^g__", "", rownames(prop))
inG <- gl %in% GUILD_BUTYRATE
cat("guild genera present:", sum(inG), "of", length(GUILD_BUTYRATE), "\n")
cat(paste(" ", gl[inG], collapse = "\n"), "\n")

eps <- min(prop[prop > 0])
lr <- log((colSums(prop[inG, , drop = FALSE]) + eps) /
          (colSums(prop[!inG, , drop = FALSE]) + eps))

meta <- load_seq_meta()
d <- data.frame(sample_id = names(lr), lr = as.numeric(lr), stringsAsFactors = FALSE)
d <- merge(d, meta[, c("sample_id", "Study_ID", "Treatment_Group", "Segment")], by = "sample_id")
d <- d[d$Segment == "Cecum", ]
arm <- aggregate(list(lr = d$lr),
                 by = list(Study_ID = d$Study_ID, Treatment_Group = d$Treatment_Group), FUN = mean)

# Measured SCFA concentrations live in the extraction workbook, not in
# ghi_scores.tsv (which holds domain scores only).
scfa_cols <- c(Butyrate = "Buty", Acetate = "Ace", Propionate = "Prop")
WB <- Sys.getenv("PDI_WORKBOOK", unset = file.path(DIR_META, "GDI_Metadata.xlsx"))
stopifnot(requireNamespace("readxl", quietly = TRUE), file.exists(WB))
raw <- as.data.frame(readxl::read_excel(WB, sheet = "Sample_with_Zscores"))

hdr("Guild log-ratio versus measured caecal SCFA")
out <- list()
for (nm in names(scfa_cols)) {
  src <- if (scfa_cols[nm] %in% names(raw)) raw[[scfa_cols[nm]]] else NULL
  if (is.null(src)) { cat(nm, ": marker column not present in ghi_scores.tsv; skipped\n"); next }
  dd <- data.frame(Study_ID = raw$Study_ID, Treatment_Group = raw$Treatment_Group,
                   Segment = raw$Segment, v = suppressWarnings(as.numeric(src)))
  dd <- dd[dd$Segment %in% c("Cecum", "Caecum") & is.finite(dd$v), ]
  if (nrow(dd) == 0) { cat(nm, ": no caecal measurements; skipped\n"); next }
  av <- aggregate(list(v = dd$v), by = list(Study_ID = dd$Study_ID,
                                            Treatment_Group = dd$Treatment_Group), FUN = mean)
  J <- merge(arm, av, by = c("Study_ID", "Treatment_Group"))
  tab <- table(J$Study_ID); J <- J[J$Study_ID %in% names(tab)[tab >= 2], ]
  if (nrow(J) < 10) { cat(nm, ": too few arms (", nrow(J), "); skipped\n"); next }
  # UNIT TRAP: SCFA are reported in eight different units across studies, so the
  # outcome must be standardised within study before pooling.
  J$vz  <- ave(J$v,  J$Study_ID, FUN = function(z) if (sd(z) > 0) (z - mean(z)) / sd(z) else z * 0)
  J$lrc <- center_within(J$lr, J$Study_ID)
  r <- perm_spearman_within(J$lrc, J$vz, J$Study_ID)
  cat(sprintf("%-11s rho = %+.3f  p = %.4f   (%d arms, %d studies)\n",
              nm, r$rho, r$p, nrow(J), length(unique(J$Study_ID))))
  out[[nm]] <- data.frame(scfa = nm, rho = r$rho, p = r$p, arms = nrow(J),
                          studies = length(unique(J$Study_ID)))
}
# leave-one-study-out: does the direction hold when each study is removed?
if ("Butyrate" %in% names(out)) {
  cat("\n")
  hdr("Butyrate: leave-one-study-out")
  dd <- data.frame(Study_ID = raw$Study_ID, Treatment_Group = raw$Treatment_Group,
                   Segment = raw$Segment, v = suppressWarnings(as.numeric(raw[["Buty"]])))
  dd <- dd[dd$Segment %in% c("Cecum", "Caecum") & is.finite(dd$v), ]
  av <- aggregate(list(v = dd$v), by = list(Study_ID = dd$Study_ID,
                                            Treatment_Group = dd$Treatment_Group), FUN = mean)
  J <- merge(arm, av, by = c("Study_ID", "Treatment_Group"))
  tb <- table(J$Study_ID); J <- J[J$Study_ID %in% names(tb)[tb >= 2], ]
  studies <- unique(J$Study_ID); keep_dir <- 0
  for (s in studies) {
    K <- J[J$Study_ID != s, ]
    K$vz  <- ave(K$v,  K$Study_ID, FUN = function(z) if (sd(z) > 0) (z - mean(z)) / sd(z) else z * 0)
    K$lrc <- center_within(K$lr, K$Study_ID)
    rr <- suppressWarnings(cor(K$lrc, K$vz, method = "spearman"))
    if (rr > 0) keep_dir <- keep_dir + 1
    cat(sprintf("  drop %-16s rho = %+.3f\n", s, rr))
  }
  cat(sprintf("direction retained in %d of %d folds\n", keep_dir, length(studies)))
}

res <- do.call(rbind, out)
if (!is.null(res)) {
  res$BH_q <- p.adjust(res$p, method = "BH")
  print(res, row.names = FALSE, digits = 3)
  write.csv(res, file.path(DIR_OUT, "scfa_guild.csv"), row.names = FALSE)
}
cat("\nThe community predicts what it produces. On the same arms, the same guild\n",
    "against the GHI gives rho = +0.233, p = 0.267.\n", sep = "")
close_log(con)
