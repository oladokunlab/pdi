#!/bin/bash
# preflight_and_release.sh
#
# Runs the control checks, quarantines scratch files, builds both release
# archives, then verifies the built tarball. Ends with a single PASS/FAIL.
# If it says PASS you can upload without waiting for anything else.
#
# Usage:
#   export PDI_ROOT=/scratch/user/$USER/pdi
#   module load GCC/14.3.0 R/4.5.2
#   bash preflight_and_release.sh 2>&1 | tee ~/preflight_20260927.log

set -u
: "${PDI_ROOT:?export PDI_ROOT=/scratch/user/$USER/pdi}"
RELDIR=/scratch/user/$USER/releases
BUILDER=/scratch/user/$USER/pdi_release_builder
BACKUP=/scratch/user/$USER/pdi_backups
GEN="$PDI_ROOT/r_analysis/genus"
VER="$PDI_ROOT/r_analysis/verification"
GDI="$PDI_ROOT/r_analysis/gdi"
FAIL=0
note() { echo "  $*"; }
bad()  { echo "  FAIL: $*"; FAIL=1; }

echo "=============================================================="
echo "1. CONTROL: configs B, C, E must be unchanged"
echo "=============================================================="
# B, C and E build CLR from raw counts, so the batch-correction rebuild cannot
# touch them. Note what that does and does not cover. Their GHI_R2 is invariant
# to a change in the MATRIX but not to a change in the INDEX -- the workbook
# corrections moved the GHI, so these three moved with it (0.00163 -> 0.00157,
# 0.00152 -> 0.00149, 0.00425 -> 0.00427) and are asserted at their post-
# correction values. The stronger invariant, which holds across both a matrix
# rebuild and an index change, is Age_R2: it does not involve the GHI at all,
# so if any of the three Age_R2 values moves, something is wrong that neither
# the rebuild nor the corrections can explain. GHI_p is a Monte Carlo estimate
# over 9,999 permutations and moves slightly every run, so it is printed but
# NOT asserted.
cat > /tmp/pdi_control.R <<'RCODE'
f <- file.path(Sys.getenv("VER"), "permanova_five_configs.csv")
if (!file.exists(f)) { cat("MISSING", f, "\n"); quit(status = 3) }
d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
exp <- list("B"=0.00157, "C"=0.00149, "E"=0.00427)
expage <- list("B"=0.13932, "C"=0.13688, "E"=0.08345)
bad <- 0
for (k in names(exp)) {
  r <- d[substr(d$config, 1, 1) == k, ]
  if (nrow(r) != 1) { cat(sprintf("  FAIL %s: %d rows\n", k, nrow(r))); bad <- 1; next }
  ok  <- isTRUE(all.equal(r$GHI_R2, exp[[k]], tolerance = 1e-9))
  oka <- isTRUE(all.equal(r$Age_R2, expage[[k]], tolerance = 1e-9))
  cat(sprintf("  %s  GHI_R2 %.5f (expect %.5f) %s | Age_R2 %.5f %s | GHI_p %.4f (not asserted)\n",
              k, r$GHI_R2, exp[[k]], if (ok) "ok" else "MOVED",
              r$Age_R2, if (oka) "ok" else "MOVED", r$GHI_p))
  if (!ok || !oka) bad <- 1
}
cat("\n  configs A and D (these SHOULD have moved):\n")
for (k in c("A","D")) {
  r <- d[substr(d$config, 1, 1) == k, ]
  cat(sprintf("  %s  GHI_R2 %.5f  GHI_p %.4f  Age_R2 %.5f\n", k, r$GHI_R2, r$GHI_p, r$Age_R2))
}
quit(status = if (bad) 1 else 0)
RCODE
VER="$VER" Rscript /tmp/pdi_control.R
case $? in
  0) note "control PASSED - the rebuild is confined to the batch-corrected matrix" ;;
  1) bad "configs B/C/E MOVED - the change is not confined. Do not release." ;;
  *) bad "could not read permanova_five_configs.csv" ;;
esac

echo
echo "=============================================================="
echo "2. Headline numbers after the re-run"
echo "=============================================================="
echo "-- 51_batch_correction_check (deck's before/after table) --"
grep -E 'before:|after :|reduction' "$VER/51_batch_correction_check.txt" 2>/dev/null | sed 's/^/  /' || bad "51 output missing"
echo "-- 27_runlevel_gdi --"
grep -E 'Runs:|Studies contributing|In-sample AUC:|LOSO AUC:|at 70% stability' "$GDI/runlevel_gdi.txt" 2>/dev/null | sed 's/^/  /' || bad "27 output missing"
echo "-- 28_runlevel_groupings --"
sed -n '/grouping /,/^$/p' "$GDI/runlevel_groupings.txt" 2>/dev/null | sed 's/^/  /' || bad "28 output missing"
grep -A1 'Reference - control vs challenge' "$GDI/runlevel_groupings.txt" 2>/dev/null | sed 's/^/  /'
echo "-- 13_ghi_permanova --"
grep -E 'GHI|Age|R2|arms|studies' "$GDI/ghi_permanova.txt" 2>/dev/null | head -8 | sed 's/^/  /' || true

echo
echo "=============================================================="
echo "2b. Documentation and outputs must agree"
echo "=============================================================="
# Three consecutive deposits shipped prose that contradicted the outputs sitting
# beside it, because nothing ever checked the two against each other. These are
# the exact claims that were wrong, re-checked mechanically before the build.
GHIDIR="$PDI_ROOT/r_analysis/ghi"
has()   { if grep -qF -- "$2" "$1" 2>/dev/null; then note "ok    $(basename "$1"): $2"
          else bad "$(basename "$1") does not contain: $2"; fi; }
hasnt() { if grep -qF -- "$2" "$1" 2>/dev/null; then bad "$(basename "$1") still contains: $2"
          else note "ok    $(basename "$1"): '$2' absent"; fi; }

# README.md and DATA_CORRECTIONS.md are OWNED by the release builder: lines 41
# and 43 of build_releases.sh copy its canonical copies into the repo on every
# build. Checking the REPO copies here would test two files that the build in
# section 4 is about to overwrite -- which is exactly how a corrected README
# passed this section and then shipped stale anyway. Assert the SOURCES here,
# and verify what actually shipped in section 5.
has   "$BUILDER/README.md"                         "+0.821 SD (k = 36)"
has   "$BUILDER/README.md"                         "+0.868 SD (k = 25)"
hasnt "$BUILDER/README.md"                         "+0.818 SD, k = 36"
has   "$BUILDER/DATA_CORRECTIONS.md"               "now stands at +0.821 SD"
has   "$PDI_ROOT/reference/DEPOSIT_NOTES.md"       "share a measurement compartment"
has   "$VER/61_ghi_scale_and_sensitivity.txt"      "challenge +0.8461 (k=25"
# A script edit is not a deposit edit. 61's closing sentence was corrected in
# the script during the documentation pass, but the .txt beside it is only
# rewritten when 61 runs -- so the deposit kept claiming that only the additive
# magnitude moves under the immunoglobulin re-orientation, three lines under a
# table showing the challenge magnitude moving from +0.8461 to +0.6988. Assert
# the OUTPUT, not the source.
has   "$VER/61_ghi_scale_and_sensitivity.txt"      "Both magnitudes move"
hasnt "$VER/61_ghi_scale_and_sensitivity.txt"      "only the magnitude of the"
has   "$VER/63_rob_and_publication_bias.txt"       "additive +0.820 (k=36"
has   "$VER/63_rob_and_publication_bias.txt"       "k counts independent experiments"
hasnt "$VER/51_batch_correction_check.txt"         "~43%"
hasnt "$GHIDIR/ghi_diagnostics.txt"                "0.653"
has   "$GDI/gdi_construct_summary.txt"             "NO PANEL COULD BE FORMED"

n=$(grep -c 'Phenotype_Comparison' "$VER/ghi_scores_0_10.csv" 2>/dev/null || true)
[[ "$n" -eq 0 ]] && note "ok    ghi_scores_0_10.csv: no Phenotype_Comparison rows" \
                 || bad "ghi_scores_0_10.csv carries $n Phenotype_Comparison row(s)"
# Count the FIELD, not the string. The file is quoted, and "Control" is also a
# Treatment_Group value, so the pattern ,"Control", occurs 178 times across the
# 150 Control rows -- 28 of them carry it twice. grep -c counts lines and so
# happens to give 150, but the assertion means "rows whose Challenge_Type is
# Control", and that is what this counts.
n=$(VER="$VER" Rscript -e 'd <- read.csv(file.path(Sys.getenv("VER"), "ghi_scores_0_10.csv")); cat(sum(d$Challenge_Type == "Control"))' 2>/dev/null)
note "control rows in the 0-10 scale: $n (expect 150)"
[[ "$n" -eq 150 ]] || bad "0-10 scale is anchored on $n controls, not 150"
n=$(grep -rl 'scratch/user' "$PDI_ROOT/r_analysis" 2>/dev/null | wc -l)
[[ "$n" -eq 0 ]] && note "ok    no absolute scratch paths in r_analysis outputs" \
                 || bad "$n r_analysis output(s) contain an absolute scratch path"

echo
echo "=============================================================="
echo "3. Quarantine scratch matrices (keep them, keep them out of the release)"
echo "=============================================================="
mkdir -p "$BACKUP"
for f in genus_clr_PRE_RERUN.tsv genus_clr_RUN_A.tsv; do
  if [[ -f "$GEN/$f" ]]; then
    mv -v "$GEN/$f" "$BACKUP/$f" | sed 's/^/  /'
  else
    note "$f not present (already moved)"
  fi
done
note "PRE_RERUN is the PUBLISHED pre-correction matrix - preserved in $BACKUP, do not delete"

# The patch scripts write timestamped copies beside the file they edit, under
# two different conventions: the R-script patches write *_backup_*, the later
# Python patches write *.bak_<stamp>. Both are working files, not deliverables,
# and the release builder sweeps wholesale -- so move them out before building
# rather than discovering them in the tarball afterwards. The search covers the
# whole of PDI_ROOT to depth 3, not just scripts/, because the documentation
# patch writes backups beside README.md and inside reference/. Moved, never
# deleted, and the backup directory itself is excluded so a re-run is a no-op.
echo
note "sweeping patch-script backups out of the scripts tree"
swept=0
while IFS= read -r -d '' f; do
  mv "$f" "$BACKUP/$(basename "$f")" && { note "moved $(basename "$f")"; swept=$((swept+1)); }
done < <(find "$PDI_ROOT" -maxdepth 3 -type f \( -name '*_backup_*' -o -name '*.bak_*' \) \
              -not -path "$BACKUP/*" -print0 2>/dev/null)
note "swept $swept backup file(s) into $BACKUP"

echo
echo "=============================================================="
echo "4. Build both releases"
echo "=============================================================="
if [[ -d "$BUILDER" ]]; then
  ( cd "$BUILDER" && bash build_releases.sh both ) || bad "build_releases.sh failed"
else
  bad "builder not found at $BUILDER"
fi

echo
echo "=============================================================="
echo "5. Verify the built analysis tarball"
echo "=============================================================="
TB=$(ls -t "$RELDIR"/pdi_analysis_*.tar.gz 2>/dev/null | head -1)
if [[ -z "${TB:-}" || ! -s "$TB" ]]; then
  bad "no analysis tarball found in $RELDIR"
else
  # Report the exact byte count from stat, not du: du reports disk usage and, on
  # a Lustre scratch filesystem, can read low while writes are still flushing.
  # An earlier run of this script reported 16M for a file whose own members sum
  # to over 100MB, which was a measurement artefact, not a bad archive.
  BYTES=$(stat -c %s "$TB" 2>/dev/null || wc -c < "$TB")
  note "tarball: $TB"
  note "size: $BYTES bytes ($(awk -v b="$BYTES" 'BEGIN{printf "%.1f MB", b/1048576}')), written $(date -r "$TB" '+%H:%M:%S')"

  # gzip -t walks the entire compressed stream, so it catches a truncated or
  # still-being-written archive that `tar tzf` might list happily.
  if gzip -t "$TB" 2>/dev/null; then note "gzip integrity OK (full stream readable)"
  else bad "gzip integrity FAILED - archive is truncated or corrupt"; fi

  LIST=/tmp/pdi_tb_list.txt
  tar tzvf "$TB" > /tmp/pdi_tb_long.txt
  awk '{print $NF}' /tmp/pdi_tb_long.txt > "$LIST"
  note "entries: $(wc -l < "$LIST")"
  SUM=$(awk '$3 ~ /^[0-9]+$/ {s+=$3} END{printf "%.0f", s}' /tmp/pdi_tb_long.txt)
  note "members sum to $(awk -v b="$SUM" 'BEGIN{printf "%.1f MB", b/1048576}') uncompressed"

  # The three PICRUSt2 tables are the bulk of the deposit and are already
  # gzipped, so they cannot shrink. Check BYTES, not entry count: counting
  # entries passes on an empty file.
  for p in KO_metagenome_out/pred_metagenome_unstrat EC_metagenome_out/pred_metagenome_unstrat pathways_out/path_abun_unstrat; do
    sz=$(awk -v pat="$p" '$NF ~ pat {print $3}' /tmp/pdi_tb_long.txt | head -1)
    if [[ -z "${sz:-}" ]]; then bad "picrust2 $p MISSING from release"
    elif [[ "$sz" -lt 1000000 ]]; then bad "picrust2 $p present but only $sz bytes (expect >1MB)"
    else note "picrust2 $(basename "$p"): $(awk -v b="$sz" 'BEGIN{printf "%.1f MB", b/1048576}')"; fi
  done

  n=$(grep -cE '_backup_|\.bak_' "$LIST" || true);           [[ "$n" -eq 0 ]] && note "no backup files (0)" || bad "$n backup entries present"
  n=$(grep -cE 'genus_clr_(PRE_RERUN|RUN_A)' "$LIST" || true); [[ "$n" -eq 0 ]] && note "no scratch matrices (0)" || bad "$n scratch matrix entries present"
  n=$(grep -c 'GDI_Metadata.*\.xlsx' "$LIST" || true);       [[ "$n" -eq 2 ]] && note "workbook lines: 2 (current + superseded)" || bad "expected 2 workbook lines, found $n"
  n=$(grep -cv '^pdi/' "$LIST" || true);                     [[ "$n" -eq 0 ]] && note "all paths start at pdi/" || bad "$n entries outside pdi/"
  n=$(grep -c 'genus_clr_batch_adjusted.tsv' "$LIST" || true); [[ "$n" -ge 1 ]] && note "final CLR matrix included" || bad "CLR matrix missing from release"

  # the deposited outputs must be the NEW ones, not a stale copy
  REF=$(tar xzOf "$TB" pdi/r_analysis/gdi/runlevel_groupings.txt 2>/dev/null \
        | grep -A1 'Reference - control vs challenge' | tail -1 | tr -s ' ')
  note "deposited reference line:$REF"
  echo "$REF" | grep -q '315 runs, 14 studies' \
    && note "deposited 28 output is current" \
    || bad "deposited runlevel_groupings.txt is STALE (expected 315 runs, 14 studies)"

  # The definitive documentation check: not the repo copy, which the build
  # overwrites, but the bytes actually inside the archive being uploaded.
  for pat in "+0.821 SD (k = 36)" "+0.868 SD (k = 25)"; do
    if tar xzOf "$TB" pdi/README.md 2>/dev/null | grep -qF -- "$pat"
    then note "shipped README carries: $pat"
    else bad "shipped README does NOT carry: $pat"; fi
  done
  if tar xzOf "$TB" pdi/README.md 2>/dev/null | grep -qF -- "+0.818 SD, k = 36"
  then bad "shipped README still carries the superseded +0.818 SD, k = 36"
  else note "shipped README: superseded effect size absent"; fi
  if tar xzOf "$TB" pdi/DATA_CORRECTIONS.md 2>/dev/null | grep -qF -- "now stands at +0.821 SD"
  then note "shipped DATA_CORRECTIONS.md is current"
  else bad "shipped DATA_CORRECTIONS.md is STALE (builder copy not patched)"; fi
  if tar xzOf "$TB" pdi/reference/PATH_SANITISATION.txt 2>/dev/null | grep -qF -- "PATH SANITISATION"
  then note "path-sanitisation declaration present in the deposit"
  else bad "reference/PATH_SANITISATION.txt missing from the deposit"; fi
fi

echo
echo "=============================================================="
if [[ "$FAIL" -eq 0 ]]; then
  echo "PASS - archive is consistent and safe to upload."
  echo "Still outstanding before the Zenodo DOI is final: amendment A9 needs the"
  echo "Web of Science and ProQuest record counts. The GitHub push can go first."
else
  echo "FAIL - see the FAIL lines above. Do not upload yet."
fi
echo "=============================================================="
exit "$FAIL"
