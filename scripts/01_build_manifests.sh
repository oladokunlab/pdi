#!/bin/bash
# =============================================================================
# 01_build_manifests.sh
#
# Builds one QIIME 2 manifest file per study from studies.tsv and the FASTQ
# files present under raw_data/<Study_ID>/. Manifests are used instead of
# Casava-style import because sequencing data retrieved from SRA/ENA/NGDC/NMDC
# does not follow Illumina's Casava filename convention.
#
# Usage (login node - no internet access required):
#   bash 01_build_manifests.sh 2>&1 | tee logs/01_manifests.log
#
# For each study, one of four manifest formats is produced, selected by the
# `layout` field in studies.tsv:
#   paired         - standard paired-end, forward/reverse manifest
#   forward_only   - paired data on disk; only the forward read is used
#   single_merged  - deposited as one pre-joined read per sample
#   ccs            - PacBio long-read, single file per sample
#
# A summary log records OK / INCOMPLETE / MISSING_DIR / HOLD for every study,
# so incomplete data can be identified before submitting stage 02.
# =============================================================================

source "$(dirname "$0")/00_config.sh"

SUMMARY="${MANIFEST_DIR}/manifest_build_summary.log"
: > "$SUMMARY"

log "Building manifests from ${STUDY_TABLE}"

# Carriage returns are stripped from studies.tsv on read. The file is edited
# both in Excel and by Python's csv writer, either of which can introduce
# CRLF line endings; a trailing \r attached to the last field would silently
# break equality tests such as [[ "$read_orientation" == "swapped" ]].
tail -n +2 "$STUDY_TABLE" | tr -d '\r' | while IFS=$'\t' read -r study n_runs region platform layout pf pr status primers_present read_orientation; do

  [[ -z "$study" ]] && continue

  sdir="${RAW_DIR}/${study}"
  man="${MANIFEST_DIR}/${study}_manifest.tsv"

  if [[ ! -d "$sdir" ]]; then
    echo -e "${study}\tMISSING_DIR\t0/${n_runs}" >> "$SUMMARY"
    continue
  fi

  # ---- paired-end studies --------------------------------------------------
  # Covers standard Illumina paired-end and DNBSEQ (MGI/BGI) data, which uses
  # the same paired FASTQ format.
  #
  # read_orientation=swapped: for studies where the reverse primer sequence
  # is found at the start of the file labelled _1, the forward and reverse
  # file assignment is corrected here by writing _2 into the forward manifest
  # column and _1 into the reverse column. This is done at the manifest
  # level rather than by swapping the primer sequences passed to cutadapt,
  # because DADA2 emits representative sequences in the orientation of
  # whichever file is treated as forward. Correcting file assignment keeps
  # every study's output sequences in one consistent orientation, which is
  # required for feature-table merging across studies (feature identity is
  # determined by an md5 hash of the sequence).
  if [[ "$layout" == "paired" ]]; then
    printf 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n' > "$man"
    count=0
    for r1 in "$sdir"/*_1.fastq.gz "$sdir"/*_1.fastq "$sdir"/*_R1*.fastq.gz "$sdir"/*_R1*.fastq; do
      [[ -e "$r1" ]] || continue
      r2="${r1/_1./_2.}"; r2="${r2/_R1/_R2}"
      [[ -e "$r2" ]] || { log "WARN ${study}: no mate for $(basename "$r1")"; continue; }
      sid=$(basename "$r1"); sid="${sid%%_1.*}"; sid="${sid%%_R1*}"
      # Strips a ".sra" filename infix left by fastq-dump (e.g.
      # SRR123.sra_1.fastq.gz), so the sample-id matches the plain accession
      # used in the metadata file.
      sid="${sid%.sra}"
      if [[ "$read_orientation" == "swapped" ]]; then
        printf '%s\t%s\t%s\n' "$sid" "$r2" "$r1" >> "$man"
      else
        printf '%s\t%s\t%s\n' "$sid" "$r1" "$r2" >> "$man"
      fi
      count=$((count+1))
    done
    [[ "$read_orientation" == "swapped" ]] && \
      log "  ${study}: R1/R2 swapped in manifest (reverse primer detected at start of _1)"

  # ---- PacBio CCS and pre-merged single-end studies ------------------------
  # Both produce one FASTQ file per sample and share a single-end manifest
  # format. They are distinguished at the denoising stage (denoise-ccs vs.
  # denoise-single) via the layout column, not here.
  #
  # single_merged studies were deposited with forward and reverse reads
  # already joined into one read per sample (approximately 400-470 bp,
  # spanning the full V3-V4 amplicon). The ".sra" filename infix is stripped
  # for the same reason as in the paired-end branch above.
  elif [[ "$layout" == "ccs" || "$layout" == "single_merged" ]]; then
    printf 'sample-id\tabsolute-filepath\n' > "$man"
    count=0
    for f in "$sdir"/*.fastq.gz "$sdir"/*.fastq; do
      [[ -e "$f" ]] || continue
      sid=$(basename "$f"); sid="${sid%%.fastq*}"; sid="${sid%.sra}"
      printf '%s\t%s\n' "$sid" "$f" >> "$man"
      count=$((count+1))
    done

  # ---- forward-only studies -------------------------------------------------
  # Paired R1/R2 files exist on disk, but forward and reverse reads do not
  # overlap sufficiently to be joined (confirmed by direct sequence
  # comparison: the terminal bases of the forward read do not appear in the
  # reverse complement of the reverse read's leading bases). Only the
  # forward read is used; the reverse read is not included in the manifest.
  # This produces a single-end manifest, but is kept as a distinct layout
  # value from single_merged because the source data and processing
  # rationale differ.
  elif [[ "$layout" == "forward_only" ]]; then
    printf 'sample-id\tabsolute-filepath\n' > "$man"
    count=0
    for f in "$sdir"/*_1.fastq.gz "$sdir"/*_1.fastq "$sdir"/*_R1*.fastq.gz "$sdir"/*_R1*.fastq; do
      [[ -e "$f" ]] || continue
      sid=$(basename "$f"); sid="${sid%%_1.*}"; sid="${sid%%_R1*}"; sid="${sid%.sra}"
      printf '%s\t%s\n' "$sid" "$f" >> "$man"
      count=$((count+1))
    done
    log "  ${study}: forward read only (reverse read excluded from manifest)"

  # ---- studies not yet ready for processing --------------------------------
  else
    echo -e "${study}\tHOLD(${status})\t0/${n_runs}" >> "$SUMMARY"
    rm -f "$man"
    continue
  fi

  if [[ "$count" -eq "$n_runs" ]]; then
    echo -e "${study}\tOK\t${count}/${n_runs}" >> "$SUMMARY"
  else
    echo -e "${study}\tINCOMPLETE\t${count}/${n_runs}" >> "$SUMMARY"
  fi
done

log "Done. Summary:"
column -t "$SUMMARY"

echo
log "Anything not OK above must be resolved before stage 02."
log "Expected run total: sum of the n_runs column in studies.tsv."
