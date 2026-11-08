#!/bin/bash
#SBATCH --job-name=pdi_tree
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=180G
#SBATCH --output=/scratch/user/%u/pdi/logs/tree_%j.log
# =============================================================================
# 03b_build_phylogeny.sh
#
# Builds the phylogenetic tree from the merged representative sequences.
#
# Usage (after 03_merge_and_taxonomy.sh has completed):
#   sbatch 03b_build_phylogeny.sh
#
# WHY THIS IS A SEPARATE SCRIPT: the tree is required only for phylogeny-aware
# diversity metrics (weighted and unweighted UniFrac). Nothing else in this
# pipeline consumes it - genus-level abundance work (the GDI and GHI
# construction) never touches it, and stage 04's rarefaction step checks for
# the tree's existence and proceeds without it if absent.
#
# Keeping it separate means the main pipeline (stages 03 and 04) finishes in a
# few hours rather than ~20, and this job can be run later, unattended, without
# blocking anything downstream. It reads only merged-rep-seqs.qza, which does
# not change once stage 03 has completed - so there is no penalty for deferring
# it, and no need to re-run anything else when it is eventually built.
#
# RUNTIME: this is by far the slowest step in the pipeline, commonly 10-20
# hours. MAFFT must align sequences of very different lengths in a single pass
# (this corpus spans ~253 bp V4, ~460 bp V3-V4, ~570 bp V3-V5, and ~1500 bp
# PacBio full-length reads). Alignment cost scales with length *variance*, not
# just sequence count, and this does not parallelise well - allocating more
# cores has limited effect. It is an algorithmic bottleneck, not a resource one.
#
# NO LOG OUTPUT DURING THE RUN: MAFFT and FastTree write nothing to this log
# while working. Hours of silence is expected and is not a hang. To confirm the
# job is genuinely progressing:
#     squeue -u $USER -o "%N" -j <jobid>
#     ssh <node> "ps aux | grep -E 'mafft|FastTree'"
# A climbing CPU-time column means it is working.
#
# INTERPRETIVE CAVEAT: aligning V3-V4, V4, V3-V5 and full-length reads in one
# alignment is imperfect - short-region ASVs overlap only part of the
# full-length sequences, so branch lengths (and therefore UniFrac distances)
# should be read with that heterogeneity in mind. Bray-Curtis, which uses no
# phylogeny, serves as a cross-check; agreement between the two suggests
# conclusions are not artefacts of the alignment.
# =============================================================================

# Under sbatch, SLURM copies this script to a spool directory and runs it from
# there, so $(dirname "$0") does not point at the scripts directory.
# SLURM_SUBMIT_DIR is the directory sbatch was invoked from, which is where the
# repository actually lives; this falls back to the script's own directory for
# interactive use on the login node.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
if [[ ! -f "${SCRIPT_DIR}/00_config.sh" ]]; then
  echo "FATAL: 00_config.sh not found in ${SCRIPT_DIR}" >&2
  echo "Submit with:  cd /scratch/user/$USER/pdi/scripts && sbatch $(basename "$0")" >&2
  exit 1
fi
source "${SCRIPT_DIR}/00_config.sh"

if [[ -z "${PDI_ROOT:-}" ]] || ! command -v log >/dev/null 2>&1; then
  echo "FATAL: 00_config.sh loaded but PDI_ROOT/log() are unset." >&2
  exit 1
fi
activate_qiime

# ---- Input check --------------------------------------------------------------
# Fails early with a clear instruction rather than partway into a long job.
if [[ ! -s "${MERGED_DIR}/merged-rep-seqs.qza" ]]; then
  log "ERROR: ${MERGED_DIR}/merged-rep-seqs.qza not found."
  log "       Run 03_merge_and_taxonomy.sh first."
  exit 1
fi

if already_done "${MERGED_DIR}/rooted-tree.qza"; then
  log "Tree already built. Nothing to do."
  exit 0
fi

# =============================================================================
# BUILD THE TREE
#
# align-to-tree-mafft-fasttree runs four steps in sequence:
#   1. MAFFT      - multiple sequence alignment of all ASVs
#   2. mask       - drops alignment columns that are mostly gaps or highly
#                   variable, which would otherwise add noise rather than signal
#   3. FastTree   - infers an unrooted tree by approximate maximum likelihood
#   4. midpoint   - roots the tree at its midpoint, since 16S data provides no
#                   outgroup to root on
# =============================================================================
log "Building phylogeny (MAFFT + FastTree) from merged-rep-seqs.qza"
log "  This is expected to take 10-20 hours with no intermediate log output."

$QIIME phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${MERGED_DIR}/merged-rep-seqs.qza" \
  --p-n-threads "$NCORES" \
  --o-alignment "${MERGED_DIR}/aligned-rep-seqs.qza" \
  --o-masked-alignment "${MERGED_DIR}/masked-aligned-rep-seqs.qza" \
  --o-tree "${MERGED_DIR}/unrooted-tree.qza" \
  --o-rooted-tree "${MERGED_DIR}/rooted-tree.qza"

log "=== PHYLOGENY COMPLETE ==="
log "Outputs:"
log "  ${MERGED_DIR}/aligned-rep-seqs.qza"
log "  ${MERGED_DIR}/masked-aligned-rep-seqs.qza"
log "  ${MERGED_DIR}/unrooted-tree.qza"
log "  ${MERGED_DIR}/rooted-tree.qza   <- this is the one downstream steps use"
log ""
log "UniFrac-based analyses (and stage 04's rarefaction curves, if re-run) can"
log "now use the tree. Nothing else needs re-running because of this."
