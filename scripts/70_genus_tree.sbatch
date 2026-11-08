#!/bin/bash
#SBATCH --job-name=pdi_genus_tree
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --output=/scratch/user/%u/pdi/logs/genus_tree_%j.log
set -euo pipefail
PDI_ROOT="${PDI_ROOT:-/scratch/user/$USER/pdi}"
QIIME2_ENV="/scratch/user/$USER/.conda/envs/qiime2"
export PATH="${QIIME2_ENV}/bin:$PATH"
NCORES="${SLURM_CPUS_PER_TASK:-8}"
EXP="${PDI_ROOT}/merged/final_exported"
OUT="${PDI_ROOT}/r_analysis/tree"
mkdir -p "$OUT" "${PDI_ROOT}/logs"
for f in "${EXP}/feature-table_w_tax.tsv" "${EXP}/dna-sequences.fasta" \
         "${PDI_ROOT}/r_analysis/genus/genus_counts.tsv"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
echo "[1/3] Choosing one representative sequence per genus"
python3 - "$EXP" "$PDI_ROOT" "$OUT" <<'PY'
import sys, re
import pandas as pd
exp, root, out = sys.argv[1:4]
genera = pd.read_csv(f"{root}/r_analysis/genus/genus_counts.tsv", sep="\t", usecols=[0]).iloc[:, 0]
genera = [g for g in genera if g.lower() != "g__unclassified"]
seqs, name = {}, None
for line in open(f"{exp}/dna-sequences.fasta"):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]; seqs[name] = []
    elif name:
        seqs[name].append(line)
seqs = {k: "".join(v) for k, v in seqs.items()}
rows = []
for chunk in pd.read_csv(f"{exp}/feature-table_w_tax.tsv", sep="\t", skiprows=1,
                         chunksize=2000, low_memory=False):
    idc, taxc = chunk.columns[0], chunk.columns[-1]
    counts = chunk.drop(columns=[idc, taxc]).apply(pd.to_numeric, errors="coerce").fillna(0)
    rows.append(pd.DataFrame({"asv": chunk[idc].astype(str),
                              "taxonomy": chunk[taxc].astype(str),
                              "total": counts.sum(axis=1).values}))
t = pd.concat(rows, ignore_index=True)
def rank(tax, p):
    m = re.search(rf"{p}__([^;]*)", tax)
    return m.group(1).strip() if m and m.group(1).strip() else ""
t["genus"] = "g__" + t.taxonomy.map(lambda x: rank(x, "g"))
t["length"] = t.asv.map(lambda a: len(seqs.get(a, "")))
tips, missing = [], []
for g in genera:
    sub = t[(t.genus == g) & (t.length > 0)]
    if sub.empty:
        missing.append(g); continue
    pref = sub[(sub.length >= 380) & (sub.length <= 480)]
    pick = (pref if not pref.empty else sub).sort_values("total", ascending=False).iloc[0]
    tax = pick.taxonomy
    tips.append(dict(tip=re.sub(r"[^A-Za-z0-9_.-]", "_", g[3:]), genus=g[3:],
                     phylum=rank(tax, "p"), cls=rank(tax, "c"), order=rank(tax, "o"),
                     family=rank(tax, "f"), rep_asv=pick.asv, length=int(pick.length),
                     v3v4_length=bool(380 <= pick.length <= 480),
                     rep_count=int(pick.total), genus_count=int(sub.total.sum()),
                     n_asvs=int(len(sub))))
tips = pd.DataFrame(tips)
tips.to_csv(f"{out}/genus_tree_tips.tsv", sep="\t", index=False)
with open(f"{out}/genus_rep_seqs.fasta", "w") as fh:
    for _, r in tips.iterrows():
        fh.write(f">{r.tip}\n{seqs[r.rep_asv]}\n")
print(f"genera requested: {len(genera)}  tips written: {len(tips)}  "
      f"V3-V4 length: {int(tips.v3v4_length.sum())}  not found: {missing}")
PY
echo "[2/3] Aligning with MAFFT L-INS-i on ${NCORES} threads"
mafft --localpair --maxiterate 1000 --thread "$NCORES" --quiet \
  "${OUT}/genus_rep_seqs.fasta" > "${OUT}/genus_rep_seqs.aln.fasta"
echo "[3/3] Maximum-likelihood tree with IQ-TREE"
IQ=$(command -v iqtree2 || command -v iqtree || true)
if [[ -n "$IQ" ]]; then
  "$IQ" -s "${OUT}/genus_rep_seqs.aln.fasta" -m MFP -B 1000 \
        -T AUTO -ntmax "$NCORES" --seed 1 -redo --prefix "${OUT}/genus_tree" --quiet
else
  echo "IQ-TREE not found; falling back to FastTree GTR+CAT"
  FastTree -nt -gtr "${OUT}/genus_rep_seqs.aln.fasta" > "${OUT}/genus_tree.treefile"
fi
echo "Done. Files in ${OUT}:"
ls -lh "$OUT"
grep -m1 "Best-fit model" "${OUT}/genus_tree.iqtree" 2>/dev/null || true
