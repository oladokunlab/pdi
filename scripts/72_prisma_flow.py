#!/usr/bin/env python3
"""72_prisma_flow.py

PRISMA 2020 flow diagram (Figure 1) in the house figure style.

Counts: identification from the three database exports (PubMed 778, Web of
Science 546, ProQuest 37); deduplication in EndNote; full-text and later
exclusions from GDI_Metadata.xlsx, sheet Excluded_Screened (98 records, primary
reason by the precedence order in the reason mapping) plus one superseded
duplicate report of an included study.

Output  figures/prisma_flow.pdf, figures/prisma_flow.png
Usage   python scripts/72_prisma_flow.py [PDI_ROOT]
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expandvars("/scratch/user/$USER/pdi")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import pdi_fig_style as S
S.apply(base=8.0)

DB = [("PubMed", 778), ("Web of Science", 546), ("ProQuest Dissertations & Theses Global", 37)]
DUP = 159
TA_EXCL = 1053
FT = 149
NOT_RETRIEVED = 0
REASONS = [  # summed over full-text screening, post-extraction, and bioinformatic QC
    ("No extractable treatment-group data", 24),
    ("Deposited reads absent, unretrievable, or unusable", 21),
    ("No gut-health biomarker", 14),
    ("No unchallenged control group", 10),
    ("Data not attributable to a named treatment group", 7),
    ("Not 16S rRNA amplicon data", 6),
    ("Fecal or non-gut sample source", 6),
    ("Biomarker only from non-gut, non-serum tissue", 6),
    ("Fewer than three sequencing runs per group", 4),
    ("Duplicate report of an included study", 1),
]
REPORTS, EXPERIMENTS = 50, 51
STAGE = (57, 37, 4, 1)          # full text, data extraction, sequence QC, duplicate
GHI_ROWS, GHI_SCORED = 566, 529
SEQ = dict(studies=44, denoised=2247, filtered=2064, genera=115, joined=1206, joined_studies=36)

IDENT = sum(n for _, n in DB)
SCREEN = IDENT - DUP
EXCL = sum(n for _, n in REASONS)
assert IDENT == 1361 and SCREEN == 1202 and SCREEN - TA_EXCL == FT
assert FT - NOT_RETRIEVED - EXCL == REPORTS, (FT, EXCL)
assert sum(STAGE) == EXCL


def f(n):
    return f"{n:,}"


W, H, Y0 = S.W2, 7.55, 0.08
fig, ax = plt.subplots(figsize=(W, H - Y0))
ax.set_xlim(0, W); ax.set_ylim(Y0, H); ax.axis("off")
fig.subplots_adjust(0, 0, 1, 1)

EDGE, FILL, LW = S.SPINE, "#F6F8FA", 0.9
FS, FS_H = 7.4, 7.8
LX, LW_BOX = 0.62, 3.05          # left column
RX, RW_BOX = 4.05, 2.92          # right column


def box(x, y, w, h, lines, edge=EDGE, bold_first=True, lw=LW):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.05",
                                facecolor=FILL, edgecolor=edge, lw=lw))
    yy = y + h - 0.13
    for i, (txt, kind) in enumerate(lines):
        ax.text(x + 0.12, yy, txt, ha="left", va="top",
                color=S.GREY if kind == "s" else S.INK,
                fontsize=FS_H if kind == "h" else (FS - 0.4 if kind == "s" else FS),
                fontweight="bold" if kind == "h" else "normal",
                fontstyle="italic" if kind == "s" else "normal")
        yy -= 0.185 if kind == "h" else (0.15 if kind == "s" else 0.165)
        if kind == "s" and i + 1 < len(lines) and lines[i + 1][1] == "t":
            yy -= 0.05
    return (x, y, w, h)


def arrow(x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=S.GREY, lw=0.9,
                                shrinkA=0, shrinkB=0, mutation_scale=8))


def tab(y0, y1, label):
    ax.add_patch(FancyBboxPatch((0.08, y0), 0.30, y1 - y0, boxstyle="round,pad=0,rounding_size=0.05",
                                facecolor=S.BLUE_D, edgecolor="none"))
    ax.text(0.23, (y0 + y1) / 2, label, rotation=90, ha="center", va="center",
            color="white", fontsize=FS_H, fontweight="bold")


ax.text(LX, H - 0.12, "Identification of studies via databases", ha="left", va="top",
        color=S.INK, fontsize=9.0, fontweight="bold")

# identification
y_id = 6.17
b1 = box(LX, y_id + 0.12, LW_BOX, 0.88,
         [(f"Records identified ($n$ = {f(IDENT)})", "h")] +
         [(f"{d} ($n$ = {f(n)})", "t") for d, n in DB])
b1r = box(RX, y_id + 0.32, RW_BOX, 0.50,
          [("Records removed before screening", "h"), (f"Duplicate records ($n$ = {f(DUP)})", "t")])
arrow(LX + LW_BOX, y_id + 0.57, RX, y_id + 0.57)

# screening
y_s1 = 5.43
box(LX, y_s1, LW_BOX, 0.36, [(f"Records screened ($n$ = {f(SCREEN)})", "h")])
box(RX, y_s1, RW_BOX, 0.36, [(f"Records excluded on title and abstract ($n$ = {f(TA_EXCL)})", "h")])
arrow(LX + LW_BOX / 2, y_id + 0.12, LX + LW_BOX / 2, y_s1 + 0.36)
arrow(LX + LW_BOX, y_s1 + 0.18, RX, y_s1 + 0.18)

y_s2 = 4.75
box(LX, y_s2, LW_BOX, 0.36, [(f"Reports sought for retrieval ($n$ = {f(FT)})", "h")])
box(RX, y_s2, RW_BOX, 0.36, [(f"Reports not retrieved ($n$ = {NOT_RETRIEVED})", "h")])
arrow(LX + LW_BOX / 2, y_s1, LX + LW_BOX / 2, y_s2 + 0.36)
arrow(LX + LW_BOX, y_s2 + 0.18, RX, y_s2 + 0.18)

y_s3 = 4.07
h_r = 0.25 + 0.165 * len(REASONS) + 2 * 0.15 + 0.10
y_r = y_s3 + 0.36 - h_r
box(LX, y_s3, LW_BOX, 0.36, [(f"Reports assessed for eligibility ($n$ = {f(FT)})", "h")])
box(RX, y_r, RW_BOX, h_r,
    [(f"Reports excluded ($n$ = {EXCL})", "h"),
     (f"At full-text reading {STAGE[0]}, data extraction {STAGE[1]},", "s"),
     (f"sequence quality control {STAGE[2]}, duplicate {STAGE[3]}", "s")] +
    [(f"{r} ($n$ = {n})", "t") for r, n in REASONS])
arrow(LX + LW_BOX / 2, y_s2, LX + LW_BOX / 2, y_s3 + 0.36)
arrow(LX + LW_BOX, y_s3 + 0.18, RX, y_s3 + 0.18)

# included
y_in = 1.62
box(LX, y_in, LW_BOX, 0.56,
    [(f"Studies included in review ($n$ = {EXPERIMENTS} experiments)", "h"),
     (f"Reports of included studies ($n$ = {REPORTS})", "t")], edge=S.MAROON, lw=1.1)
arrow(LX + LW_BOX / 2, y_s3, LX + LW_BOX / 2, y_in + 0.56)

# analysis tier: what the included studies supplied
y_an = 0.18
h_an = 0.90
box(LX, y_an, LW_BOX, h_an,
    [("Biomarker synthesis (GHI)", "h"),
     (f"{EXPERIMENTS} experiments, {GHI_ROWS} treatment groups extracted", "t"),
     (f"{GHI_SCORED} groups scored in at least one domain", "t"),
     ("18 markers in 7 domains", "t")])
box(RX, y_an, RW_BOX, h_an,
    [("Sequence synthesis (GDI)", "h"),
     (f"{SEQ['studies']} studies, {f(SEQ['denoised'])} runs denoised", "t"),
     (f"{f(SEQ['filtered'])} runs and {SEQ['genera']} genera after filtering", "t"),
     (f"{f(SEQ['joined'])} runs from {SEQ['joined_studies']} studies joined to a GHI value", "t")])
xm = LX + LW_BOX / 2
ax.plot([xm, xm], [y_in, y_an + h_an + 0.17], color=S.GREY, lw=0.9)
ax.plot([xm, RX + RW_BOX / 2], [y_an + h_an + 0.17] * 2, color=S.GREY, lw=0.9)
arrow(xm, y_an + h_an + 0.17, xm, y_an + h_an)
arrow(RX + RW_BOX / 2, y_an + h_an + 0.17, RX + RW_BOX / 2, y_an + h_an)
tab(y_an, y_an + h_an, "Analysis")

tab(y_id + 0.12, y_id + 1.00, "Identification")
tab(y_s3, y_s1 + 0.36, "Screening")
tab(y_in, y_in + 0.56, "Included")

os.makedirs(os.path.join(ROOT, "figures"), exist_ok=True)
fig.savefig(os.path.join(ROOT, "figures", "prisma_flow.pdf"))
fig.savefig(os.path.join(ROOT, "figures", "prisma_flow.png"), dpi=600)
print(f"stages {STAGE} | ghi {GHI_SCORED}/{GHI_ROWS} | seq {SEQ}")
print(f"identified {IDENT} | duplicates {DUP} | screened {SCREEN} | excluded TA {TA_EXCL} | "
      f"full text {FT} | excluded {EXCL} | reports {REPORTS} | experiments {EXPERIMENTS}")
