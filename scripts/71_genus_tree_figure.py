#!/usr/bin/env python3
"""71_genus_tree_figure.py

Circular maximum-likelihood tree of the 114 genera retained after the 10%
prevalence filter (one representative V3-V4 ASV per genus; IQ-TREE, model by
BIC, 1,000 ultrafast bootstraps; built by scripts/70_genus_tree.sbatch).

Rings, inside to outside: MaAsLin2 direction of the GHI coefficient in the
primary model (q <= 0.25), then membership of the 14-genus butyrate-producer
guild (GUILD_BUTYRATE in R_scripts/57_scfa_guild.R).

Inputs  r_analysis/tree/genus_tree.treefile, genus_tree_tips.tsv
        r_analysis/gdi/maaslin2_all_results.tsv
Output  figures/genus_tree.pdf, figures/genus_tree.png

Usage   python scripts/71_genus_tree_figure.py [PDI_ROOT]
"""
import os, re, sys, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.lines import Line2D

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expandvars("/scratch/user/$USER/pdi")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import pdi_fig_style as S
S.apply(base=8.0)

TREE = os.path.join(ROOT, "r_analysis/tree/genus_tree.treefile")
TIPS = os.path.join(ROOT, "r_analysis/tree/genus_tree_tips.tsv")
MAAS = os.path.join(ROOT, "r_analysis/gdi/maaslin2_all_results.tsv")
OUT = os.path.join(ROOT, "figures")
Q = 0.25
GUILD_BUTYRATE = ["Faecalibacterium", "Butyricicoccus", "Roseburia", "Anaerostipes",
                  "Blautia", "Coprococcus", "Subdoligranulum", "Eubacterium",
                  "Butyrivibrio", "Anaerobutyricum", "Agathobacter", "Oscillibacter",
                  "Pseudobutyrivibrio", "Shuttleworthia"]

# validly published phylum names (LPSN), matching the article text
PHYLUM_NAME = {"Firmicutes": "Bacillota", "Bacteroidota": "Bacteroidota",
               "Proteobacteria": "Pseudomonadota", "Actinobacteriota": "Actinomycetota"}
PHYLUM_COL = {"Bacillota": "#2E5A7D", "Bacteroidota": "#B07C3C",
              "Pseudomonadota": "#6E5A8C", "Actinomycetota": "#4F7A5E",
              "Other phyla": "#9AA6B1"}
POS, NEG, NS = S.DIVERGING[-1], S.DIVERGING[0], "#E3E8EC"
GUILD_COL = S.SAGE
BAR_COL = "#8A99A8"
COUNTS = os.path.join(ROOT, "r_analysis/genus/genus_counts.tsv")


# ---------------------------------------------------------------- newick
class Node:
    __slots__ = ("name", "length", "support", "children", "parent")

    def __init__(self):
        self.name, self.length, self.support = None, 0.0, None
        self.children, self.parent = [], None


def parse_newick(s):
    s = s.strip().rstrip(";")
    root, cur, i = Node(), None, 0
    cur = root
    tok = re.compile(r"[^(),:;]+")
    while i < len(s):
        c = s[i]
        if c == "(":
            n = Node(); n.parent = cur; cur.children.append(n); cur = n; i += 1
        elif c == ",":
            n = Node(); n.parent = cur.parent; cur.parent.children.append(n); cur = n; i += 1
        elif c == ")":
            cur = cur.parent; i += 1
            m = tok.match(s, i)
            if m:
                lab = m.group(0); i = m.end()
                try:
                    cur.support = float(lab.split("/")[-1])
                except ValueError:
                    cur.name = lab
        elif c == ":":
            m = re.compile(r"[-0-9.eE]+").match(s, i + 1)
            cur.length = float(m.group(0)); i = m.end()
        else:
            m = tok.match(s, i); cur.name = m.group(0); i = m.end()
    return root


def leaves(n):
    return [n] if not n.children else [l for c in n.children for l in leaves(c)]


def reroot_midpoint(root):
    """Return a new root placed at the midpoint of the longest tip-to-tip path."""
    # undirected adjacency
    adj = {}
    def walk(n):
        for c in n.children:
            adj.setdefault(n, []).append((c, c.length)); adj.setdefault(c, []).append((n, c.length)); walk(c)
    walk(root)
    def far(src):
        dist, prev, stack = {src: 0.0}, {src: None}, [src]
        while stack:
            u = stack.pop()
            for v, w in adj.get(u, []):
                if v not in dist:
                    dist[v] = dist[u] + w; prev[v] = u; stack.append(v)
        tipd = {k: d for k, d in dist.items() if not k.children}
        t = max(tipd, key=tipd.get)
        return t, tipd[t], prev
    a, _, _ = far(leaves(root)[0])
    b, L, prev = far(a)
    path = [b]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    # walk from b toward a until half the length
    half, acc = L / 2.0, 0.0
    for u, v in zip(path[:-1], path[1:]):
        w = dict(adj[u])[v]
        if acc + w >= half:
            cut = half - acc
            break
        acc += w
    # support of an edge is held by the child end in the original rooting
    def edge_support(x, y):
        return y.support if y.parent is x else x.support if x.parent is y else None
    sup_uv = edge_support(u, v)
    new = Node()
    def build(src, frm, length, sup):
        n = Node(); n.name, n.length, n.support = src.name, length, sup
        for nb, w in adj[src]:
            if nb is not frm:
                ch = build(nb, src, w, edge_support(src, nb)); ch.parent = n; n.children.append(ch)
        return n
    cu = build(u, v, cut, sup_uv); cv = build(v, u, w - cut, sup_uv)
    for c in (cu, cv):
        c.parent = new; new.children.append(c)
    # collapse unary nodes left where the old root had degree 2
    def squash(n):
        for i, c in enumerate(list(n.children)):
            while len(c.children) == 1:
                g = c.children[0]; g.length += c.length; g.parent = n; n.children[i] = g; c = g
            squash(c)
    squash(new)
    return new


def ladderize(n):
    for c in n.children:
        ladderize(c)
    n.children.sort(key=lambda c: len(leaves(c)))


# ---------------------------------------------------------------- labels
CODE = re.compile(r"\d|group$|^UCG|^RF39$")


def it(w):
    return r"$\mathit{" + w + "}$"


def label(raw):
    t = raw.strip("_")
    m = re.match(r"_?([A-Za-z]+)__(.+)", raw)          # [Eubacterium] hallii group
    if m:
        rest = m.group(2).split("_")
        ep = rest[0] if rest[0] != "group" else ""
        return "[" + it(m.group(1)) + "] " + (it(ep) + " " if ep else "") + "group"
    if t.startswith("Candidatus_"):
        return it("Candidatus") + " " + t.split("_", 1)[1]
    if t == "Clostridium_sensu_stricto_1":
        return it("Clostridium") + " sensu stricto 1"
    if t == "Escherichia-Shigella":
        return it("Escherichia") + "-" + it("Shigella")
    if t in ("uncultured", "Incertae_Sedis"):
        return t.replace("_", " ")
    parts = t.split("_")
    if len(parts) > 1 and re.match(r"^[A-Z][a-z]+(aceae|ia|ales)$", parts[0]):
        return it(parts[0]) + " " + " ".join(parts[1:])   # Lachnospiraceae UCG-010
    if parts[0] == "Family":
        return " ".join(parts)
    if CODE.search(t) or t.isupper():
        return t.replace("_", " ")
    return it(t)


def norm(s):
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]", "_", re.sub(r"^g__", "", s))).strip("_")


# ---------------------------------------------------------------- data
tree = parse_newick(open(TREE).read())
tree = reroot_midpoint(tree)
ladderize(tree)
tips = pd.read_csv(TIPS, sep="\t").set_index("tip")
cnt = pd.read_csv(COUNTS, sep="\t", index_col=0)
abund = (cnt / cnt.sum(0)).mean(1) * 100          # mean relative abundance, %
abund.index = abund.index.map(norm)
abund = abund.to_dict()
ma = pd.read_csv(MAAS, sep="\t")
ma = ma[ma["metadata"] == "GHI"].copy()
ma["k"] = ma["feature"].map(norm)
mres = ma.set_index("k")[["coef", "qval"]]

L = leaves(tree)
n = len(L)
assert n == len(tips) == 114, (n, len(tips))


def depth(nd):
    d = 0.0
    while nd.parent is not None:
        d += nd.length; nd = nd.parent
    return d


maxd = max(depth(l) for l in L)
GAP = 0.06                                   # fraction of the circle left open
ang = {}
for i, l in enumerate(L):
    ang[l] = 2 * math.pi * (GAP / 2 + (1 - GAP) * (i + 0.5) / n) + math.pi / 2


def set_ang(nd):
    if nd.children:
        for c in nd.children:
            set_ang(c)
        ang[nd] = (ang[nd.children[0]] + ang[nd.children[-1]]) / 2
set_ang(tree)

R_TREE = 1.0
rad = lambda nd: R_TREE * depth(nd) / maxd

fig, ax = plt.subplots(figsize=(S.W2, S.W2))
ax.set_aspect("equal"); ax.axis("off")
LW = 0.55
EDGE = "#3B4A59"


def draw(nd):
    r0 = rad(nd)
    if nd.children:
        a = [ang[c] for c in nd.children]
        th = np.linspace(min(a), max(a), 60)
        ax.plot(r0 * np.cos(th), r0 * np.sin(th), color=EDGE, lw=LW, solid_capstyle="butt")
        for c in nd.children:
            r1 = rad(c)
            ax.plot([r0 * math.cos(ang[c]), r1 * math.cos(ang[c])],
                    [r0 * math.sin(ang[c]), r1 * math.sin(ang[c])], color=EDGE, lw=LW)
            draw(c)
        if nd.support is not None and nd.support >= 95 and nd.parent is not None:
            ax.plot(r0 * math.cos(ang[nd]), r0 * math.sin(ang[nd]), "o", ms=2.0,
                    color=S.INK, mec="none", zorder=4)
draw(tree)

# dotted leaders from each tip to the rings
R_RING0, R_RING1, R_RING2, W_RING = 1.02, 1.08, 1.14, 0.05
R_BAR, BAR_LEN, AB_MIN, AB_MAX = 1.20, 0.17, 3e-3, 10.0   # bars on log10 % scale
for l in L:
    a = ang[l]
    ax.plot([rad(l) * math.cos(a), (R_RING0 - 0.005) * math.cos(a)],
            [rad(l) * math.sin(a), (R_RING0 - 0.005) * math.sin(a)],
            color=S.GREY_L, lw=0.35, ls=(0, (1, 1.2)))

step = 360 * (1 - GAP) / n
phyl_of = {}
n_pos = n_neg = n_ns = n_nt = n_guild = 0
for l in L:
    g = tips.loc[l.name]
    ph = PHYLUM_NAME.get(g["phylum"], "Other phyla")
    phyl_of[l] = ph
    a = math.degrees(ang[l])
    ax.add_patch(Wedge((0, 0), R_RING0 + W_RING, a - step / 2, a + step / 2,
                       width=W_RING, facecolor=PHYLUM_COL[ph], edgecolor="none"))
    key = norm(g["genus"])
    if key in mres.index:
        c, q = mres.loc[key, "coef"], mres.loc[key, "qval"]
        col = (POS if c > 0 else NEG) if q <= Q else NS
        n_pos += (q <= Q) and c > 0; n_neg += (q <= Q) and c < 0; n_ns += q > Q
        ax.add_patch(Wedge((0, 0), R_RING1 + W_RING, a - step / 2 + 0.15, a + step / 2 - 0.15,
                           width=W_RING, facecolor=col, edgecolor="none"))
    else:
        n_nt += 1
    if g["genus"] in GUILD_BUTYRATE:
        n_guild += 1
        ax.add_patch(Wedge((0, 0), R_RING2 + W_RING, a - step / 2 + 0.15, a + step / 2 - 0.15,
                           width=W_RING, facecolor=GUILD_COL, edgecolor="none"))
    # tip label, coloured by phylum
    ab = max(abund.get(key, AB_MIN), AB_MIN)
    blen = BAR_LEN * math.log10(ab / AB_MIN) / math.log10(AB_MAX / AB_MIN)
    ax.add_patch(Wedge((0, 0), R_BAR + blen, a - step * 0.26, a + step * 0.26,
                       width=blen, facecolor=BAR_COL, edgecolor="none"))
    r_lab = R_BAR + BAR_LEN + 0.02
    deg = a % 360
    flip = 90 < deg < 270
    ax.text(r_lab * math.cos(ang[l]), r_lab * math.sin(ang[l]), label(l.name),
            rotation=deg + 180 if flip else deg, rotation_mode="anchor",
            ha="right" if flip else "left", va="center", fontsize=5.6,
            color=S.INK)


# abundance key in the open gap at the top
for v_, lab_ in ((0.01, "0.01"), (0.1, "0.1"), (1, "1"), (10, "10%")):
    rr_ = R_BAR + BAR_LEN * math.log10(v_ / AB_MIN) / math.log10(AB_MAX / AB_MIN)
    ax.plot([-0.012, 0.012], [rr_, rr_], color=S.INK, lw=0.5)
    ax.text(0.025, rr_, lab_, ha="left", va="center", fontsize=5.4, color=S.INK)
ax.plot([0, 0], [R_BAR, R_BAR + BAR_LEN], color=S.INK, lw=0.5)
ax.text(0, R_BAR + BAR_LEN + 0.03, "Mean relative\nabundance", ha="center", va="bottom",
        fontsize=5.6, color=S.INK, linespacing=1.1)

# scale bar, lower left
sb = 0.1
y0 = -1.80
x0 = -1.88
ax.plot([x0, x0 + R_TREE * sb / maxd], [y0, y0], color=S.INK, lw=0.8)
ax.text(x0, y0 + 0.04, f"{sb:g} substitutions per site",
        ha="left", va="bottom", fontsize=6.4, color=S.INK)

handles = [Line2D([], [], ls="none", marker="s", ms=6, color=PHYLUM_COL[p], label=p)
           for p in ("Bacillota", "Bacteroidota", "Pseudomonadota", "Actinomycetota", "Other phyla")]
leg1 = ax.legend(handles=handles, title="Phylum (inner ring)", loc="upper left", bbox_to_anchor=(0.0, 1.0),
                 fontsize=6.8, title_fontsize=7.2, handletextpad=0.3, labelspacing=0.35,
                 alignment="left")
ax.add_artist(leg1)
h2 = [Line2D([], [], ls="none", marker="s", ms=6, color=POS, label=f"Higher with GHI, $q$ ≤ {Q:g}"),
      Line2D([], [], ls="none", marker="s", ms=6, color=NEG, label=f"Lower with GHI, $q$ ≤ {Q:g}"),
      Line2D([], [], ls="none", marker="s", ms=6, color=NS, label=f"Tested, $q$ > {Q:g}"),
      Line2D([], [], ls="none", marker="s", ms=6, color=GUILD_COL, label="Butyrate guild"),
      Line2D([], [], ls="none", marker="s", ms=6, color=BAR_COL, label="Mean relative abundance (bars)"),
      Line2D([], [], ls="none", marker="o", ms=3, color=S.INK, label="Bootstrap ≥ 95")]
ax.legend(handles=h2, title="Rings and nodes", loc="upper right", bbox_to_anchor=(1.0, 1.0),
          fontsize=6.8, title_fontsize=7.2, handletextpad=0.3, labelspacing=0.35, alignment="left")

lim = 1.97
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
fig.subplots_adjust(0, 0, 1, 1)
os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, "genus_tree.pdf"))
fig.savefig(os.path.join(OUT, "genus_tree.png"), dpi=600)
# phylogenetic signal of the GHI coefficient: Moran's I, inverse patristic weights
def _anc(nd):
    out, d = {}, 0.0
    while nd is not None:
        out[id(nd)] = d; d += nd.length; nd = nd.parent
    return out
tested = [l for l in L if norm(tips.loc[l.name, "genus"]) in mres.index]
x = np.array([mres.loc[norm(tips.loc[l.name, "genus"]), "coef"] for l in tested])
A = [_anc(l) for l in tested]
D = np.array([[min(a[k] + b[k] for k in a if k in b) for b in A] for a in A])
Wm = np.zeros_like(D); Wm[D > 0] = 1 / D[D > 0]
def moran(v):
    z = v - v.mean(); return len(v) / Wm.sum() * (z @ Wm @ z) / (z @ z)
rng = np.random.default_rng(1)
I0 = moran(x); perm = np.array([moran(rng.permutation(x)) for _ in range(9999)])
p_up = (1 + (perm >= I0).sum()) / (len(perm) + 1)
with open(os.path.join(ROOT, "r_analysis/tree/phylo_signal_moran.txt"), "w") as fh:
    fh.write(f"Moran's I of the MaAsLin2 GHI coefficient on the genus tree\n"
             f"weights 1/patristic distance (midpoint-rooted IQ-TREE ML tree)\n"
             f"genera {len(x)}  I = {I0:.4f}  E[I] = {-1/(len(x)-1):.4f}  "
             f"P (one-sided, 9,999 permutations, seed 1) = {p_up:.4f}\n")
print(f"Moran I {I0:.4f} (E {-1/(len(x)-1):.4f}) P {p_up:.4f}")
print(f"tips {n} | tested {n - n_nt} (higher {n_pos}, lower {n_neg}, q>{Q:g} {n_ns}) | "
      f"not tested {n_nt} | guild {n_guild}")
