#!/usr/bin/env python3
"""Draw the PDI figures, one file per figure.

Run from the directory holding pdi_fig_style.py and its fonts/ folder:

    export PDI_ROOT=/scratch/user/$USER/pdi
    python make_figures.py                       # writes into $PDI_ROOT/figures
    python make_figures.py --outdir /tmp/figs    # or somewhere else

Input is r_analysis/ghi/ghi_scores_pca.tsv under PDI_ROOT. If that file does not
already carry Intervention_Class (it does once 04b has been run with the patched
keep list), pass the classification table with --labels.

The meta-analysis estimates in DOMAIN_ESTIMATES are transcribed from the 04e
output and must be updated by hand if 04e is re-run. Every other number in every
figure is computed from the scores file at draw time.
"""
import argparse, os, sys, importlib.util
import numpy as np
import pandas as pd
import scipy.stats as st
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("pdi_fig_style",
                                               os.path.join(HERE, "pdi_fig_style.py"))
S = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(S)
S.apply(base=11.0)

INK, MAROON, CLARET = S.INK, S.MAROON, S.CLARET
GREY, GREY_L, SPINE, SAGE = S.GREY, S.GREY_L, S.SPINE, S.SAGE
ADD_C, TRT_C = S.TEAL, S.BLUE          # additive green, challenged-and-treated blue
COLHEAD = "#5A6B78"                    # column heads, darker than the GREY annotations
GRIDC = "#D6DEE6"                      # gridlines, one step lighter than the leader lines
FUNNEL_S, FUNNEL_A = 26, 0.80
CMAP = LinearSegmentedColormap.from_list("pdi", S.DIVERGING)

CONTROL = "Unchallenged control"
ADDITIVE = "Non-disruptive additive"
TREATED = "Challenged and treated"
CHALLENGE = ["Toxin challenge", "Environmental challenge", "Infectious challenge",
             "Dietary challenge", "Antibiotic challenge"]
NO_CONCURRENT_CONTROL = ["Massacci_2019"]
ANCHOR_MEAN, ANCHOR_SD = 0.222050, 0.902578      # from the 150 published control groups

CLASS_ORDER = ["Toxin challenge", "Environmental challenge", "Infectious challenge",
               "Dietary challenge", CONTROL, "Antibiotic challenge", TREATED, ADDITIVE]

DOMAINS = ["Morphology", "Barrier", "Inflammation", "Anti_Inflammatory", "SCFA",
           "Performance", "Immunoglobulin_Systemic"]
DOMAIN_LABELS = ["Morphology", "Barrier", "Inflammation", "Anti-inflam.", "SCFA",
                 "Performance", "Systemic Ig"]
MIN_PAIRS = 12                                   # below this a matrix cell reads nd

# Domains pooled exactly as 04e does it: per-study paired differences, then a
# DerSimonian-Laird random-effects model weighted by the number of group means.
# Reproduces 04e's Morphology (+0.334 additive, +0.517 challenge), Performance
# (+1.086, +2.390) and Combined (+0.820, +0.819) rows to three decimals.
DOMAIN_COLS = [("Morphology", "GHI_dom_Morphology"),
               ("Barrier", "GHI_dom_Barrier"),
               ("Inflammation", "GHI_dom_Inflammation"),
               ("Anti-inflammatory", "GHI_dom_Anti_Inflammatory"),
               ("SCFA", "GHI_dom_SCFA"),
               ("Performance", "GHI_dom_Performance"),
               ("Systemic Ig", "GHI_dom_Immunoglobulin_Systemic"),
               ("Combined index", "GHI_Combined_raw")]


def pool_re(y, w0):
    y = np.asarray(y, float); w0 = np.asarray(w0, float); k = len(y)
    if k < 3:
        return None
    mu_f = (w0 * y).sum() / w0.sum()
    Q = (w0 * (y - mu_f) ** 2).sum()
    C = w0.sum() - (w0 ** 2).sum() / w0.sum()
    tau2 = max(0.0, (Q - (k - 1)) / C)
    wi = 1 / (1 / w0 + tau2)
    mu = (wi * y).sum() / wi.sum()
    se = np.sqrt(1 / wi.sum())
    return dict(k=k, est=mu, lo=mu - 1.96 * se, hi=mu + 1.96 * se,
                I2=(100 * (Q - (k - 1)) / Q if Q > k - 1 else 0.0))


def study_diffs(g, col, better, worse, skip=()):
    y, w = [], []
    for s, d in g.groupby("Study_ID"):
        if s in skip:
            continue
        d = d[d[col].notna()]
        a = d.loc[d.Class.isin(better), col]
        b = d.loc[d.Class.isin(worse), col]
        if len(a) and len(b):
            y.append(a.mean() - b.mean())
            w.append(len(a) + len(b))
    return y, w


def study_label(s):
    """Lucke_2018 -> Lucke et al. 2018;  Ma_2024_Layers -> Ma et al. 2024 (layers)."""
    import re
    parts = s.split("_")
    yi = next((i for i, p in enumerate(parts)
               if re.fullmatch(r"(19|20)\d{2}[a-z]?", p)), None)
    if yi is None:
        return s.replace("_", " ")
    out = f"{' '.join(parts[:yi])} et al. {parts[yi]}"
    if parts[yi + 1:]:
        out += " (" + " ".join(w.lower() for w in parts[yi + 1:]) + ")"
    return out


def pretty_taxon(s):
    """Render a SILVA label with per-token italics, as mathtext.

    The bacterial code italicises validly published names at every rank from
    class upward, so Clostridia, Christensenellaceae, Oscillospiraceae and
    Muribaculaceae are all italic. SILVA placeholder codes are not validly
    published names and stay upright: UCG-014, R-7, NK4A214, DTU089, CHKCI002,
    GCA-900066575. So do the words "group" and "clade", which are English.

    Deciding per token rather than per label is what makes composite SILVA
    labels come out right, e.g. Clostridia UCG-014 or Ruminococcus torques
    group. Hyphens are kept OUTSIDE math mode, because mathtext renders a
    hyphen as a minus sign.
    """
    import re
    UPRIGHT = {"group", "clade", "unclassified", "sensu", "stricto",
               "incertae", "sedis", "and", "or", "cluster", "complex"}
    s = re.sub(r"^[a-z]__", "", s).replace("[", "").replace("]", "")
    s = re.sub(r"\s+", " ", s.replace("_", " ").replace(".", " ")).strip()
    if s.lower().startswith("incertae sedis") or not s:
        return "Unclassified"

    def tok(t):
        # Decide on the whole space-delimited token: a placeholder code keeps
        # its hyphen (UCG-014, R-7, GCA-900066575), so test for a digit first.
        if not t or any(c.isdigit() for c in t) or t.lower() in UPRIGHT:
            return t
        # Fully alphabetic: italicise, splitting on any hyphen so that the
        # hyphen itself stays out of math mode (Escherichia-Shigella).
        return "-".join(r"$\mathit{%s}$" % w if w else "" for w in t.split("-"))

    # R's make.names turns "GCA-900066575" into "GCA.900066575", so a code and
    # its digit tail arrive as separate tokens. Rejoin a pure-digit token onto a
    # preceding token that carries no lowercase letters, i.e. a placeholder code
    # (GCA, UCG, R, NK4A214) rather than a Latin name.
    parts = s.split(" ")
    # make.names also prefixes a leading bracket with "X":
    # "[Eubacterium] coprostanoligenes group" -> "X.Eubacterium..coprostanoligenes.group"
    if len(parts) > 1 and parts[0] == "X":
        parts = parts[1:]
    toks = []
    for t in parts:
        if (toks and t.isdigit() and toks[-1] and not any(c.islower() for c in toks[-1])):
            toks[-1] = toks[-1] + "-" + t
        else:
            toks.append(t)
    return " ".join(tok(t) for t in toks)



def grid_on(ax, axis="x"):
    ax.grid(axis=axis, color=GRIDC, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    return ax


def load(root, labels=None):
    scores = os.path.join(root, "r_analysis", "ghi", "ghi_scores_pca.tsv")
    if not os.path.exists(scores):
        sys.exit(f"scores file not found: {scores}\nset PDI_ROOT or pass --root")
    g = pd.read_csv(scores, sep="\t")
    CLASS_COLS = ("Intervention_Class", "Class")
    if not any(c in g.columns for c in CLASS_COLS):
        if not labels:
            sys.exit("ghi_scores_pca.tsv carries no Intervention_Class column; "
                     "re-run 04b with the patched keep list, or pass --labels")
        key = ["Study_ID", "Treatment_Group", "Segment", "Age_d"]
        lab = pd.read_csv(labels)
        lcol = next((c for c in CLASS_COLS if c in lab.columns), None)
        if lcol is None:
            sys.exit(f"{labels} carries none of {CLASS_COLS}")
        for d in (g, lab):
            for k in key:
                d[k] = d[k].astype(str).str.strip()
        g = g.join(lab.drop_duplicates(key).set_index(key)[[lcol]].rename(
            columns={lcol: "Class"}), on=key)
    g = g.rename(columns={"Intervention_Class": "Class"})
    g = g[g.GHI_Combined_raw.notna()].copy()
    missing = g.Class.isna().sum()
    if missing:
        sys.exit(f"{missing} scored groups carry no intervention class")
    return g


def contrasts(g, treated, invert=False):
    arms = g.groupby(["Study_ID", "Treatment_Group", "Class"], as_index=False).GHI_Combined_raw.mean()
    if invert:
        arms = arms[~arms.Study_ID.isin(NO_CONCURRENT_CONTROL)]
    rows = []
    for s, d in arms.groupby("Study_ID"):
        c = d.loc[d.Class == CONTROL, "GHI_Combined_raw"]
        t = d.loc[d.Class.isin(treated), "GHI_Combined_raw"]
        if len(c) and len(t):
            e = t.mean() - c.mean()
            rows.append((s, -e if invert else e, len(c) + len(t)))
    out = pd.DataFrame(rows, columns=["study", "effect", "n_arms"])
    out["prec"] = 1 / np.sqrt(out.n_arms)
    return out


def gauge(x):
    """Raw composite to the 0-10 scale, clipped as script 61 clips it
    (529 scored, 3 clipped). Without the clip the additive class mean prints
    5.84 instead of 5.83."""
    return np.clip(5 + 1.5 * (np.asarray(x) + ANCHOR_MEAN) / ANCHOR_SD, 0, 10)


# ---------------------------------------------------------------- figures ---
def fig_class_means(g, out):
    tab = (g.groupby("Class")
             .agg(mean=("GHI_Combined_raw", lambda v: float(np.mean(gauge(v)))),
                  groups=("GHI_Combined_raw", "size"),
                  studies=("Study_ID", "nunique"))
             .reindex(CLASS_ORDER))
    colour = {CONTROL: GREY, TREATED: TRT_C, ADDITIVE: ADD_C}
    fig, ax = plt.subplots(figsize=(7.0, 4.6), layout="constrained")
    yy = np.arange(len(tab))[::-1]
    ax.axvline(5.00, color=INK, lw=1.0, ls=(0, (5, 3)), zorder=1)
    hi = 6.92                                   # counts sit past the last gridline
    for (nm, r), yv in zip(tab.iterrows(), yy):
        col = colour.get(nm, CLARET)
        ax.hlines(yv, 5.00, r["mean"], color=GREY_L, lw=1.2, zorder=2)
        ax.scatter([r["mean"]], [yv], s=44 + r["groups"] * 0.30, color=col, zorder=3, linewidths=0)
        ax.text(r["mean"] + (0.075 if r["mean"] >= 5 else -0.075), yv, f"{r['mean']:.2f}",
                fontsize=9.6, va="center", ha=("left" if r["mean"] >= 5 else "right"), color=INK)
        ax.text(hi - 0.40, yv, f"{int(r['groups'])}", fontsize=9.6, va="center", ha="right", color=INK)
        ax.text(hi, yv, f"{int(r['studies'])}", fontsize=9.6, va="center", ha="right", color=INK)
    ax.text(hi - 0.40, len(tab) - 0.26, "Groups", fontsize=9.2, ha="right", color=COLHEAD)
    ax.text(hi, len(tab) - 0.26, "Studies", fontsize=9.2, ha="right", color=COLHEAD)
    ax.set_yticks(yy); ax.set_yticklabels(tab.index, fontsize=10.2)
    ax.set_xticks([4.0, 4.5, 5.0, 5.5, 6.0])
    ax.set_xlim(float(tab["mean"].min()) - 0.32, hi + 0.12)
    ax.set_ylim(-0.75, len(tab) - 0.05)
    ax.set_xlabel("Gut Health Index  (0\u201310 scale)", labelpad=6)
    S.bare(ax); grid_on(ax, "x")
    return S.save(fig, "ghi_class_means", outdir=out)


def fig_domain_estimates(g, out):
    rows = []
    for lbl, col in DOMAIN_COLS:
        ya, wa = study_diffs(g, col, [ADDITIVE], [CONTROL])
        yc, wc = study_diffs(g, col, [CONTROL], CHALLENGE, skip=NO_CONCURRENT_CONTROL)
        rows.append((lbl, pool_re(ya, wa), pool_re(yc, wc)))
    fig, ax = plt.subplots(figsize=(6.2, 6.0), layout="constrained")
    h = 0.20
    for i, (lbl, a, c) in enumerate(rows):
        for r, col, off, mk in ((a, ADD_C, +h, "o"), (c, CLARET, -h, "s")):
            if r is None:
                continue
            if lbl == "Combined index":
                col = MAROON
            yv = len(rows) - 1 - i + off
            ax.plot([r["lo"], r["hi"]], [yv, yv], color=col, lw=1.8,
                    solid_capstyle="round", zorder=2)
            ax.scatter([r["est"]], [yv], s=38, color=col, marker=mk, zorder=3,
                       linewidths=0)
    ax.axvline(0, color=INK, lw=1.0, ls=(0, (5, 3)), zorder=1)
    ax.set_yticks(range(len(rows))[::-1])
    ax.set_yticklabels([r[0] for r in rows], fontsize=10.2)
    ax.set_xlabel("Standardized mean difference", labelpad=6)
    ax.scatter([], [], s=38, color=ADD_C, marker="o", label="Additive minus control")
    ax.scatter([], [], s=38, color=CLARET, marker="s", label="Control minus challenged")
    ax.set_ylim(-0.6, len(rows) + 0.55)          # a clear band at the top for the key
    leg = ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), fontsize=9.4,
                    handletextpad=0.4, borderaxespad=0, frameon=True, borderpad=0.6)
    leg.get_frame().set(edgecolor=GREY_L, facecolor="white", linewidth=0.9)
    S.bare(ax); grid_on(ax, "x")
    return S.save(fig, "ghi_domain_estimates", outdir=out)


def fig_roc(g, out):
    ctrl = g.loc[g.Class == CONTROL, "GHI_Combined_raw"].values
    chal = g.loc[g.Class.isin(CHALLENGE), "GHI_Combined_raw"].values
    y = np.r_[np.zeros(len(chal)), np.ones(len(ctrl))]
    y = y[np.argsort(-np.r_[chal, ctrl])]
    tpr = np.cumsum(y) / y.sum(); fpr = np.cumsum(1 - y) / (1 - y).sum()
    auc = np.trapezoid(np.r_[0, tpr], np.r_[0, fpr])
    fig, ax = plt.subplots(figsize=(4.6, 4.6), layout="constrained")
    ax.plot([0, 1], [0, 1], color=GREY_L, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.plot(np.r_[0, fpr], np.r_[0, tpr], color=MAROON, lw=2.0, zorder=3)
    ax.fill_between(np.r_[0, fpr], np.r_[0, tpr], color=MAROON, alpha=0.06, zorder=2)
    ax.text(0.98, 0.10, f"AUC {auc:.3f}", fontsize=11.5, ha="right", color=MAROON)
    ax.text(0.98, 0.028, "",
            fontsize=9.0, ha="right", color=GREY)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.set_xticks([0, 0.5, 1.0]); ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel("False positive rate", labelpad=6); ax.set_ylabel("True positive rate")
    S.bare(ax, spines=("bottom", "left")); grid_on(ax, "both")
    return S.save(fig, "ghi_roc", outdir=out)


def fig_funnel(g, out):
    add = contrasts(g, [ADDITIVE]); chal = contrasts(g, CHALLENGE, invert=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.6), layout="constrained")
    for d, col, mk, lb in ((add, ADD_C, "o", "Additive"),
                           (chal, CLARET, "s", "Challenge")):
        ax.scatter(d.effect, d.prec, s=FUNNEL_S, color=col, marker=mk, alpha=FUNNEL_A,
                   linewidths=0, zorder=3, label=lb)
        ax.axvline(float(d.effect.mean()), color=col, ls=(0, (4, 3)), lw=1.0, zorder=2)
    # Ring the INFLUENTIAL study, the one whose removal moves the funnel slope
    # most, not the one with the largest effect. Those are different studies here.
    base = st.linregress(chal.prec, chal.effect).slope
    infl = max(chal.study, key=lambda s: abs(
        st.linregress(chal.loc[chal.study != s, "prec"],
                      chal.loc[chal.study != s, "effect"]).slope - base))
    out_s = chal.loc[chal.study == infl].iloc[0]
    ax.scatter([out_s.effect], [out_s.prec], s=96, facecolors="none", edgecolors=MAROON,
               lw=1.4, zorder=4)
    ax.annotate(study_label(out_s.study), xy=(out_s.effect, out_s.prec),
                xytext=(out_s.effect + 0.95, out_s.prec + 0.055),
                fontsize=9.4, color=MAROON, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.9, color=MAROON, shrinkA=2, shrinkB=8))
    ax.invert_yaxis()
    span = float(max(abs(chal.effect.min()), abs(add.effect.max()), abs(chal.effect.max())))
    ax.set_xlim(-span - 0.8, span + 0.8)
    ax.set_xlabel("Within-study difference in the index", labelpad=6)
    ax.set_ylabel("Precision", fontsize=10.2)
    leg = ax.legend(loc="upper right", bbox_to_anchor=(4.0, 0.33),
                    bbox_transform=ax.transData, fontsize=9.4, handletextpad=0.4,
                    borderaxespad=0, frameon=True, borderpad=0.6)
    leg.get_frame().set(edgecolor=GREY_L, facecolor="white", linewidth=0.9)
    S.bare(ax, spines=("bottom", "left")); grid_on(ax, "both")
    return S.save(fig, "ghi_funnel", outdir=out)


def fig_domain_matrix(g, out):
    M = g[[f"GHI_dom_{d}" for d in DOMAINS]]
    C = np.full((7, 7), np.nan)
    for i in range(7):
        for j in range(7):
            a, b = M.iloc[:, i], M.iloc[:, j]
            ok = a.notna() & b.notna()
            if ok.sum() >= MIN_PAIRS:
                C[i, j] = st.spearmanr(a[ok], b[ok]).statistic
    fig, ax = plt.subplots(figsize=(6.8, 5.2), layout="constrained")
    im = ax.imshow(np.where(np.eye(7, dtype=bool), np.nan, C), cmap=CMAP,
                   vmin=-0.6, vmax=0.6, aspect="equal")
    for i in range(7):
        for j in range(7):
            if i == j:
                continue
            if np.isnan(C[i, j]):
                ax.text(j, i, "nd", fontsize=8.2, ha="center", va="center", color=GREY)
            else:
                ax.text(j, i, f"{C[i, j]:.2f}".replace("0.", "."), fontsize=8.6,
                        ha="center", va="center",
                        color=("white" if abs(C[i, j]) > 0.38 else INK))
    ax.set_xticks(range(7)); ax.set_xticklabels(DOMAIN_LABELS, rotation=90, ha="center", fontsize=9.6)
    ax.set_yticks(range(7)); ax.set_yticklabels(DOMAIN_LABELS, fontsize=9.6)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, ticks=[-0.6, 0, 0.6])
    cb.outline.set_visible(False); cb.ax.tick_params(length=0, labelsize=9.4)
    cb.set_label(r"Spearman $\rho$", fontsize=10.0)
    return S.save(fig, "ghi_domain_matrix", outdir=out)


def fig_distributions(g, out):
    """Every scored group on the 0-10 scale, with the class mean.

    Absolute scores, not deviations: the 0-10 scale anchored on the control mean
    is the instrument, and a negative axis would not be it. The spread is real
    within-study variation between arms, and the overlap between classes is a
    property of the published literature rather than of the index.
    """
    groups = [("Challenged", g[g.Class.isin(CHALLENGE)], CLARET),
              ("Control", g[g.Class == CONTROL], GREY),
              ("Treated", g[g.Class == TREATED], TRT_C),
              ("Additive", g[g.Class == ADDITIVE], ADD_C)]
    rng = np.random.default_rng(42)
    fig, ax = plt.subplots(figsize=(5.8, 4.3), layout="constrained")
    for i, (nm, d, col) in enumerate(groups):
        v = gauge(d.GHI_Combined_raw.values)
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=13, color=col,
                   alpha=0.42, linewidths=0, zorder=2)
        m = float(np.mean(v))
        ax.plot([i - 0.32, i + 0.32], [m, m], color=col, lw=2.4, zorder=4,
                solid_capstyle="round")
        ax.text(i, 10.75, f"{m:.2f}", fontsize=9.6, ha="center", color=col)
    ax.axhline(5.00, color=INK, lw=1.0, ls=(0, (5, 3)), zorder=1)
    ax.set_xticks(range(4)); ax.set_xticklabels([x[0] for x in groups], fontsize=10.2)
    ax.set_xlim(-0.62, 3.62); ax.set_ylim(-0.4, 11.9); ax.set_yticks([0, 5, 10])
    ax.set_ylabel("Gut Health Index")
    S.bare(ax, spines=("left",)); grid_on(ax, "y")
    return S.save(fig, "ghi_distributions", outdir=out)



def fig_pcoa(root, out):
    """Aitchison ordination of 1,578 runs, coloured by intervention family.

    Points are grouped into the four families the primary contrast uses, not
    all eight classes: with one class at n = 14 an eight-colour scatter of
    1,578 overlapping points communicates nothing. The eight-class PERMANOVA
    is reported in the text, not drawn here.

    Centroids carry a cross of +/- 1 standard error on each axis, which is the
    quantity a reader needs to judge whether the families separate at all.
    """
    st = pd.read_csv(_path(root, "gdi", "pcoa_strata.tsv"), sep="\t")
    ve = pd.read_csv(_path(root, "gdi", "pcoa_varexp.tsv"), sep="\t").set_index("axis").varexp_pct
    fam = {CONTROL: ("Control", GREY, "o"), ADDITIVE: ("Additive", ADD_C, "o"),
           TREATED: ("Challenged and treated", TRT_C, "s")}
    st["fam"] = st.Class.map(lambda c: fam[c][0] if c in fam else "Challenged")
    order = ["Challenged", "Control", "Challenged and treated", "Additive"]
    style = {"Challenged": (CLARET, "s"), "Control": (GREY, "o"),
             "Challenged and treated": (TRT_C, "s"), "Additive": (ADD_C, "o")}

    fig, ax = plt.subplots(figsize=(5.9, 5.0), layout="constrained")
    for nm in order:
        d = st[st.fam == nm]
        col, mk = style[nm]
        ax.scatter(d.PC1, d.PC2, s=9, color=col, alpha=0.30, marker=mk,
                   linewidths=0, zorder=2)
    for nm in order:
        d = st[st.fam == nm]
        col, mk = style[nm]
        cx, cy = d.PC1.mean(), d.PC2.mean()
        ex, ey = d.PC1.std()/np.sqrt(len(d)), d.PC2.std()/np.sqrt(len(d))
        ax.plot([cx-ex, cx+ex], [cy, cy], color=col, lw=1.6, zorder=5,
                solid_capstyle="round")
        ax.plot([cx, cx], [cy-ey, cy+ey], color=col, lw=1.6, zorder=5,
                solid_capstyle="round")
        ax.scatter([cx], [cy], s=84, color=col, marker=mk, zorder=6,
                   edgecolors="white", linewidths=1.3,
                   label=f"{nm}  ({len(d)})")
    ax.axhline(0, color=GRIDC, lw=0.8, zorder=1)
    ax.axvline(0, color=GRIDC, lw=0.8, zorder=1)
    ax.set_xlabel(f"PC1  ({ve['PC1']:.1f}% of variance)")
    ax.set_ylabel(f"PC2  ({ve['PC2']:.1f}% of variance)")
    leg = ax.legend(loc="lower left", frameon=True, borderpad=0.6, fontsize=8.6,
                    handletextpad=0.5, labelspacing=0.5)
    leg.get_frame().set(edgecolor=GREY_L, facecolor="white", linewidth=0.9)
    S.clean(ax)
    return S.save(fig, "gdi_pcoa", outdir=out)


def _path(root, *parts):
    return os.path.join(root, "r_analysis", *parts)


def fig_domain_coverage(g, out):
    rows = []
    for lbl, col in DOMAIN_COLS[:-1]:
        d = g[g[col].notna()]
        rows.append((lbl, len(d), d.Study_ID.nunique()))
    rows.sort(key=lambda r: r[1])
    fig, ax = plt.subplots(figsize=(6.4, 4.0), layout="constrained")
    yy = np.arange(len(rows))
    ax.barh(yy, [r[1] for r in rows], height=0.52, color=S.BLUE_L, zorder=2)
    for i, (lbl, ng, ns) in enumerate(rows):
        ax.text(ng + 6, i, f"{ng}", fontsize=9.6, va="center", ha="left", color=INK)
        ax.text(ng + 46, i, f"{ns} studies", fontsize=9.0, va="center", ha="left",
                color=COLHEAD)
    ax.set_yticks(yy); ax.set_yticklabels([r[0] for r in rows], fontsize=10.2)
    ax.set_xlim(0, max(r[1] for r in rows) * 1.42)
    ax.set_xlabel("Groups with the domain measured", labelpad=6)
    S.bare(ax); grid_on(ax, "x")
    return S.save(fig, "ghi_domain_coverage", outdir=out)


def fig_contrast_forest(g, out):
    arms = g.groupby(["Study_ID", "Treatment_Group", "Class"],
                     as_index=False).GHI_Combined_raw.mean()
    panels = [("Additive minus control", [ADDITIVE], [CONTROL], ADD_C, ()),
              ("Control minus challenged", [CONTROL], CHALLENGE, CLARET,
               NO_CONCURRENT_CONTROL)]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 7.4), layout="constrained")
    for ax, (ttl, better, worse, col, skip) in zip(axes, panels):
        rows = []
        for st_id, d in arms.groupby("Study_ID"):
            if st_id in skip:
                continue
            a = d.loc[d.Class.isin(better), "GHI_Combined_raw"]
            b = d.loc[d.Class.isin(worse), "GHI_Combined_raw"]
            if len(a) and len(b):
                rows.append((st_id, a.mean() - b.mean(), len(a) + len(b)))
        t = pd.DataFrame(rows, columns=["study", "effect", "n"]).sort_values("effect")
        yy = np.arange(len(t))
        ax.axvline(0, color=GREY_L, lw=1.1, zorder=1)
        ax.hlines(yy, 0, t.effect, color=GREY_L, lw=1.1, zorder=2)
        ax.scatter(t.effect, yy, s=np.clip(t.n * 7, 18, 70), color=col, zorder=3,
                   linewidths=0)
        m = float(t.effect.mean())
        ax.axvline(m, color=MAROON, ls=(0, (4, 3)), lw=1.3, zorder=4)
        ax.set_yticks(yy); ax.set_yticklabels([study_label(x) for x in t.study], fontsize=8.0)
        ax.set_ylim(-1, len(t))
        ax.set_xlabel("Difference in the index", labelpad=6)
        ax.text(0.0, 1.015, ttl, transform=ax.transAxes, fontsize=9.4, color=INK,
                va="bottom", ha="left")
        S.bare(ax); grid_on(ax, "x")
    return S.save(fig, "ghi_contrast_forest", outdir=out)


def fig_loso_scatter(root, out):
    f = _path(root, "gdi", "loso_predictions.tsv")
    if not os.path.exists(f):
        return []
    d = pd.read_csv(f, sep="\t").dropna(subset=["ghi_centered", "gdi_centered"])
    r = st.spearmanr(d.gdi_centered, d.ghi_centered)
    fig, ax = plt.subplots(figsize=(5.2, 4.6), layout="constrained")
    ax.axhline(0, color=GREY_L, lw=1.0, zorder=1)
    ax.axvline(0, color=GREY_L, lw=1.0, zorder=1)
    ax.scatter(d.gdi_centered, d.ghi_centered, s=26, color=CLARET, alpha=0.75,
               linewidths=0, zorder=3)
    ax.text(0.015, 0.975, f"Spearman $\\rho$ = {r.statistic:+.3f}",
            transform=ax.transAxes, fontsize=10.4, ha="left", va="top", color=MAROON)
    ax.set_xlabel("Predicted index", labelpad=6)
    ax.set_ylabel("Measured index")
    S.bare(ax, spines=("bottom", "left")); grid_on(ax, "both")
    return S.save(fig, "gdi_loso_scatter", outdir=out)


def fig_loso_per_study(root, out):
    f = _path(root, "gdi", "loso_per_study.tsv")
    if not os.path.exists(f):
        return []
    d = pd.read_csv(f, sep="\t").dropna(subset=["rho"]).sort_values("rho")
    fig, ax = plt.subplots(figsize=(6.4, 7.0), layout="constrained")
    yy = np.arange(len(d))
    ax.axvline(0, color=INK, lw=1.0, ls=(0, (5, 3)), zorder=1)
    ax.hlines(yy, 0, d.rho, color=GREY_L, lw=1.1, zorder=2)
    cols = [CLARET if v < 0 else ADD_C for v in d.rho]
    ax.scatter(d.rho, yy, s=np.clip(d.n * 5, 18, 70), color=cols, zorder=3, linewidths=0)
    ax.set_yticks(yy); ax.set_yticklabels([study_label(x) for x in d.study], fontsize=8.0)
    ax.set_xlim(-1.08, 1.08); ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_ylim(-1, len(d))
    ax.set_xlabel("Within-study Spearman $\\rho$", labelpad=6)
    S.bare(ax); grid_on(ax, "x")
    return S.save(fig, "gdi_loso_per_study", outdir=out)


def fig_guild_scfa(root, out):
    f = _path(root, "verification", "scfa_guild.csv")
    if not os.path.exists(f):
        return []
    d = pd.read_csv(f)
    fig, ax = plt.subplots(figsize=(5.8, 3.0), layout="constrained")
    yy = np.arange(len(d))[::-1]
    ax.axvline(0, color=GREY_L, lw=1.1, zorder=1)
    for (_, r), yv in zip(d.iterrows(), yy):
        col = SAGE if r["BH_q"] < 0.05 else GREY
        ax.hlines(yv, 0, r["rho"], color=GREY_L, lw=1.1, zorder=2)
        ax.scatter([r["rho"]], [yv], s=46, color=col, zorder=3, linewidths=0)
        ax.text(0.62, yv, f"$\\rho$ = {r['rho']:+.3f}    $q$ = {r['BH_q']:.3f}",
                fontsize=9.2, va="center", ha="left", color=col)   # fixed column, so a
        # negative rho never pushes its label back across the category name
    ax.set_yticks(yy); ax.set_yticklabels([x.capitalize() for x in d.scfa], fontsize=10.2)
    ax.set_xlim(-0.25, 0.92); ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_xlabel("Spearman $\\rho$", labelpad=6)
    S.bare(ax)
    return S.save(fig, "gdi_guild_scfa", outdir=out)


def fig_permanova(root, out):
    f = _path(root, "verification", "permanova_five_configs.csv")
    if not os.path.exists(f):
        return []
    d = pd.read_csv(f)
    r = d.iloc[0]                       # the eligible set; the rest are sensitivity runs
    fig, ax = plt.subplots(figsize=(5.0, 4.0), layout="constrained")
    vals = [("Gut Health Index", r.GHI_R2 * 100, r.GHI_p, CLARET),
            ("Age", r.Age_R2 * 100, r.Age_p, S.BLUE)]
    for i, (nm, v, pv, col) in enumerate(vals):
        ax.bar(i, v, width=0.46, color=col, zorder=2)
        ptxt = "$P$ < 0.001" if pv < 0.001 else f"$P$ = {pv:.3f}"
        ax.text(i, v + 0.35, f"{v:.2f}%", fontsize=10.0, ha="center", color=col)
        ax.text(i, v + 1.05, ptxt, fontsize=9.0, ha="center", color=COLHEAD)
    ax.set_xticks(range(2)); ax.set_xticklabels([v[0] for v in vals], fontsize=10.2)
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, max(v[1] for v in vals) * 1.32)
    ax.set_ylabel("Community variation explained (%)")
    S.bare(ax); grid_on(ax, "y")
    return S.save(fig, "gdi_permanova", outdir=out)


def fig_stability(root, out):
    f = _path(root, "gdi", "gdi_stability.tsv")
    if not os.path.exists(f):
        return []
    d = pd.read_csv(f, sep="\t").sort_values("stability")
    d = d.assign(genus=d.genus.str.replace(r"^g__", "", regex=True))
    v = d.stability * (100 if d.stability.max() <= 1 else 1)
    fig, ax = plt.subplots(figsize=(6.0, 6.4), layout="constrained")
    yy = np.arange(len(d))
    ax.barh(yy, v, height=0.56, color=S.BLUE_L, zorder=2)
    labs = [pretty_taxon(x) for x in d.genus]
    ax.axvline(70, color=MAROON, lw=1.2, ls=(0, (5, 3)), zorder=4)
    ax.text(70, len(d) - 0.2, "70% selection threshold", fontsize=9.2, color=MAROON,
            ha="right", va="bottom")
    ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=8.2)
    ax.set_xlim(0, 100); ax.set_ylim(-0.8, len(d) + 0.4)
    ax.set_xlabel("Selection frequency (%)", labelpad=6)
    S.bare(ax); grid_on(ax, "x")
    return S.save(fig, "gdi_stability", outdir=out)


def fig_maaslin(root, out, top=20):
    f = _path(root, "gdi", "maaslin2_all_results.tsv")
    if not os.path.exists(f):
        return []
    d = pd.read_csv(f, sep="\t").dropna(subset=["coef", "qval"])
    n_all, n_sig = len(d), int((d.qval < 0.25).sum())
    d = d.assign(feature=d.feature.str.replace(r"^g__", "", regex=True)
                                  .str.replace(".", " ", regex=False).str.strip())
    d = d.reindex(d.coef.abs().sort_values(ascending=False).index).head(top)
    d = d.sort_values("coef")
    fig, ax = plt.subplots(figsize=(6.4, 6.0), layout="constrained")
    yy = np.arange(len(d))
    ax.axvline(0, color=INK, lw=1.0, ls=(0, (5, 3)), zorder=1)
    ax.hlines(yy, 0, d.coef, color=GREY_L, lw=1.1, zorder=2)
    cols = [MAROON if q < 0.25 else GREY for q in d.qval]
    ax.scatter(d.coef, yy, s=36, color=cols, zorder=3, linewidths=0)
    labs = [pretty_taxon(x) for x in d.feature]
    ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=8.2)
    ax.set_ylim(-1, len(d))
    ax.set_xlabel("MaAsLin2 coefficient", labelpad=6)
    S.bare(ax); grid_on(ax, "x")
    return S.save(fig, "gdi_maaslin_ranked", outdir=out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.environ.get("PDI_ROOT", os.path.expandvars("/scratch/user/$USER/pdi")))
    ap.add_argument("--labels", default=None)
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    out = a.outdir or os.path.join(a.root, "figures")
    g = load(a.root, a.labels)
    written = []
    for fn in (lambda: fig_class_means(g, out), lambda: fig_domain_estimates(g, out),
               lambda: fig_roc(g, out), lambda: fig_funnel(g, out),
               lambda: fig_domain_matrix(g, out), lambda: fig_distributions(g, out),
               lambda: fig_domain_coverage(g, out), lambda: fig_contrast_forest(g, out),
               lambda: fig_loso_scatter(a.root, out), lambda: fig_loso_per_study(a.root, out),
               lambda: fig_guild_scfa(a.root, out), lambda: fig_permanova(a.root, out),
               lambda: fig_stability(a.root, out), lambda: fig_maaslin(a.root, out),
               lambda: fig_pcoa(a.root, out)):
        try:
            written += fn()
        except FileNotFoundError as e:
            print(f"SKIP: {os.path.basename(str(getattr(e, 'filename', '') or ''))} not found")
        plt.close("all")
    print(f"{len(g)} scored groups, {g.Study_ID.nunique()} studies, "
          f"{g.Class.nunique()} classes")
    for p in written:
        print(p)


if __name__ == "__main__":
    main()
