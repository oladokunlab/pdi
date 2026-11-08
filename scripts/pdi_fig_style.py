"""PDI figure style. Import at the top of every figure script.

    import pdi_fig_style as S
    S.apply()
    fig, ax = plt.subplots(figsize=(S.W1, 3.4))
    ...
    S.clean(ax, grid="x")
    S.save(fig, "fig03_forest_contrasts")

Fonts ship in ./fonts next to this file, so the figures render identically on
Grace and on a laptop with no system font installation.
"""
import os
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager as fm

_HERE = os.path.dirname(os.path.abspath(__file__))
FONT_DIR = next((d for d in (os.path.join(_HERE, "fonts"), os.path.join(_HERE, os.pardir, "fonts")) if os.path.isdir(d)), os.path.join(_HERE, "fonts"))
FAMILY = "DejaVu Sans"
SERIF = "DejaVu Serif"
if os.path.isdir(FONT_DIR):
    for _f in sorted(os.listdir(FONT_DIR)):
        if _f.lower().endswith((".ttf", ".otf")):
            try:
                fm.fontManager.addfont(os.path.join(FONT_DIR, _f))
            except Exception:
                pass
    _names = {f.name for f in fm.fontManager.ttflist}
    for _candidate in ("Nimbus Sans", "Source Sans 3"):   # Nimbus Sans carries
        if _candidate in _names:                          # Helvetica's metrics and
            FAMILY = _candidate                           # ships under AGPL v3 with
            break                                         # a PDF-embedding exemption
    if "Source Serif 4" in _names:
        SERIF = "Source Serif 4"

# ---- palette -------------------------------------------------------------
INK      = "#002147"   # Oxford blue. ALL text, every label, every number.
MAROON   = "#500000"   # Aggie maroon. One element per panel, never a series.
CLARET   = "#8E4A50"   # maroon mid, challenge series
CLARET_L = "#C9A3A5"
BLUE_D   = "#2E5A7D"
BLUE     = "#4A7BA0"   # challenged and treated
BLUE_L   = "#A8C2D6"
TEAL     = "#2E6E75"   # non-disruptive additive
SAGE     = "#4F7A5E"   # the positive sequence result
AMBER    = "#B07C3C"   # fourth category only
GREY     = "#7E8C99"   # control, and all secondary text
GREY_L   = "#C3CDD5"   # leader lines, reference rules
GRID     = "#EAEFF3"
SPINE    = "#CFD8DF"
CHANCE   = "#C3CDD5"   # anything at chance

#: concept -> colour. Never reassign; a concept keeps its hue in every figure.
SEMANTIC = {
    "control": GREY, "challenge": CLARET, "treated": BLUE,
    "additive": TEAL, "sequence": SAGE, "chance": CHANCE, "emphasis": MAROON,
}
DIVERGING = ["#7A3438", "#C08E90", "#F2F5F7", "#8FB2CA", "#1E4A6B"]

# ---- type ladder, points -------------------------------------------------
T_PANEL, T_LABEL, T_TICK, T_ANNOT, T_SECOND = 11.0, 10.0, 9.2, 8.6, 8.2

# ---- column widths, inches (single and double column) --------------------
W1, W15, W2 = 3.46, 5.12, 7.09


def apply(base=10.0):
    matplotlib.rcParams.update({
        "font.family": FAMILY, "font.size": base, "text.color": INK, "mathtext.fontset": "custom", "mathtext.rm": FAMILY, "mathtext.it": FAMILY + ":italic", "mathtext.bf": FAMILY + ":bold",
        "axes.labelcolor": INK, "axes.edgecolor": SPINE, "axes.linewidth": 0.8,
        "xtick.color": INK, "ytick.color": INK,
        "axes.labelsize": base, "axes.titlesize": base + 1.0,
        "xtick.labelsize": base - 0.8, "ytick.labelsize": base - 0.8,
        "legend.fontsize": base - 1.0, "legend.frameon": False,
        "axes.unicode_minus": True, "figure.dpi": 110, "savefig.dpi": 600, "axes.axisbelow": True,
        "pdf.fonttype": 42, "ps.fonttype": 42,  # embed as TrueType, not outlines
        "svg.fonttype": "none",
    })


def clean(ax, grid="y", spines=("bottom",)):
    """Strip the frame to the axes that carry information."""
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in spines)
        if s in spines:
            ax.spines[s].set_color(SPINE)
    if grid in ("x", "both"):
        ax.grid(axis="x", color=GRID, lw=0.9, zorder=0)
    if grid in ("y", "both"):
        ax.grid(axis="y", color=GRID, lw=0.9, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, color=SPINE)
    if "left" not in spines:
        ax.tick_params(axis="y", length=0)
    if "bottom" not in spines:
        ax.tick_params(axis="x", length=0)
    return ax


RAIL_W = 0.40  # rail width as a fraction of the plot width


def rail(fig, cell, letter, head, note=None, wrap=19):
    """Split a gridspec cell into a left rail and a plot.

    The rail carries the panel letter, the head and the method note, and closes
    with a maroon hairline. Nothing is set above the plot, so every panel in the
    paper has the same silhouette and the rails line up down the page.

    Category labels must be placed INSIDE the returned axes (flush left, in data
    coordinates) and y tick labels switched off, or they will be drawn leftward
    across the rail.
    """
    import textwrap
    sub = cell.subgridspec(1, 2, width_ratios=[RAIL_W, 1.0], wspace=0.0)
    axl = fig.add_subplot(sub[0])
    axl.axis("off")
    axl.set_xlim(0, 1)
    axl.set_ylim(0, 1)
    axl.text(0.0, 0.985, letter, fontsize=T_PANEL + 1.5, color=MAROON,
             weight="semibold", va="top", ha="left")
    axl.text(0.0, 0.885, "\n".join(textwrap.wrap(head, wrap)), fontsize=T_LABEL,
             color=INK, va="top", ha="left", linespacing=1.32)
    if note:
        axl.text(0.0, 0.0, "\n".join(textwrap.wrap(note, wrap + 3)),
                 fontsize=T_SECOND, color=GREY, va="bottom", ha="left",
                 linespacing=1.38)
    axl.plot([0.955, 0.955], [0.0, 1.0], transform=axl.transAxes, color=MAROON,
             lw=0.9, clip_on=False, zorder=5)
    return fig.add_subplot(sub[1])


def bare(ax, spines=("bottom",)):
    """No gridlines at all. Values are printed at the marks instead."""
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in spines)
        if s in spines:
            ax.spines[s].set_color(SPINE)
    ax.grid(False)
    ax.tick_params(length=3, color=SPINE)
    if "left" not in spines:
        ax.tick_params(axis="y", length=0)
    if "bottom" not in spines:
        ax.tick_params(axis="x", length=0)
    return ax


def letter(ax, ch, y=1.035, size=13.0):
    """Panel letter alone. No head, no rule: the description belongs in the
    manuscript caption, not inside the image."""
    ax.text(0.0, y, ch, fontsize=size, weight="bold", color=INK,
            transform=ax.transAxes, va="baseline", ha="left")


def panel(ax, letter, head=None, y=1.075, rule=True, gap_pt=15):
    """Plate head: serif letter, serif head set a fixed 15 pt clear of it, and a
    hairline rule across the panel. The fixed point gap is deliberate; an axes
    fraction offset closes up in narrow panels and the letter collides."""
    from matplotlib.transforms import ScaledTranslation
    ax.text(0.0, y, letter, fontsize=T_PANEL, color=INK, transform=ax.transAxes,
            va="baseline", weight="semibold")
    if head:
        tr = ax.transAxes + ScaledTranslation(gap_pt / 72.0, 0,
                                              ax.figure.dpi_scale_trans)
        ax.text(0.0, y, head, fontsize=T_LABEL, color=INK, transform=tr,
                va="baseline")
    if rule:
        ax.plot([0, 1], [y - 0.038, y - 0.038], transform=ax.transAxes,
                color=SPINE, lw=0.7, clip_on=False, zorder=5)


def reference(ax, x=None, y=None, label=None, color=None):
    """A reference rule (control anchor, chance, zero) drawn behind the data."""
    color = color or GREY_L
    if x is not None:
        ax.axvline(x, color=color, lw=1.1, zorder=1)
    if y is not None:
        ax.axhline(y, color=color, lw=1.1, zorder=1)
    if label:
        lo, hi = ax.get_ylim() if x is not None else ax.get_xlim()
        if x is not None:
            ax.text(x, hi, label, fontsize=T_ANNOT, ha="center", va="bottom",
                    color=MAROON)
        else:
            ax.text(hi, y, label, fontsize=T_ANNOT, ha="right", va="bottom",
                    color=MAROON)


def frame(fig, color=None, lw=0.9):
    """Thin grey rule around the whole figure, so it reads as its own object
    rather than as part of the page it is placed on."""
    from matplotlib.patches import Rectangle
    fig.add_artist(Rectangle((0.002, 0.002), 0.996, 0.996, transform=fig.transFigure,
                             fill=False, edgecolor=(color or GREY_L), lw=lw, zorder=10))
    return fig


def save(fig, stem, outdir="figures", formats=("png", "pdf")):
    os.makedirs(outdir, exist_ok=True)
    paths = []
    for ext in formats:
        p = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(p, facecolor="white")   # no tight bbox: the frame() rule
        # must enclose the labels, so the figure area has to contain them
        paths.append(p)
    return paths
