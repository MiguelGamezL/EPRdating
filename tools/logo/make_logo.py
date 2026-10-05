"""EPRdating logo: the CO2- spectrum of tooth enamel with the name.

    python tools/logo/make_logo.py            # writes docs/img/logo/

The spectrum (spectrum.csv) is the largest of an additive-dose series,
digitised without smoothing. "EPR" ends where the peak starts to rise and
"dating" starts where its falling flank crosses the baseline; the double
minimum hangs under "dating". Colours follow the EPRAYA logo; the mark of the
Grupo de Física Aplicada (UNAL) sits in the bottom-right corner.

Font: Barlow Medium (SIL Open Font License), from
https://github.com/google/fonts/tree/main/ofl/barlow; put Barlow-Medium.ttf
next to this script.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "docs" / "img" / "logo"
NAVY, BLUE = "#002060", "#00a2e8"  # EPRAYA's colours
DARK_BG = "#0d1117"
plt.rcParams["svg.fonttype"] = "path"

X, Y = np.loadtxt(HERE / "spectrum.csv", delimiter=",", unpack=True)
FONT = FontProperties(fname=str(HERE / "Barlow-Medium.ttf"))
MARK = plt.imread(str(HERE / "group_mark.png"))

SIZE, SPAN, AMP, LW, PAD, MARGIN = 1.25, 11.5, 3.8, 5.5, 0.14, 0.45  # inches
MARK_H = 0.95


def _glyphs(s):
    tp = TextPath((0, 0), s, prop=FONT, size=1.0)
    return tp.transformed(Affine2D().scale(SIZE)), tp.get_extents()


def _word(ax, s, x, y, color, bg):
    path, ext = _glyphs(s)
    patch = PathPatch(path.transformed(Affine2D().translate(x - ext.x0 * SIZE, y)), facecolor=color,
                      edgecolor="none", zorder=3)
    if bg:
        patch.set_path_effects([pe.Stroke(linewidth=10, foreground=bg), pe.Normal()])
    ax.add_patch(patch)


def logo(name, ink=NAVY, bg="white", mark=True):
    ip = int(np.argmax(Y))
    rise = X[np.where(Y[:ip] < 0.12)[0][-1]]  # foot of the rising flank
    fall = X[ip + np.where(Y[ip:] < 0.0)[0][0]]  # falling flank crosses the baseline
    w1 = _glyphs("EPR")[1].width * SIZE
    w2 = _glyphs("dating")[1].width * SIZE
    x0 = MARGIN
    epr_right = x0 + rise * SPAN - PAD
    dat_left = x0 + fall * SPAN + PAD

    def under(a, b):
        return Y[(X >= (a - x0) / SPAN) & (X <= (b - x0) / SPAN)]

    lift = max(under(epr_right - w1, epr_right).max(), under(dat_left, dat_left + w2).max(), 0) * AMP + 0.16
    W = max(x0 + SPAN + MARGIN, dat_left + w2 + MARGIN)
    base_y = -Y.min() * AMP + 0.4
    H = base_y + AMP + 0.5

    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    if bg:
        fig.patch.set_facecolor(bg)
    else:
        fig.patch.set_alpha(0)
    ax.plot(x0 + X * SPAN, base_y + Y * AMP, color=BLUE, lw=LW, solid_capstyle="round", solid_joinstyle="round",
            zorder=2)
    if mark:
        mw = MARK_H * MARK.shape[1] / MARK.shape[0]
        ax.imshow(MARK, extent=(x0 + SPAN - mw, x0 + SPAN, 0.35, 0.35 + MARK_H), zorder=5,
                  interpolation="lanczos")
    _word(ax, "EPR", epr_right - w1, base_y + lift, ink, bg)
    _word(ax, "dating", dat_left, base_y + lift, ink, bg)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("auto")
    ax.axis("off")
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=150, transparent=not bg, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)


def icon(name="icon"):
    """Square icon: the spectrum alone (favicon, avatars)."""
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.06, 0.06, 0.88, 0.88])
    fig.patch.set_alpha(0)
    ax.plot(X, Y, color=BLUE, lw=14, solid_capstyle="round", solid_joinstyle="round")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(Y.min() - 0.08, Y.max() + 0.08)
    ax.axis("off")
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=128, transparent=True)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    logo("eprdating_logo")
    logo("eprdating_logo_dark", ink="white", bg=DARK_BG)
    logo("eprdating_logo_transparent", bg=None)
    icon()
    print("wrote", *sorted(p.name for p in OUT.iterdir()))
