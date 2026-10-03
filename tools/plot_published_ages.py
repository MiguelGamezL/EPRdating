"""Figure of the published-age benchmark: EPRdating vs the published ages.

    python tools/plot_published_ages.py        # writes docs/img/published_ages.png

Each point is one published age (US-ESR, CS-US, EU), recomputed from the
published inputs with the conventions of the program the authors used
(tests/validation/published_cases.py). Hollow points: Khok Sung US-ESR ages
whose published dentine p does not reproduce the measured dentine ratio.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from eprdating.plot import CATEGORICAL, GRID, INK, INK_2, INK_3, SURFACE

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "validation"))
from published_cases import benchmark

COLOR = {"DATA": CATEGORICAL[0], "USESR": CATEGORICAL[1], "ROSY": CATEGORICAL[2]}
MARK = {"US": "o", "CSUS": "s", "EU": "^", "LU": "v"}


def main(out=ROOT / "docs" / "img" / "published_ages.png"):
    warnings.filterwarnings("ignore")
    rows = [r for r in benchmark() if r.get("ours")]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axhspan(-5, 5, color=GRID, alpha=0.6, lw=0, zorder=0)
    ax.axhline(0, color=INK_3, lw=1, zorder=1)
    for r in rows:
        dev = 100 * (r["ours"] / r["published"] - 1)
        unconv = "dentine p" in (r.get("excluded") or "")
        c = COLOR[r["program"]]
        ax.scatter(r["published"], dev, s=34, marker=MARK[r["model"].upper()], zorder=3,
                   facecolor="none" if unconv else c, edgecolor=c, linewidth=1.4 if unconv else 0.6)
    ax.set_xscale("log")
    ticks = [2, 5, 10, 20, 50, 100, 200, 500, 1000]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel("Published age (ka)", color=INK_2)
    ax.set_ylabel("EPRdating − published (%)", color=INK_2)
    ax.tick_params(colors=INK_2)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_ylim(-15, 15)
    handles = [plt.Line2D([], [], ls="", marker="o", color=c, label=p) for p, c in COLOR.items()]
    handles += [plt.Line2D([], [], ls="", marker=m, color=INK_3, label=k) for k, m in MARK.items()]
    handles.append(plt.Line2D([], [], ls="", marker="o", markerfacecolor="none", markeredgecolor=INK_3,
                              label="DATA dentine p\nnot converged"))
    leg = ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8.5)
    for t in leg.get_texts():
        t.set_color(INK)
    ax.text(0.01, 0.97, "shaded: ±5 %", transform=ax.transAxes, va="top", fontsize=8, color=INK_3)
    fig.tight_layout()
    fig.savefig(out, dpi=160, facecolor=SURFACE)
    print("wrote", out)


if __name__ == "__main__":
    main()
