"""Figures of the EPRdating article (PDF), from the validation data.

    python tools/paper_figures.py OUTDIR [--m18 FOLDER]

fig2_m18.pdf needs the M18 spectra (not in the repository); the others use
only the reference data under tests/validation.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from eprdating import Sediment, ToothLayers, ToothSample, USModel, fit_dose_response, plot
from eprdating.dose_rate import cosmic_dose_rate_sea_level
from eprdating.plot import CATEGORICAL, GRID, INK, INK_2, INK_3, SURFACE
from eprdating.usesr import th230_u234

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / "tests" / "validation"
sys.path.insert(0, str(VAL))
sys.path.insert(0, str(ROOT / "examples"))

plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8, "pdf.fonttype": 42})


def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_3)
    ax.tick_params(colors=INK_2)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def _logticks(ax, ticks):
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.minorticks_off()


def _save(fig, out):
    fig.savefig(out, facecolor=SURFACE)
    fig.savefig(out.with_suffix(".png"), dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print("wrote", out)


# --- Fig. 2: M18 ----------------------------------------------------------------
def fig_m18(folder: Path, out: Path):
    import dose_series_dat as ds

    rows, shown = ds.measure_series(folder)  # default window: 100 G around g = 2.0023
    used = [r for r in rows if r[0] not in ds.EXCLUDED]
    drc = fit_dose_response([r[1] for r in used], [r[2] for r in used], "LIN", sigma=[r[3] for r in used],
                            De_min=-np.inf)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.6), facecolor=SURFACE, gridspec_kw={"width_ratios": [1, 1.15]})
    keep = [f for f in ds.SERIES if f not in ds.EXCLUDED]
    spectra = [(shown[f][0], shown[f][1]) for f in keep]
    fits = [(shown[f][0], shown[f][2]) for f in keep]
    plot.plot_spectra(spectra, [f"{ds.SERIES[f]:.0f} Gy" for f in keep], fits=fits, ax=a1)
    a1.set_title("a  Spectra and template fits", loc="left", fontsize=9, color=INK)
    ex = [(r[1], r[2], r[3]) for r in rows if r[0] in ds.EXCLUDED]
    plot.plot_dose_response(drc, excluded=ex, ax=a2, intensity_label="CO$_2^-$ amplitude (a.u.)")
    a2.set_xlim(-140, 205)  # keeps the De label inside the axes
    a2.set_title("b  Dose response (linear)", loc="left", fontsize=9, color=INK)
    fig.tight_layout(w_pad=3.0)
    _save(fig, out)
    return drc


# --- Fig. 3: ROSY --------------------------------------------------------------------
def fig_rosy(out: Path):
    import test_rosy_reference as tr

    pts = []
    for cid in tr.AGE_CASES:
        p = tr.CASES[cid]
        geo = ToothLayers(enamel_um=tr._v(p["thick_en"]), dentine_um=tr._v(p["thick_den"]),
                          strip_outer_um=tr._v(p["strip_out"]), strip_inner_um=tr._v(p["strip_in"]),
                          enamel_density=tr._v(p["density_en"]), dentine_density=tr._v(p["density_den"]),
                          sediment_density=tr._v(p["density_sed"]))
        if p["env_option"] == 2:
            env = {"gamma": None, "cosmic": cosmic_dose_rate_sea_level(tr._v(p["depth"]), tr._v(p["overburden_density"]))}
        else:
            env = {"gamma": tr._v(p["gamma_cosmic"]) / 1000, "cosmic": 0.0}
        cu = (-1.0 if p["uptake_en"] == 1 else 0.0, -1.0 if p["uptake_den"] == 1 else 0.0)
        for mode, (pe, pd) in (("EU", (-1.0, -1.0)), ("LU", (0.0, 0.0)), ("CU", cu)):
            s = ToothSample(De=tr._v(p["De"]), enamel_U=tr._v(p["U_en"]), dentine_U=tr._v(p["U_den"]),
                            sediment=Sediment(U=tr._v(p["U_sed"]), Th=tr._v(p["Th_sed"]), K=tr._v(p["K_sed"]),
                                              water=tr._v(p["water_sed"]) / 100),
                            beta=geo, k_alpha=tr._v(p["alpha_eff"]), alpha_efficiency="energy",
                            uptake_enamel=USModel(pe), uptake_dentine=USModel(pd), factors="adamiec_aitken_1998",
                            **env)
            rosy = tr.RES[cid][mode]["age"] / 1000
            pts.append((mode, rosy, 100 * (s.age().age / rosy - 1), cid.startswith("brennan")))
    fig, ax = plt.subplots(figsize=(3.6, 3.0), facecolor=SURFACE)
    _style(ax)
    ax.axhspan(-1, 1, color=GRID, alpha=0.7, lw=0)
    ax.axhline(0, color=INK_3, lw=0.8)
    for mode, mk, c in (("EU", "o", CATEGORICAL[0]), ("LU", "s", CATEGORICAL[1]), ("CU", "^", CATEGORICAL[2])):
        xs = [p[1] for p in pts if p[0] == mode]
        ys = [p[2] for p in pts if p[0] == mode]
        ax.scatter(xs, ys, s=16, marker=mk, color=c, label=mode, zorder=3, lw=0)
    _logticks(ax, [20, 50, 100, 200, 500])
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel("ROSY 2.0 age (ka)")
    ax.set_ylabel("EPRdating − ROSY (%)")
    ax.legend(frameon=False, loc="lower right", ncols=3)
    fig.tight_layout()
    _save(fig, out)
    d = [p[2] for p in pts]
    print(f"  ROSY: n = {len(d)}, {min(d):+.2f} to {max(d):+.2f} %")


# --- Fig. 4: beta per segment vs one factor, against DATA -----------------------------
def fig_beta(out: Path):
    import test_data_reference as td

    seg, one = [], []
    for name, up in td.RUNS:
        d = td.RES[name].get("beta_de1", {}).get(up, [0])[0]
        if d < 20:
            continue
        T = td.RES[name]["age"][up][0]
        seg.append((T, td._ours(name, up, True)["beta_de1"] / d))
        one.append((T, td._ours(name, up, False)["beta_de1"] / d))
    fig, ax = plt.subplots(figsize=(3.6, 3.0), facecolor=SURFACE)
    _style(ax)
    ax.axhline(1, color=INK_3, lw=0.8)
    ax.scatter(*zip(*seg), s=16, color=CATEGORICAL[0], label="per U-series segment (EPRdating default)", lw=0, zorder=3)
    ax.scatter(*zip(*one), s=16, marker="s", facecolor="none", edgecolor=CATEGORICAL[1], lw=1,
               label="one factor for the chain (as DATA)", zorder=3)
    _logticks(ax, [5, 10, 20, 50, 100, 200, 500, 2000])
    ax.set_ylim(0.95, 1.6)
    ax.set_xlabel("DATA age (ka)")
    ax.set_ylabel("Dentine beta dose, EPRdating / DATA")
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    _save(fig, out)


# --- Fig. 5: US-ESR, DATA's p against the measured ratios ------------------------------
def fig_useries(out: Path):
    import test_data_useries as tu
    from published_cases import _th, _v, load_studies

    pts = []  # (measured, predicted from the program's own (T, p), tissue, source, flagged)
    for name, r in tu.RES.items():
        u = r.get("usesr", {})
        if "age" not in u:
            continue
        T = u["age"][0]
        for tissue, key in (("enamel", "p_internal"), ("de1", "p_beta_de1")):
            if key in u:
                (r48, _), (th, _) = tu.CASES[name]["useries"][tissue]
                pred = th230_u234(T, u[key][0], r48)
                pts.append((th, pred, tissue, "DATA, synthetic", name in tu.DATA_UNCONVERGED and tissue == "de1"))
    ks = load_studies()["duval2019_khok_sung"]
    for s in ks["samples"]:
        pub = s["published"]
        a = (pub.get("ages_ka") or {}).get("US")
        if not a or not pub.get("p_dentine"):
            continue
        for tissue, key in (("enamel", "p_enamel"), ("dentine", "p_dentine")):
            t = s[tissue]
            pred = th230_u234(a[0], pub[key][0], _v(t["u234_u238"]))
            meas = _v(_th(t))
            pts.append((meas, pred, tissue, "Khok Sung (Duval et al. 2019)", abs(pred - meas) > 0.02))
    fig, ax = plt.subplots(figsize=(3.6, 3.2), facecolor=SURFACE)
    _style(ax)
    ax.plot([0, 1], [0, 1], color=INK_3, lw=0.8)
    for src, mk in (("DATA, synthetic", "o"), ("Khok Sung (Duval et al. 2019)", "D")):
        for tissue, c in (("enamel", CATEGORICAL[0]), ("dentine", CATEGORICAL[1])):
            sel = [p for p in pts if p[3] == src and (p[2] == tissue or (tissue == "dentine" and p[2] == "de1"))]
            ok = [p for p in sel if not p[4]]
            bad = [p for p in sel if p[4]]
            lab = f"{tissue}, {'synthetic' if src.startswith('DATA') else 'Khok Sung'}"
            if ok:
                ax.scatter([p[0] for p in ok], [p[1] for p in ok], s=15, marker=mk, color=c, lw=0, label=lab, zorder=3)
            if bad:
                ax.scatter([p[0] for p in bad], [p[1] for p in bad], s=22, marker=mk, facecolor="none", edgecolor=c,
                           lw=1.1, zorder=3, label=f"{lab}, not converged")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Measured $^{230}$Th/$^{234}$U")
    ax.set_ylabel("$^{230}$Th/$^{234}$U from the published (age, p)")
    ax.legend(frameon=False, loc="upper left", fontsize=7)
    fig.tight_layout()
    _save(fig, out)


def fig_published(out: Path):
    import plot_published_ages as pp

    pp.main(out.with_suffix(".png"))
    plt.close("all")
    import matplotlib.pyplot as _plt  # same figure as PDF

    _plt.rcParams["pdf.fonttype"] = 42
    pp.plt.savefig = _plt.savefig
    pp.main(out)


def main():
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path)
    ap.add_argument("--m18", type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(ROOT / "tools"))
    if a.m18:
        drc = fig_m18(a.m18, a.out / "fig2_m18.pdf")
        print(f"  M18 linear fit, all points: De = {drc.De:.1f} ± {drc.De_sigma:.1f} Gy, chi2_red = {drc.chi2_red:.1f}")
    fig_rosy(a.out / "fig3_rosy.pdf")
    fig_beta(a.out / "fig4_beta.pdf")
    fig_useries(a.out / "fig5_useries.pdf")
    fig_published(a.out / "fig6_published.pdf")


if __name__ == "__main__":
    main()
