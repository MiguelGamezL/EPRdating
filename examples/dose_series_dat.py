"""Additive-dose series from .dat/.par spectra: intensities, DRC and De.

Usage::

    python examples/dose_series_dat.py FOLDER DESIGN.csv [--window-G 100] [--center-g 2.0023]
                                       [--bridge NARROW WIDE] [--no-repeats]
                                       [--exclude-undetected ALPHA] [--epraya] [--plot out.png]

``DESIGN.csv`` lists the spectra of the series, one per row::

    file,dose,aliquot,exclude
    S_N.dat,0,,
    S_50Gy.dat,50,a50,
    S_50Gy_repeat.dat,50,a50,
    S_400Gy.dat,400,,tube cracked during the measurement

``aliquot`` (optional) groups the spectra of one aliquot: they are averaged,
the first one listed is the main spectrum. ``exclude`` (optional) gives the
reason to leave an aliquot out of the fits; it is still measured and drawn.
When some aliquots (e.g. the natural) were measured on a wider sweep,
``--bridge NARROW WIDE`` names one aliquot measured both ways, which puts
the wide sweeps on the scale of the narrow ones. The data are not part of
the repository.

Steps
-----
1. intensity window: 100 G around g = 2.0023 by default (``--window-G``,
   ``--center-g``); the sweep outside it gives the noise;
2. line-shape template: the field-aligned average of the strongest spectra
   or, with ``--epraya``, the orthorhombic CO2- radical simulated with EPRAYA
   and broadened (linewidth, modulation, lock-in time constant, field offset
   fitted to that average);
3. each aliquot: its repeated spectra are averaged (normalised to one power,
   weighted by their noise) and the template amplitude is fitted inside the
   window (linear baseline, common field shift within ±0.6 mT), with errors
   by noise injection; the repeats are checked against each other; a
   detection test gives the probability that noise alone produces the
   amplitude (``p_noise``). ``--exclude-undetected 0.01`` leaves out the
   aliquots above that level, for comparison only: dropping undetected
   aliquots biases De upwards (see ``tools/synthetic_dose_series.py``);
4. with ``--bridge``, wide-sweep aliquots are put on the scale of the narrow
   sweeps with an aliquot measured both ways;
5. linear dose-response fits (De allowed to be negative, errors inflated by
   the Birge ratio) for several point selections.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import warnings
from pathlib import Path

import numpy as np

from eprdating import fit_dose_response
from eprdating.spectra import (
    DEFAULT_WINDOW,
    IntensityWindow,
    aligned_average,
    combine_spectra,
    combined_intensity,
    intensity,
    pseudo_modulation,
    read_epr,
    time_constant_filter,
)

WIDE_SWEEP_KEEP = (305.0, 370.0)  # mT kept from the wide sweeps
SCALE = 1e8  # intensities in readable units
MAX_SHIFT = 0.6  # mT


class Design:
    """The series as read from ``DESIGN.csv``.

    ``series``: main file of each aliquot -> added dose (Gy); ``repeats``:
    further spectra of that aliquot; ``excluded``: main file -> reason."""

    def __init__(self, path: Path, bridge: tuple[str, str] | None = None):
        groups: dict[str, list[tuple[str, float]]] = {}
        self.excluded: dict[str, str] = {}
        with open(path, newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                name = row["file"].strip()
                if not name:
                    continue
                key = (row.get("aliquot") or "").strip() or name
                groups.setdefault(key, []).append((name, float(row["dose"])))
                if (row.get("exclude") or "").strip():
                    self.excluded[groups[key][0][0]] = row["exclude"].strip()
        self.series = {g[0][0]: g[0][1] for g in groups.values()}
        self.repeats = {g[0][0]: [n for n, _ in g[1:]] for g in groups.values() if len(g) > 1}
        self.bridge = bridge

    def strong(self) -> list[str]:
        """Main files of the upper half of the added doses (template)."""
        doses = [d for f, d in self.series.items() if d > 0 and f not in self.excluded]
        cut = float(np.median(doses)) if doses else 0.0
        return [f for f, d in self.series.items() if d >= cut and d > 0 and f not in self.excluded]


def load(folder: Path, name: str):
    s = read_epr(folder / name)
    return s.window(*WIDE_SWEEP_KEEP) if s.B.size > 1000 else s


def empirical_template(spectra, window: IntensityWindow, ref_power_mW: float | None = None):
    """Field-aligned average of strong spectra (same sweep), baseline removed."""
    comb = [combine_spectra([s], window, ref_power_mW=ref_power_mW).spectrum for s in spectra]
    B = comb[0].B
    m = window.mask(B, comb[0].freq_GHz)
    ys = []
    for c in comb:
        p = np.polyfit(B[~m], c.y[~m], 3)
        ys.append(c.y - np.polyval(p, B))
    avg, _ = aligned_average(B, ys, (B[m].min(), B[m].max()), max_shift=MAX_SHIFT)
    return B, avg


def epraya_template(B_ref, avg, window: IntensityWindow, freq_GHz: float):
    """Orthorhombic CO2- (Callens et al.) simulated with EPRAYA; broadening
    and field offset fitted to the empirical average."""
    from eprdating.spectra.epraya_backend import Species, simulate_species

    center = window.center(freq_GHz)
    m = window.mask(B_ref, freq_GHz)
    best = None
    for hpp in (0.2, 0.3, 0.45):
        Bs, ys = simulate_species(Species(g=[2.0031, 1.9973, 2.0019], Hpp=[0, hpp]), freq_GHz,
                                  (center - 10, center + 10), points=4001, grid=60)
        for ma, tau, sh in itertools.product((0.0, 0.2, 0.4), (2, 3, 4, 5, 6), np.arange(-0.8, 0.81, 0.05)):
            t = time_constant_filter(np.interp(B_ref, Bs + sh, pseudo_modulation(Bs, ys, ma)), tau)[m]
            A = np.c_[t, np.ones_like(t), B_ref[m] - center]
            x, *_ = np.linalg.lstsq(A, avg[m], rcond=None)
            chi = float(np.sum((avg[m] - A @ x) ** 2))
            if best is None or chi < best[0]:
                best = (chi, hpp, ma, tau, sh, Bs, ys)
    chi, hpp, ma, tau, sh, Bs, ys = best
    r2 = 1 - chi / np.sum((avg[m] - avg[m].mean()) ** 2)
    print(f"EPRAYA template: Hpp = {hpp} mT, modulation = {ma} mT, time constant = {tau} points, "
          f"field offset = {sh:+.2f} mT, R² vs average = {r2:.3f}")

    def make(B):  # the time constant acts per point, so apply it on each grid
        return time_constant_filter(np.interp(B, Bs + sh, pseudo_modulation(Bs, ys, ma)), tau)
    return make


def measure_series(folder: Path, design: Design, window: IntensityWindow = DEFAULT_WINDOW,
                   repeats: bool = True, epraya: bool = False):
    """Intensity of every aliquot. Returns ``(rows, shown)``: rows are
    ``(file, dose, value, sigma, Intensity)`` in ``SCALE`` units; ``shown``
    maps each file to ``(B, y, fitted)`` inside the window, for plotting."""
    groups = {f: [f, *(design.repeats.get(f, []) if repeats else [])] for f in design.series}
    spectra = {n: load(folder, n) for g in groups.values() for n in g}
    strong = [spectra[f] for f in design.strong()]
    ref_power = strong[0].power_mW
    B_ref, avg = empirical_template(strong, window, ref_power)
    template = epraya_template(B_ref, avg, window, strong[0].freq_GHz) if epraya else (B_ref, avg)
    kw = {"template": template, "window": window, "max_shift": MAX_SHIFT, "ref_power_mW": ref_power}

    res, shown = {}, {}
    for f, names in groups.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            res[f] = combined_intensity([spectra[n] for n in names], **kw)
            c = combine_spectra([spectra[n] for n in names], window, ref_power_mW=ref_power).spectrum
        m = window.mask(c.B, c.freq_GHz)
        shown[f] = (c.B[m], c.y[m] * SCALE, res[f].fit.fitted * SCALE)
    print(f"intensity window: {window.describe(strong[0].freq_GHz)}")

    k, sk, step0 = 1.0, 0.0, None
    if design.bridge:
        narrow, wide_name = design.bridge
        wide = intensity(load(folder, wide_name), **kw)
        k = res[narrow].value / wide.value
        sk = k * np.hypot(res[narrow].sigma / res[narrow].value, wide.sigma / wide.value)
        step0 = float(np.median(np.diff(spectra[narrow].B)))
        print(f"narrow/wide sweep factor = {k:.3f} ± {sk:.3f}")
    print()

    rows = []
    for f, d in design.series.items():
        r = res[f]
        A, s = r.value * SCALE, r.sigma * SCALE
        if step0 and float(np.median(np.diff(spectra[f].B))) > 1.1 * step0:  # wide -> narrow scale
            A, s = A * k, np.hypot(s * k, A * sk)  # right side uses the unscaled A
            B, y, fit = shown[f]
            shown[f] = (B, y * k, fit * k)
        rows.append((f, d, A, s, r))
        note = f"  excluded: {design.excluded[f]}" if f in design.excluded else ""
        rep = (f"  ({r.n_repeats} spectra, chi2_red = {r.chi2_red:.2f})" if r.n_repeats > 1 else "")
        print(f"{f:26s} {d:7.1f} Gy  I = {A:6.2f} ± {s:5.2f}  shift {r.shift_mT:+.2f} mT  "
              f"p_noise {r.p_noise:.3f}{rep}{note}")
    return rows, shown


def fit_selections(rows, excluded=(), alpha: float | None = None):
    """Linear fits for several selections. With ``alpha``, aliquots whose
    signal is not detected at that false-alarm level are left out too."""
    ok = [r for r in rows if r[0] not in excluded and (alpha is None or r[4].detected(alpha) is not False)]
    half = float(np.median([r[1] for r in ok])) if ok else 0.0
    sel = {
        "all aliquots": ok,
        "irradiated only": [r for r in ok if r[1] > 0],
        f"up to {half:g} Gy": [r for r in ok if r[1] <= half],
    }
    sel = {k: v for k, v in sel.items() if len(v) >= 3}
    return {name: fit_dose_response([r[1] for r in rr], [r[2] for r in rr], "LIN", sigma=[r[3] for r in rr],
                                    De_min=-np.inf) for name, rr in sel.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", type=Path)
    ap.add_argument("design", type=Path, help="CSV with columns file, dose[, aliquot][, exclude]")
    ap.add_argument("--window-G", type=float, default=100.0, help="intensity window width (G)")
    ap.add_argument("--center-g", type=float, default=2.0023, help="g-value at the window centre")
    ap.add_argument("--bridge", nargs=2, metavar=("NARROW", "WIDE"),
                    help="one aliquot measured on the narrow and on the wide sweep")
    ap.add_argument("--no-repeats", action="store_true", help="use only the main spectrum of each aliquot")
    ap.add_argument("--exclude-undetected", type=float, metavar="ALPHA",
                    help="leave out aliquots whose signal is not detected at this false-alarm level")
    ap.add_argument("--epraya", action="store_true", help="simulated template (needs EPRAYA)")
    ap.add_argument("--plot", type=Path, help="save a figure")
    a = ap.parse_args()

    window = IntensityWindow(a.window_G, "G", center_g=a.center_g)
    design = Design(a.design, tuple(a.bridge) if a.bridge else None)
    rows, shown = measure_series(a.folder, design, window, repeats=not a.no_repeats, epraya=a.epraya)
    print()
    if a.exclude_undetected is not None:
        out = [f"{r[1]:.0f} Gy" for r in rows if r[4].detected(a.exclude_undetected) is False]
        print(f"not detected at {a.exclude_undetected:g}, left out: {', '.join(out) or 'none'}")
    fits = fit_selections(rows, design.excluded, a.exclude_undetected)
    for name, f in fits.items():
        print(f"{name:24s} De = {f.De:6.1f} ± {f.De_sigma:5.1f} Gy   chi2_red = {f.chi2_red:.1f}")

    if a.plot:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        gap = 1.2 * max(np.ptp(shown[f][1]) for f, *_ in rows)
        for i, (f, d, *_rest) in enumerate(rows):
            B, y, fit = shown[f]
            ax1.plot(B, y - y.mean() + gap * i, color="0.7", lw=0.8)
            ax1.plot(B, fit - y.mean() + gap * i, color="C3" if f in design.excluded else "C0", lw=1.3)
            ax1.text(B.max() + 0.2, gap * i, f"{d:g} Gy", va="center", fontsize=8)
        ax1.set_xlabel("B (mT)")
        ax1.set_yticks([])
        ax1.set_title("spectra and template fits")
        for f, d, A, s, _ in rows:
            ax2.errorbar(d, A, s, fmt="o", color="C3" if f in design.excluded else "k", ms=4)
        dmax = max(r[1] for r in rows)
        x = np.linspace(-0.35 * dmax, 1.1 * dmax, 10)
        for (name, fit), c in zip(fits.items(), ("C0", "C1", "C2"), strict=False):
            ax2.plot(x, fit.predict(x), color=c, label=f"{name}: De = {fit.De:.0f} ± {fit.De_sigma:.0f} Gy")
        ax2.axhline(0, color="0.5", lw=0.5)
        ax2.set_xlabel("added dose (Gy)")
        ax2.set_ylabel("CO2- amplitude (a.u.)")
        ax2.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=150)


if __name__ == "__main__":
    main()
