"""Additive-dose series from .dat/.par spectra: intensities, DRC and De.

Usage::

    python examples/dose_series_dat.py FOLDER [--window-G 100] [--center-g 2.0023]
                                              [--no-repeats] [--exclude-undetected ALPHA]
                                              [--epraya] [--plot out.png]

The series below is the M18 enamel sample (aliquots irradiated in 20 Gy
steps, 4 scans at 19 mW; the natural aliquot measured on a wide sweep, 1 scan
at 18 mW). Aliquots 1 and 9 were also measured with a single scan; those
repeats are averaged with the 4-scan spectra. Edit ``SERIES`` and
``REPEATS`` for other samples. The data are not part of the repository.

Steps
-----
1. intensity window: 100 G around g = 2.0023 by default (``--window-G``,
   ``--center-g``); the sweep outside it gives the noise;
2. line-shape template: the field-aligned average of the strongest spectra
   or, with ``--epraya``, the orthorhombic CO2- radical simulated with EPRAYA
   and broadened (linewidth, modulation, lock-in time constant, field offset
   fitted to that average);
3. each aliquot: its repeated spectra are averaged (normalised to 19 mW,
   weighted by their noise) and the template amplitude is fitted inside the
   window (linear baseline, common field shift within ±0.6 mT), with errors
   by noise injection; the repeats are checked against each other; a
   detection test gives the probability that noise alone produces the
   amplitude (``p_noise``). ``--exclude-undetected 0.01`` leaves out the
   aliquots above that level, for comparison only: dropping undetected
   aliquots biases De upwards (see ``tools/synthetic_dose_series.py``);
4. the wide-sweep natural aliquot is put on the scale of the narrow sweeps
   with an aliquot measured both ways;
5. linear dose-response fits (De allowed to be negative, errors inflated by
   the Birge ratio) for several point selections.
"""

from __future__ import annotations

import argparse
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

# main file of each aliquot -> added dose (Gy)
SERIES = {"M18_0_18mW.dat": 0.0, **{f"M18_{i}_19mW_4SCAN.dat": 20.0 * i for i in range(1, 10)}}
# further spectra of the same aliquot, same sweep (averaged with the main one)
REPEATS = {"M18_1_19mW_4SCAN.dat": ["M18_1_19mW.dat"], "M18_9_19mW_4SCAN.dat": ["M18_9_19mW.dat"]}
EXCLUDED = {"M18_5_19mW_4SCAN.dat": "no detectable signal at 100 Gy (failed measurement)"}
# same aliquot measured with the narrow and the wide sweep (for the natural)
BRIDGE = ("M18_9_19mW_4SCAN.dat", "M18_9_19mW_TOTAL.dat")
WIDE_SWEEP_KEEP = (305.0, 370.0)  # mT kept from the wide sweeps
REF_POWER = 19.0  # mW
SCALE = 1e8  # intensities in readable units
MAX_SHIFT = 0.6  # mT


def load(folder: Path, name: str):
    s = read_epr(folder / name)
    return s.window(*WIDE_SWEEP_KEEP) if s.B.size > 1000 else s


def empirical_template(spectra, window: IntensityWindow):
    """Field-aligned average of strong spectra (same sweep), baseline removed."""
    comb = [combine_spectra([s], window, ref_power_mW=REF_POWER).spectrum for s in spectra]
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


def measure_series(folder: Path, window: IntensityWindow = DEFAULT_WINDOW, repeats: bool = True,
                   epraya: bool = False):
    """Intensity of every aliquot. Returns ``(rows, shown)``: rows are
    ``(file, dose, value, sigma, Intensity)`` in ``SCALE`` units; ``shown``
    maps each file to ``(B, y, fitted)`` inside the window, for plotting."""
    groups = {f: [f, *(REPEATS.get(f, []) if repeats else [])] for f in SERIES}
    spectra = {n: load(folder, n) for g in groups.values() for n in g}
    spectra[BRIDGE[1]] = load(folder, BRIDGE[1])
    strong = [spectra[f] for f, d in SERIES.items() if d >= 60 and f not in EXCLUDED]
    B_ref, avg = empirical_template(strong, window)
    template = epraya_template(B_ref, avg, window, strong[0].freq_GHz) if epraya else (B_ref, avg)

    res, shown = {}, {}
    for f, names in groups.items():
        kw = {"template": template, "window": window, "max_shift": MAX_SHIFT}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            res[f] = combined_intensity([spectra[n] for n in names], ref_power_mW=REF_POWER, **kw)
            c = combine_spectra([spectra[n] for n in names], window, ref_power_mW=REF_POWER).spectrum
        m = window.mask(c.B, c.freq_GHz)
        shown[f] = (c.B[m], c.y[m] * SCALE, res[f].fit.fitted * SCALE)
    wide = intensity(spectra[BRIDGE[1]], template=template, window=window, ref_power_mW=REF_POWER,
                     max_shift=MAX_SHIFT)
    k = res[BRIDGE[0]].value / wide.value
    sk = k * np.hypot(res[BRIDGE[0]].sigma / res[BRIDGE[0]].value, wide.sigma / wide.value)
    print(f"intensity window: {window.describe(strong[0].freq_GHz)}")
    print(f"narrow/wide sweep factor = {k:.3f} ± {sk:.3f}\n")

    step0 = float(np.median(np.diff(spectra[BRIDGE[0]].B)))
    rows = []
    for f, d in SERIES.items():
        r = res[f]
        A, s = r.value * SCALE, r.sigma * SCALE
        if float(np.median(np.diff(spectra[f].B))) > 1.1 * step0:  # wide sweep -> narrow-sweep scale
            A, s = A * k, np.hypot(s * k, A * sk)  # right side uses the unscaled A
            B, y, fit = shown[f]
            shown[f] = (B, y * k, fit * k)
        rows.append((f, d, A, s, r))
        note = f"  excluded: {EXCLUDED[f]}" if f in EXCLUDED else ""
        rep = (f"  ({r.n_repeats} spectra, chi2_red = {r.chi2_red:.2f})" if r.n_repeats > 1 else "")
        print(f"{f:26s} {d:5.0f} Gy  I = {A:6.2f} ± {s:5.2f}  shift {r.shift_mT:+.2f} mT  "
              f"p_noise {r.p_noise:.3f}{rep}{note}")
    return rows, shown


def fit_selections(rows, alpha: float | None = None):
    """Linear fits for several selections. With ``alpha``, aliquots whose
    signal is not detected at that false-alarm level are left out too."""
    ok = [r for r in rows if r[0] not in EXCLUDED and (alpha is None or r[4].detected(alpha) is not False)]
    sel = {
        "all aliquots": ok,
        "irradiated only": [r for r in ok if r[1] > 0],
        "up to 80 Gy": [r for r in ok if r[1] <= 80],
    }
    sel = {k: v for k, v in sel.items() if len(v) >= 3}
    return {name: fit_dose_response([r[1] for r in rr], [r[2] for r in rr], "LIN", sigma=[r[3] for r in rr],
                                    De_min=-np.inf) for name, rr in sel.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", type=Path)
    ap.add_argument("--window-G", type=float, default=100.0, help="intensity window width (G)")
    ap.add_argument("--center-g", type=float, default=2.0023, help="g-value at the window centre")
    ap.add_argument("--no-repeats", action="store_true", help="use only the main spectrum of each aliquot")
    ap.add_argument("--exclude-undetected", type=float, metavar="ALPHA",
                    help="leave out aliquots whose signal is not detected at this false-alarm level")
    ap.add_argument("--epraya", action="store_true", help="simulated template (needs EPRAYA)")
    ap.add_argument("--plot", type=Path, help="save a figure")
    a = ap.parse_args()

    window = IntensityWindow(a.window_G, "G", center_g=a.center_g)
    rows, shown = measure_series(a.folder, window, repeats=not a.no_repeats, epraya=a.epraya)
    print()
    if a.exclude_undetected is not None:
        out = [f"{r[1]:.0f} Gy" for r in rows if r[4].detected(a.exclude_undetected) is False]
        print(f"not detected at {a.exclude_undetected:g}, left out: {', '.join(out) or 'none'}")
    fits = fit_selections(rows, a.exclude_undetected)
    for name, f in fits.items():
        print(f"{name:24s} De = {f.De:6.1f} ± {f.De_sigma:5.1f} Gy   chi2_red = {f.chi2_red:.1f}")

    if a.plot:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        for i, (f, d, *_rest) in enumerate(rows):
            B, y, fit = shown[f]
            ax1.plot(B, y - y.mean() + 40 * i, color="0.7", lw=0.8)
            ax1.plot(B, fit - y.mean() + 40 * i, color="C0" if f not in EXCLUDED else "C3", lw=1.3)
            ax1.text(B.max() + 0.2, 40 * i, f"{d:.0f} Gy", va="center", fontsize=8)
        ax1.set_xlabel("B (mT)")
        ax1.set_yticks([])
        ax1.set_title("spectra and template fits")
        for f, d, A, s, _ in rows:
            ax2.errorbar(d, A, s, fmt="o", color="C3" if f in EXCLUDED else "k", ms=4)
        x = np.linspace(-60, 200, 10)
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
