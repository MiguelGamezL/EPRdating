"""Additive-dose series from .dat/.par spectra: intensities, DRC and De.

Usage::

    python examples/dose_series_dat.py FOLDER [--epraya] [--plot out.png]

The series below is the M18 enamel sample (aliquots irradiated in 20 Gy
steps, 4 scans at 19 mW; the natural aliquot measured on a wide sweep, 1 scan
at 18 mW). Edit ``SERIES`` for other samples. The data are not part of the
repository.

Steps
-----
1. read every spectrum, average its scans and normalise to 19 mW (√P);
2. build a line-shape template: the field-aligned average of the strongest
   spectra or, with ``--epraya``, the orthorhombic CO2- radical simulated with
   EPRAYA and broadened (linewidth, lock-in time constant, field offset
   fitted to that average);
3. fit the template amplitude in each spectrum (±4 mT window, linear
   baseline, common field shift within ±0.6 mT), with errors by noise
   injection from the signal-free part of the same spectrum;
4. put the wide-sweep natural aliquot on the scale of the narrow sweeps with
   an aliquot measured both ways;
5. linear dose-response fits (De allowed to be negative, errors inflated by
   the Birge ratio) for several point selections.
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np

from eprdating import fit_dose_response
from eprdating.spectra import (
    ComponentBasis,
    aligned_average,
    normalise,
    pseudo_modulation,
    read_epr,
    subtract_baseline,
    time_constant_filter,
)

# file -> added dose (Gy)
SERIES = {"M18_0_18mW.dat": 0.0, **{f"M18_{i}_19mW_4SCAN.dat": 20.0 * i for i in range(1, 10)}}
EXCLUDED = {"M18_5_19mW_4SCAN.dat": "no detectable signal at 100 Gy (failed measurement)"}
# same aliquot measured with the narrow and the wide sweep (for the natural)
BRIDGE = ("M18_9_19mW_4SCAN.dat", "M18_9_19mW_TOTAL.dat")
REF_POWER = 19.0  # mW
CENTER, HALF = 337.1, 4.0  # mT (actual field), fit window
SIGNAL = (331.4, 342.9)  # mT, excluded from baseline/noise estimates
MAX_SHIFT = 0.6  # mT


def prepared(path: Path):
    s = read_epr(path)
    if s.B.size > 1000:  # wide sweep: keep the region of the narrow sweeps
        s = s.window(305.0, 370.0)
    y = normalise(s.y, power_mW=s.power_mW, ref_power_mW=REF_POWER) * 1e8
    yb = subtract_baseline(s.B, y, exclude=SIGNAL, order=3)
    left, right = yb[s.B < SIGNAL[0]], yb[s.B > SIGNAL[1]]
    return s.B, y, (left if left.size >= right.size else right)


def empirical_template(spectra):
    B = spectra[0][0]
    ys = [subtract_baseline(B, y, exclude=SIGNAL, order=3) for _, y, _ in spectra]
    avg, _ = aligned_average(B, ys, (CENTER - 3, CENTER + 3), max_shift=MAX_SHIFT)
    return B, avg


def epraya_template(B_ref, avg):
    """Orthorhombic CO2- (Callens et al.) simulated with EPRAYA; broadening
    and field offset fitted to the empirical average."""
    from eprdating.spectra.epraya_backend import Species, simulate_species

    m = np.abs(B_ref - CENTER) <= HALF
    best = None
    for hpp in (0.2, 0.3, 0.45):
        Bs, ys = simulate_species(Species(g=[2.0031, 1.9973, 2.0019], Hpp=[0, hpp]), 9.43,
                                  (CENTER - 10, CENTER + 10), points=4001, grid=60)
        for ma, tau, sh in itertools.product((0.0, 0.2, 0.4), (2, 3, 4, 5, 6), np.arange(-0.8, 0.81, 0.05)):
            t = time_constant_filter(np.interp(B_ref, Bs + sh, pseudo_modulation(Bs, ys, ma)), tau)[m]
            A = np.c_[t, np.ones_like(t), B_ref[m] - CENTER]
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


def amplitude(B, y, noise, template):
    m = np.abs(B - CENTER) <= HALF
    shape = template(B)[m] if callable(template) else template
    basis = ComponentBasis(B[m], {"CO2-": shape}, baseline_order=1)
    r = basis.fit(y[m], nonnegative=False, max_shift=MAX_SHIFT, noise=noise, n_noise=300)
    return r.amplitudes["CO2-"], r.errors["CO2-"], r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", type=Path)
    ap.add_argument("--epraya", action="store_true", help="simulated template (needs EPRAYA)")
    ap.add_argument("--plot", type=Path, help="save a figure")
    a = ap.parse_args()

    data = {f: prepared(a.folder / f) for f in [*SERIES, BRIDGE[1]]}
    strong = [data[f] for f, d in SERIES.items() if d >= 60 and f not in EXCLUDED]
    B_ref, avg = empirical_template(strong)
    template = epraya_template(B_ref, avg) if a.epraya else (B_ref, avg)

    res = {f: amplitude(*data[f], template) for f in data}
    k = res[BRIDGE[0]][0] / res[BRIDGE[1]][0]
    sk = k * np.hypot(res[BRIDGE[0]][1] / res[BRIDGE[0]][0], res[BRIDGE[1]][1] / res[BRIDGE[1]][0])
    print(f"narrow/wide sweep factor = {k:.3f} ± {sk:.3f}\n")

    step = lambda f: float(np.median(np.diff(data[f][0])))
    rows = []
    for f, d in SERIES.items():
        A, s, r = res[f]
        if step(f) > 1.1 * step(BRIDGE[0]):  # wide sweep -> narrow-sweep scale
            A, s = A * k, np.hypot(s * k, A * sk)
        rows.append((f, d, A, s, r.shift))
        note = f"  excluded: {EXCLUDED[f]}" if f in EXCLUDED else ""
        print(f"{f:26s} {d:5.0f} Gy  I = {A:6.2f} ± {s:5.2f}  shift {r.shift:+.2f} mT{note}")

    print()
    sel = {
        "all (natural included)": [r for r in rows if r[0] not in EXCLUDED],
        "irradiated only": [r for r in rows if r[0] not in EXCLUDED and r[1] > 0],
        "natural + 20-80 Gy": [r for r in rows if r[0] not in EXCLUDED and r[1] <= 80],
    }
    fits = {}
    for name, rr in sel.items():
        f = fit_dose_response([r[1] for r in rr], [r[2] for r in rr], "LIN", sigma=[r[3] for r in rr],
                              De_min=-np.inf)
        fits[name] = f
        print(f"{name:24s} De = {f.De:6.1f} ± {f.De_sigma:5.1f} Gy   chi2_red = {f.chi2_red:.1f}")

    if a.plot:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        for i, (f, d, *_rest) in enumerate(rows):
            B, y, _ = data[f]
            m = np.abs(B - CENTER) <= HALF
            r = res[f][2]
            ax1.plot(B[m], y[m] - y[m].mean() + 40 * i, color="0.7", lw=0.8)
            ax1.plot(B[m], r.fitted - y[m].mean() + 40 * i, color="C0" if f not in EXCLUDED else "C3", lw=1.3)
            ax1.text(CENTER + HALF + 0.2, 40 * i, f"{d:.0f} Gy", va="center", fontsize=8)
        ax1.set_xlabel("B (mT)")
        ax1.set_yticks([])
        ax1.set_title("spectra and template fits")
        for f, d, A, s, _ in rows:
            ax2.errorbar(d, A, s, fmt="o", color="C3" if f in EXCLUDED else "k", ms=4)
        x = np.linspace(-60, 200, 10)
        for (name, fit), c in zip(fits.items(), ("C0", "C1", "C2")):
            ax2.plot(x, fit.predict(x), color=c, label=f"{name}: De = {fit.De:.0f} ± {fit.De_sigma:.0f} Gy")
        ax2.axhline(0, color="0.5", lw=0.5)
        ax2.set_xlabel("added dose (Gy)")
        ax2.set_ylabel("CO2- amplitude (a.u.)")
        ax2.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=150)


if __name__ == "__main__":
    main()
