"""U, Th and K of a sediment from HPGe spectra, and its dose rate.

Usage::

    python examples/norm_gamma.py FOLDER --sample corte0 --mass 500 [--water 0.10]

FOLDER holds ``U_24h.txt``, ``Th_24h.txt``, ``K_24h.txt`` (IAEA RGU-1, RGTh-1,
RGK-1, 500 g each), ``fondo_24h.txt`` (background), ``<sample>_24h.txt`` and
a calibration-source spectrum. The data are not part of the repository.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from eprdating import matrix_dose_rates, water_correction
from eprdating.gamma import (
    CALIBRATION_LINES,
    IAEA_RGK_1,
    IAEA_RGTH_1,
    IAEA_RGU_1,
    Calibration,
    analyse,
    calibrate_natural,
    read_spectrum_txt,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", type=Path)
    ap.add_argument("--sample", default="corte0")
    ap.add_argument("--mass", type=float, required=True, help="dry mass of the sample in the container (g)")
    ap.add_argument("--water", type=float, default=0.0, help="water content (water mass / dry mass)")
    ap.add_argument("--source", default="calibracion-57Co-22Na-137Cs-88Y_300s.txt")
    a = ap.parse_args()
    d = a.folder

    src = read_spectrum_txt(d / a.source)
    lines = [e for k in ("57Co", "22Na", "137Cs", "88Y") for e in CALIBRATION_LINES[k]]
    cal = src.calibrate(lines, guess=Calibration((1.86, 0.3517), (0.5, 0.0005)))
    print(f"source calibration: E = {cal.coef[0]:.3f} + {cal.coef[1]:.6f} n keV, "
          f"FWHM(1332 keV) = {2.355 * cal.sigma_keV(1332):.2f} keV")

    S = {k: read_spectrum_txt(d / f"{k}_24h.txt") for k in ("U", "Th", "K", "fondo", a.sample)}
    for s in S.values():
        c = calibrate_natural(s, cal, min_significance=5)
        print(f"  {s.name:14s} recalibrated on {len(c.residuals_keV)} natural lines "
              f"(1460.8 keV at channel {c.channel(1460.82):.1f})")

    refs = {"U": (S["U"], IAEA_RGU_1), "Th": (S["Th"], IAEA_RGTH_1), "K": (S["K"], IAEA_RGK_1)}
    r = analyse(S[a.sample], a.mass, refs, background=S["fondo"])
    print()
    print(r.summary())

    dry = matrix_dose_rates(r.U.value, r.Th.value, r.K.value)
    print(f"\ninfinite-matrix dose rates (Gy/ka), dry → water {a.water:.0%}:")
    for rad in ("alpha", "beta", "gamma"):
        print(f"  {rad:5s} {dry[rad]:.3f} → {water_correction(dry[rad], a.water, rad):.3f}")


if __name__ == "__main__":
    main()
