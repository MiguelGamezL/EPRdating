"""Synthetic additive-dose series that mimic M18, to test the analysis.

Each series: 10 aliquots (0-180 Gy, 20 Gy steps), true intensity
a*(D + De) with the real CO2- line shape (empirical template of M18), the
real noise of M18 (same power spectrum: time-constant correlation and slow
baseline wander), a random field position per aliquot (sd 0.25 mT, as seen
in M18), and optionally a random multiplicative factor per aliquot
(inter-aliquot scatter: mass, filling, orientation).

The natural is measured either as in M18 (one scan, noise 2.5x that of the
4-scan spectra) or as well as the others (4 scans).

Analysis as in examples/dose_series_dat.py: template fit inside the default
window (100 G around g = 2.0023) with a free shift (+-0.6 mT), noise-injection
errors, detection test, linear dose response (De may be negative), keeping
all aliquots, leaving out the natural, or leaving out the aliquots whose
signal is not detected at 1 %.

    python tools/synthetic_dose_series.py M18_FOLDER [N_SERIES]

The M18 spectra (not in the repository) provide the line shape and the noise.
"""

import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.signal import welch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples"))
import dose_series_dat as ds

from eprdating import fit_dose_response
from eprdating.spectra import DEFAULT_WINDOW, Spectrum, intensity

warnings.simplefilter("ignore")
M18 = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("M18")
SCALE = 1e8
DOSES = np.arange(0, 181, 20.0)
SLOPE = 0.010  # intensity per Gy (x1e-8), M18-like
SHIFT_SD = 0.25  # mT
GAIN = 30.0  # receiver gain of the M18 spectra (intensities are per unit gain)

# --- real ingredients from M18 ---------------------------------------------------
strong = [ds.load(M18, f) for f, d in ds.SERIES.items() if d >= 60 and f not in ds.EXCLUDED]
B_t, y_t = ds.empirical_template(strong, DEFAULT_WINDOW)
y_t = y_t / np.ptp(y_t)
ref = ds.load(M18, "M18_7_19mW_4SCAN.dat")
B, F = ref.B, ref.freq_GHz
m = DEFAULT_WINDOW.mask(B, F)


def noise_psd(files):
    """Average power spectrum of the signal-free sides of the 4-scan spectra."""
    ps = []
    for f in files:
        s = ds.load(M18, f)
        y = s.y * SCALE
        lo, hi = B[m].min(), B[m].max()
        for side in (y[s.B < lo], y[s.B > hi]):
            side = side - np.polyval(np.polyfit(np.arange(side.size), side, 3), np.arange(side.size))
            fr, p = welch(side, nperseg=128, detrend=False)
            ps.append(p)
    return fr, np.mean(ps, axis=0)


FR, PSD = noise_psd([f for f in ds.SERIES if f.endswith("4SCAN.dat") and f not in ds.EXCLUDED])


def make_noise(rng, n=B.size, factor=1.0):
    f = np.fft.rfftfreq(n)
    amp = np.sqrt(np.interp(f, FR, PSD) * n / 2)
    z = amp * np.exp(2j * np.pi * rng.random(f.size))
    z[0] = 0
    return factor * np.fft.irfft(z, n=n)


def series(rng, De, cv, natural_noise):
    out = []
    for D in DOSES:
        A = GAIN * SLOPE * (D + De) * (1 + cv * rng.standard_normal())
        shape = np.interp(B, B_t + rng.normal(0, SHIFT_SD), y_t)
        y = 3.6e2 + A * shape + make_noise(rng, factor=natural_noise if D == 0 else 1.0)
        out.append(Spectrum(B=B, scans=y[None, :] / SCALE, name=f"{D:.0f} Gy", freq_GHz=F, power_mW=19.0))
    return out


def analyse(spectra):
    res = [intensity(s, template=(B_t, y_t), max_shift=0.6, n_noise=100, n_null=300) for s in spectra]
    v = np.array([r.value for r in res]) * SCALE / GAIN
    e = np.array([r.sigma for r in res]) * SCALE / GAIN
    p = np.array([r.p_noise for r in res])
    return v, e, p


def de(v, e, keep):
    if keep.sum() < 3:
        return np.nan, np.nan, np.nan
    f = fit_dose_response(DOSES[keep], v[keep], "LIN", sigma=e[keep], De_min=-np.inf)
    return f.De, f.De_sigma, f.chi2_red


def run(De, cv, natural_noise, n, seed):
    rng = np.random.default_rng(seed)
    rows = {k: [] for k in ("all", "no natural", "detected 1%")}
    nat_gt_20, nat_detected, nat_bias = 0, 0, []
    for _ in range(n):
        v, e, p = analyse(series(rng, De, cv, natural_noise))
        nat_gt_20 += v[0] > v[1]
        nat_detected += p[0] < 0.01
        nat_bias.append(v[0] / (SLOPE * De) - 1)
        rows["all"].append(de(v, e, np.ones(10, bool)))
        rows["no natural"].append(de(v, e, DOSES > 0))
        rows["detected 1%"].append(de(v, e, p < 0.01))
    return rows, nat_gt_20 / n, nat_detected / n, np.median(nat_bias)


if __name__ == "__main__":
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    print(f"{n} synthetic series per case; De estimates: median, 16-84 % range, "
          "coverage of the 1-sigma error, median chi2_red\n")
    for De in (20.0, 50.0):
        for cv in (0.0, 0.2):
            for nn, label in ((2.5, "natural as in M18 (1 scan)"), (1.0, "natural with 4 scans")):
                rows, gt, det, nb = run(De, cv, nn, n, seed=int(De * 10 + cv * 100 + nn * 7))
                print(f"De = {De:.0f} Gy, inter-aliquot scatter {cv:.0%}, {label}: "
                      f"natural > 20 Gy in {gt:.0%} of series, natural detected in {det:.0%}, "
                      f"natural intensity bias {nb:+.0%}")
                for k, r in rows.items():
                    r = np.array(r)
                    r = r[np.isfinite(r[:, 0])]
                    if not len(r):
                        print(f"    {k:12s} no fit")
                        continue
                    q16, q50, q84 = np.quantile(r[:, 0], [0.16, 0.5, 0.84])
                    cover = np.mean(np.abs(r[:, 0] - De) <= r[:, 1])
                    print(f"    {k:12s} De {q50:6.1f} Gy  [{q16:6.1f}, {q84:6.1f}]   coverage {cover:.0%}   "
                          f"chi2 {np.median(r[:, 2]):.1f}")
                print()
