"""Synthetic additive-dose series that mimic a measured series, to test the analysis.

Each series has the added doses of the measured one; true intensity
a*(D + De) with the measured CO2- line shape (empirical template of the
strongest spectra), the measured noise (same power spectrum: time-constant
correlation and slow baseline wander), a random field position per aliquot
(sd 0.25 mT), and optionally a random multiplicative factor per aliquot
(inter-aliquot scatter: mass, filling, orientation).

The natural is measured either with one scan (noise 2.5x that of the
irradiated 4-scan spectra) or as well as the others.

Analysis as in examples/dose_series_dat.py: template fit inside the default
window (100 G around g = 2.0023) with a free shift (+-0.6 mT), noise-injection
errors, detection test, linear dose response (De may be negative), keeping
all aliquots, leaving out the natural, or leaving out the aliquots whose
signal is not detected at 1 %.

    python tools/synthetic_dose_series.py FOLDER DESIGN.csv [N_SERIES]
    python tools/synthetic_dose_series.py FOLDER DESIGN.csv [N_SERIES] --strategies

``DESIGN.csv`` as in examples/dose_series_dat.py. ``--strategies`` compares
measurement designs (more scans, less inter-aliquot scatter, more aliquots or
naturals, a wider dose range, a single aliquot) by the spread of the De they
give for a true De of 30 Gy.

The measured spectra provide the line shape and the noise.
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
ARGS = [a for a in sys.argv[1:] if a != "--strategies"]
if len(ARGS) < 2:
    sys.exit(__doc__)
FOLDER, DESIGN = Path(ARGS[0]), ds.Design(Path(ARGS[1]))
SCALE = 1e8
DOSES = np.array(sorted(DESIGN.series.values()))
SHIFT_SD = 0.25  # mT

# --- measured ingredients ----------------------------------------------------------
strong = [ds.load(FOLDER, f) for f in DESIGN.strong()]
B_t, y_t = ds.empirical_template(strong, DEFAULT_WINDOW)
y_t = y_t / np.ptp(y_t)
ref = strong[0]
B, F, POWER = ref.B, ref.freq_GHz, ref.power_mW
GAIN = ref.gain or 1.0  # intensities are per unit gain
m = DEFAULT_WINDOW.mask(B, F)
# slope of the measured series (intensity per Gy, x1e-8)
_I = [intensity(s, template=(B_t, y_t), max_shift=0.6, n_noise=5, n_null=5).value * SCALE / GAIN for s in strong]
SLOPE = float(np.polyfit([DESIGN.series[f] for f in DESIGN.strong()], _I, 1)[0]) if len(_I) > 1 else 0.01


def noise_psd(files):
    """Average power spectrum of the signal-free sides of the 4-scan spectra."""
    ps = []
    for f in files:
        s = ds.load(FOLDER, f)
        y = s.y * SCALE
        lo, hi = B[m].min(), B[m].max()
        for side in (y[s.B < lo], y[s.B > hi]):
            side = side - np.polyval(np.polyfit(np.arange(side.size), side, 3), np.arange(side.size))
            fr, p = welch(side, nperseg=128, detrend=False)
            ps.append(p)
    return fr, np.mean(ps, axis=0)


FR, PSD = noise_psd([f for f in DESIGN.series if DESIGN.series[f] > 0 and f not in DESIGN.excluded
                     and ds.load(FOLDER, f).B.size == B.size])


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
        out.append(Spectrum(B=B, scans=y[None, :] / SCALE, name=f"{D:g} Gy", freq_GHz=F, power_mW=POWER))
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
    nat_gt_first, nat_detected, nat_bias = 0, 0, []
    for _ in range(n):
        v, e, p = analyse(series(rng, De, cv, natural_noise))
        nat_gt_first += v[0] > v[1]
        nat_detected += p[0] < 0.01
        nat_bias.append(v[0] / (SLOPE * De) - 1)
        rows["all"].append(de(v, e, np.ones(DOSES.size, bool)))
        rows["no natural"].append(de(v, e, DOSES > 0))
        rows["detected 1%"].append(de(v, e, p < 0.01))
    return rows, nat_gt_first / n, nat_detected / n, np.median(nat_bias)


def design(rng, De, doses, cv, noise, nat_noise):
    """Series with any doses; ``noise`` scales the noise of the irradiated
    aliquots (0.5 = 16 scans against 4), ``nat_noise`` that of the naturals."""
    out = []
    for D in doses:
        A = GAIN * SLOPE * (D + De) * (1 + cv * rng.standard_normal())
        shape = np.interp(B, B_t + rng.normal(0, SHIFT_SD), y_t)
        y = 3.6e2 + A * shape + make_noise(rng, factor=nat_noise if D == 0 else noise)
        out.append(Spectrum(B=B, scans=y[None, :] / SCALE, freq_GHz=F, power_mW=POWER))
    return out


def study(name, doses, cv, noise=1.0, nat_noise=2.5, De=30.0, n=150, seed=1):
    rng = np.random.default_rng(seed)
    doses = np.asarray(doses, float)
    est, err = [], []
    for _ in range(n):
        res = [intensity(s, template=(B_t, y_t), max_shift=0.6, n_noise=100, n_null=5)
               for s in design(rng, De, doses, cv, noise, nat_noise)]
        v = np.array([r.value for r in res]) * SCALE / GAIN
        e = np.array([r.sigma for r in res]) * SCALE / GAIN
        f = fit_dose_response(doses, v, "LIN", sigma=e, De_min=-np.inf)
        est.append(f.De)
        err.append(f.De_sigma)
    q16, q50, q84 = np.quantile(est, [0.16, 0.5, 0.84])
    print(f"{name:52s} De {q50:5.1f}  spread ±{(q84 - q16) / 2:5.1f} Gy  "
          f"typical error {np.median(err):5.1f} Gy", flush=True)


def strategies(n):
    d = DOSES
    step = float(np.median(np.diff(np.unique(d))))
    print(f"true De = 30 Gy, {n} series per design; spread = half the 16-84 % range of the De\n")
    study("as measured (30 % scatter, natural 1 scan)", d, 0.30, n=n)
    study("natural with 4 scans", d, 0.30, nat_noise=1.0, n=n)
    study("3 naturals", np.r_[0, 0, 0, d[d > 0]], 0.30, n=n)
    study("16 scans for every aliquot", d, 0.30, noise=0.5, nat_noise=0.5, n=n)
    study("twice the aliquots, half the step", np.arange(0, d.max() + step / 4, step / 2), 0.30, n=n)
    study("twice the dose range", 2 * d, 0.30, n=n)
    study("scatter 15 %", d, 0.15, n=n)
    study("scatter 5 %", d, 0.05, n=n)
    study("single aliquot (3 %), 4 scans", d, 0.03, nat_noise=1.0, n=n)
    study("scatter 15 % + 16 scans + 3 naturals", np.r_[0, 0, 0, d[d > 0]], 0.15, noise=0.5, nat_noise=0.5, n=n)
    study("scatter 5 % + 16 scans", d, 0.05, noise=0.5, nat_noise=0.5, n=n)
    study("single aliquot (3 %) + 16 scans", d, 0.03, noise=0.5, nat_noise=0.5, n=n)


if __name__ == "__main__":
    if "--strategies" in sys.argv:
        strategies(int(ARGS[2]) if len(ARGS) > 2 else 150)
        sys.exit()
    n = int(ARGS[2]) if len(ARGS) > 2 else 200
    print(f"{n} synthetic series per case; De estimates: median, 16-84 % range, "
          "coverage of the 1-sigma error, median chi2_red\n")
    for De in (20.0, 50.0):
        for cv in (0.0, 0.2):
            for nn, label in ((2.5, "natural with 1 scan"), (1.0, "natural with 4 scans")):
                rows, gt, det, nb = run(De, cv, nn, n, seed=int(De * 10 + cv * 100 + nn * 7))
                print(f"De = {De:.0f} Gy, inter-aliquot scatter {cv:.0%}, {label}: "
                      f"natural > first irradiated in {gt:.0%} of series, natural detected in {det:.0%}, "
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
