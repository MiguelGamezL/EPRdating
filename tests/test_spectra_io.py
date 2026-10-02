import numpy as np
import pytest
from scipy.signal import lfilter

from eprdating import fit_dose_response
from eprdating.spectra import (
    ComponentBasis,
    aligned_average,
    normalise,
    pseudo_modulation,
    read_dat,
    subtract_baseline,
    time_constant_filter,
)


def lorentz_derivative(B, B0, dBpp):
    x = (B - B0) / (dBpp * np.sqrt(3) / 2)
    return -x / (1 + x**2) ** 2


def _write(tmp_path, stem, scans, B_G, par=None, separator=True, trailing=0):
    lines, k = [], 0
    for s, y in enumerate(scans):
        if s and separator:
            lines.append(f"{k} 0 {B_G[0]:.6f} NaN 1")
        for i, v in enumerate(y):
            lines.append(f"{k} {i} {B_G[i]:.6f} {v:.6E} 1")
            k += 1
    for i in range(trailing):
        lines.append(f"{k} {i} {B_G[i]:.6f} {scans[0][i]:.6E} 1")
        k += 1
    (tmp_path / f"{stem}.dat").write_text("\r\n".join(lines) + "\r\n")
    p = {"N": len(B_G), "CF": 3350.0, "CF_": 3354.0, "SW": 499.1453, "Nscans": len(scans), "Freq": 9.43,
         "TC": 15, "MA": 40, "RG": 30}
    p.update(par or {})
    (tmp_path / f"{stem}.par").write_text("\r\n".join(f"{k} : {v}" for k, v in p.items()))
    return tmp_path / f"{stem}.dat"


def test_read_dat_par(tmp_path):
    B = np.linspace(3100, 3600, 64)
    rng = np.random.default_rng(0)
    scans = rng.normal(size=(3, 64))
    f = _write(tmp_path, "S_7_19mW_4SCAN", scans, B, trailing=5)
    s = read_dat(f)
    assert s.n_scans == 3 and s.dropped_points == 5
    # G -> mT, moved by CF_ - CF = +4 G to the actual field
    assert s.B[0] == pytest.approx(310.4) and s.B[-1] == pytest.approx(360.4)
    assert read_dat(f, actual_field=False).B[0] == pytest.approx(310.0)
    np.testing.assert_allclose(s.y, scans.mean(0), rtol=1e-5)
    assert s.freq_GHz == 9.43 and s.power_mW == 19.0 and s.gain == 30
    w = s.window(330, 340)
    assert w.B.min() >= 330 and w.B.max() <= 340 and w.scans.shape[0] == 3


def test_read_dat_without_separator_lines(tmp_path):
    B = np.linspace(3100, 3600, 32)
    f = _write(tmp_path, "x_18.5mW", np.ones((2, 32)), B, separator=False)
    s = read_dat(f)
    assert s.n_scans == 2 and s.power_mW == 18.5


def test_subtract_baseline_and_normalise():
    B = np.linspace(330, 345, 300)
    sig = lorentz_derivative(B, 337.5, 0.5)
    y = sig + 3 + 0.2 * (B - 337) + 0.01 * (B - 337) ** 3
    out = subtract_baseline(B, y, exclude=(335.5, 339.5), order=3)
    assert np.max(np.abs(out - sig)) < 0.02
    assert normalise(np.ones(3), power_mW=4.0, ref_power_mW=1.0, gain=2.0, mass_mg=10.0)[0] == pytest.approx(
        1 / 2 / 2 / 10)


def test_pseudo_modulation_limits_and_lorentzian_optimum():
    B = np.linspace(320, 360, 8001)
    dBpp = 0.5
    y = lorentz_derivative(B, 340, dBpp)
    np.testing.assert_allclose(pseudo_modulation(B, y, 1e-4), y, atol=1e-6)
    # the recorded signal (∝ Bm × derivative-normalised output) peaks at
    # Bm ≈ 3.5 ΔBpp for a Lorentzian (Poole, Electron Spin Resonance)
    ratios = np.linspace(1, 8, 57)
    amp = [r * np.ptp(pseudo_modulation(B, y, r * dBpp)) for r in ratios]
    assert ratios[int(np.argmax(amp))] == pytest.approx(3.5, abs=0.3)
    # overmodulation broadens the line
    yb = pseudo_modulation(B, y, 5 * dBpp)
    assert B[np.argmin(yb)] - B[np.argmax(yb)] > 2 * dBpp


def test_time_constant_filter_delays_line():
    B = np.linspace(330, 350, 2001)
    y = lorentz_derivative(B, 340, 0.5)
    yf = time_constant_filter(y, 20)
    c = lambda s: 0.5 * (B[np.argmax(s)] + B[np.argmin(s)])
    assert c(yf) > c(y)


def test_aligned_average_recovers_shifts():
    B = np.linspace(330, 345, 600)
    step = B[1] - B[0]
    shifts = np.array([0, 3, -2, 5]) * step
    ys = [lorentz_derivative(B, 337 + s, 0.4) for s in shifts]
    avg, found = aligned_average(B, ys, (335, 339), max_shift=0.3, reference=0)
    np.testing.assert_allclose(found, shifts, atol=1e-9)
    np.testing.assert_allclose(avg, ys[0], atol=1e-3)


def test_fit_shift_and_noise_injection_errors():
    B = np.linspace(332.7, 340.7, 82)
    t = lorentz_derivative(B, 336.7, 0.8)
    basis = ComponentBasis(B, {"s": t}, baseline_order=1)
    rng = np.random.default_rng(3)
    a = np.exp(-1 / 3)

    def noise(n):
        return lfilter([1 - a], [1, -a], rng.standard_normal(n + 50))[50:]

    amps, errs = [], []
    for _ in range(40):
        y = 30 * lorentz_derivative(B, 336.9, 0.8) / np.ptp(t) + 3 * noise(B.size)
        r = basis.fit(y, nonnegative=False, max_shift=0.5, noise=3 * noise(400), n_noise=80)
        amps.append(r.amplitudes["s"])
        errs.append(r.errors["s"])
        assert abs(r.shift - 0.2) < 0.15
    assert np.mean(amps) == pytest.approx(30, rel=0.05)
    assert np.mean(errs) == pytest.approx(np.std(amps), rel=0.4)


def test_negative_De_and_birge_scaling():
    D = np.array([0, 20, 40, 60, 80.0])
    I = 2.0 * (D - 10)  # extrapolates to De = -10
    f = fit_dose_response(D, I, "LIN", sigma=np.ones(5), De_min=-np.inf)
    assert f.De == pytest.approx(-10, abs=1e-6)
    noisy = I + np.array([3, -3, 3, -3, 3.0])
    a = fit_dose_response(D, noisy, "LIN", sigma=np.ones(5), De_min=-np.inf)
    b = fit_dose_response(D, noisy, "LIN", sigma=np.ones(5), De_min=-np.inf, scale_errors=False)
    assert a.chi2_red > 1
    assert a.De_sigma == pytest.approx(b.De_sigma * np.sqrt(a.chi2_red))
