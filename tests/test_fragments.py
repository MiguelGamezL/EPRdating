"""Enamel fragments measured at several angles (eprdating.spectra.fragments)."""

import numpy as np
import pytest

from eprdating import fit_dose_response
from eprdating.spectra import (
    AngularSpectrum,
    IntensityWindow,
    Spectrum,
    XrayCalibration,
    angular_profile,
    angular_set,
    field_for_g,
    fragment_intensity,
    merge_angular,
    parse_angular_name,
)

F = 9.39
B = np.linspace(313.0, 363.0, 2361)
WIN = IntensityWindow(40, center_g=2.000)


def line(Bx, g, amp, w=0.35):
    x = (Bx - field_for_g(g, F)) / (w * np.sqrt(3) / 2)
    return amp * -2 * x / (1 + x**2) ** 2 / 0.77


def fragment(amp, rng, orientations="XYZ", angles=range(0, 181, 20), noise=2.0, name="S"):
    """An anisotropic signal: position and height change with the angle and
    differ between orientations."""
    out = []
    for k, o in enumerate(orientations):
        for a in angles:
            t = np.radians(a + 30 * k)
            g = 2.0006 + 0.0012 * np.cos(2 * t)
            y = line(B, g, amp * (1 + 0.3 * np.sin(t) ** 2)) + noise * rng.standard_normal(B.size)
            out.append(Spectrum(B=B, scans=y[None, :], name=f"{name}_{o}_gon_{a}dg_result", freq_GHz=F))
    return out


def test_names():
    assert parse_angular_name("BF_ESR_01_360s_Z_gon_60dg_result") == ("Z", 60.0)
    assert parse_angular_name("AV_ESR_05_Frag_X_14400s_gon_0dg_result") == ("X", 0.0)
    assert parse_angular_name("BF_ESR_01_7202s_y_gon_120dg_result") == ("Y", 120.0)
    assert parse_angular_name("BF_ESR_06_Nat_Z_2_gon_100dg_result") == ("Z", 100.0)
    assert parse_angular_name("tooth_or-b_angle15", r"or-(?P<orientation>\w)_angle(?P<angle>\d+)") == ("B", 15.0)
    with pytest.raises(ValueError):
        parse_angular_name("M18_3_19mW_4SCAN")


def test_merge_balances_orientations():
    rng = np.random.default_rng(0)
    items = angular_set(fragment(100.0, rng))
    # drop most angles of one orientation: the balanced merge still weighs it as the others
    items = [a for a in items if a.orientation != "Y" or a.angle_deg < 40]
    m = merge_angular(items)
    assert m.n_scans == 3 and m.params["orientations"] == ["X", "Y", "Z"] and m.params["n_spectra"] == 22
    per = {o: np.mean([a.spectrum.y for a in items if a.orientation == o], axis=0) for o in "XYZ"}
    assert np.allclose(m.y, np.mean(list(per.values()), axis=0))
    flat = merge_angular(items, balance=False)
    assert flat.n_scans == 22 and not np.allclose(flat.y, m.y)
    prof = angular_profile(items, window=WIN)
    assert list(prof["X"][0]) == list(range(0, 181, 20)) and prof["Y"][0].size == 2
    assert np.ptp(prof["X"][1]) > 0.1 * prof["X"][1].mean()  # the fragment is anisotropic


def test_merge_puts_other_frequencies_on_the_g_scale():
    rng = np.random.default_rng(1)
    a = fragment(100.0, rng, "X", (0,), noise=0.0)[0]
    b = Spectrum(B=a.B * 9.45 / F, scans=a.scans, name="b_X_gon_20dg", freq_GHz=9.45)
    m = merge_angular([AngularSpectrum(a, "X", 0), AngularSpectrum(b, "X", 20)])
    # b recorded at another frequency lands on a's g scale: both angles give a's spectrum
    assert np.allclose(m.scans[0], np.interp(m.B, a.B, a.y), atol=1e-6 * np.ptp(a.y))


def test_fragment_series_recovers_De_in_seconds():
    rng = np.random.default_rng(2)
    times = np.array([0, 90, 360, 900, 1800, 3600, 7200, 14400.0])
    De_s, Dsat = 6000.0, 20000.0
    I, S = [], []
    for t in times:
        amp = 100 * (1 - np.exp(-(t + De_s) / Dsat))
        r = fragment_intensity(angular_set(fragment(amp, rng)), "peak_to_peak", window=WIN, n_noise=30)
        I.append(r.value)
        S.append(r.sigma)
    f = fit_dose_response(times, I, "SSE", sigma=S)
    assert abs(f.De - De_s) < 3 * f.De_sigma + 0.05 * De_s
    cal = XrayCalibration(0.235, 0.012)
    De, sD = cal.De(f.De, f.De_sigma)
    assert De == pytest.approx(0.235 * f.De)
    assert sD >= De * 0.012 / 0.235
    assert cal.dose([0, 100]).tolist() == pytest.approx([0.0, 23.5])
    with pytest.raises(ValueError):
        XrayCalibration(0.0)
    # with the orientation scatter the error grows (anisotropy)
    r2 = fragment_intensity(angular_set(fragment(50.0, rng)), "peak_to_peak", window=WIN, n_noise=30,
                            orientation_scatter=True)
    r1 = fragment_intensity(angular_set(fragment(50.0, rng)), "peak_to_peak", window=WIN, n_noise=30)
    assert r2.sigma > r1.sigma and r2.n_repeats == 3
