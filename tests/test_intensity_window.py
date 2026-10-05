"""Intensity window and intensity of a single spectrum."""

import warnings

import numpy as np
import pytest

from eprdating.spectra import (
    DEFAULT_WINDOW,
    IntensityWindow,
    Spectrum,
    double_integral,
    field_for_g,
    intensity,
    peak_to_peak,
    time_constant_filter,
)

B = np.linspace(310.39, 360.31, 512)  # the UNAL X-band sweep (50 mT, 512 points)
F = 9.43


def ld(x, x0, w):
    g = w * np.sqrt(3) / 2
    u = (x - x0) / g
    return -2 * u / (1 + u**2) ** 2


SHAPE = ld(B, field_for_g(2.0006, F), 0.6)
SHAPE /= np.ptp(SHAPE)


def spectrum(A, sigma=0.0, rng=None, extra=None, nscans=1, **kw):
    rng = rng or np.random.default_rng(0)
    rows = []
    for _ in range(nscans):
        y = A * SHAPE + 0.3 + 0.002 * (B - 335)
        if extra is not None:
            y = y + extra
        if sigma:
            y = y + time_constant_filter(sigma * 2.0 * rng.standard_normal(B.size), 2)
        rows.append(y)
    return Spectrum(B=B, scans=np.array(rows), name=kw.pop("name", "s"), freq_GHz=kw.pop("freq_GHz", F),
                    power_mW=kw.pop("power_mW", 19.0), **kw)


def test_default_window_is_100_G_around_g_2_0023():
    w = DEFAULT_WINDOW
    assert (w.width, w.unit, w.center_g) == (100.0, "G", 2.0023)
    assert w.half_width_mT == pytest.approx(5.0)
    assert w.center(F) == pytest.approx(field_for_g(2.0023, F))
    lo, hi = w.bounds(F)
    assert hi - lo == pytest.approx(10.0)
    assert "100 G around g = 2.0023" in w.describe(F)


def test_window_options_and_errors():
    assert IntensityWindow(8, unit="mT").half_width_mT == 4.0
    assert IntensityWindow(60).half_width_mT == pytest.approx(3.0)
    assert IntensityWindow(80, center_mT=337.1).bounds(None) == pytest.approx((333.1, 341.1))
    with pytest.raises(ValueError):
        IntensityWindow(100, unit="T")
    with pytest.raises(ValueError):
        IntensityWindow(0)
    with pytest.raises(ValueError):
        DEFAULT_WINDOW.center(None)  # a g-centred window needs the frequency
    with pytest.warns(UserWarning, match="beyond the sweep"):
        IntensityWindow(800).mask(B, F)
    # at Q band the same window follows the resonance field
    assert DEFAULT_WINDOW.center(34.0) == pytest.approx(field_for_g(2.0023, 34.0))


def test_window_keeps_other_signals_out():
    other = 3.0 * ld(B, field_for_g(2.0023, F) + 6.5, 0.5)  # a line 6.5 mT above the centre
    s = spectrum(1.0, extra=other)
    assert peak_to_peak(s.B, s.y, DEFAULT_WINDOW, F) == pytest.approx(1.0, rel=0.02)
    assert peak_to_peak(s.B, s.y, IntensityWindow(150), F) > 2.0  # a wide window takes it in
    # tuples in mT still work
    lo, hi = DEFAULT_WINDOW.bounds(F)
    assert peak_to_peak(s.B, s.y, (lo, hi)) == peak_to_peak(s.B, s.y, DEFAULT_WINDOW, F)
    m = DEFAULT_WINDOW.mask(s.B, F)
    assert double_integral(s.B, s.y, window=DEFAULT_WINDOW, freq_GHz=F) == double_integral(s.B[m], s.y[m])


def test_template_intensity_and_its_error():
    rng = np.random.default_rng(3)
    vals, errs = [], []
    for _ in range(150):
        r = intensity(spectrum(5.0, sigma=0.5, rng=rng), template=(B, SHAPE), n_noise=60)
        vals.append(r.value)
        errs.append(r.sigma)
    assert np.mean(vals) == pytest.approx(5.0, rel=0.03)
    # the noise-injection error matches the scatter over noise realisations
    assert np.mean(errs) == pytest.approx(np.std(vals, ddof=1), rel=0.2)
    assert r.method == "template" and r.window_mT[1] - r.window_mT[0] == pytest.approx(10.0, abs=0.2)


def test_other_methods_and_options():
    s = spectrum(5.0, sigma=0.3)
    for method in ("peak_to_peak", "t1_b2", "double_integral"):
        r = intensity(s, method)
        assert np.isfinite(r.value) and r.sigma > 0, method
    with pytest.raises(ValueError):
        intensity(s)  # the template method needs a template
    with pytest.raises(ValueError):
        intensity(s, "area")
    # power normalisation (square-root law)
    low = spectrum(5.0 * np.sqrt(4.75 / 19.0), power_mW=4.75)
    assert intensity(low, template=(B, SHAPE), ref_power_mW=19.0).value == pytest.approx(5.0, rel=1e-3)
    # a narrow sweep leaves no signal-free stretch as long as the window
    narrow = Spectrum(B=B[200:330], scans=spectrum(5.0, sigma=0.3).scans[:, 200:330], freq_GHz=F)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = intensity(narrow, "peak_to_peak")
    assert np.isnan(r.sigma)


def test_detection_test_false_alarm_rate():
    rng = np.random.default_rng(11)
    p_noise = [intensity(spectrum(0.0, sigma=0.5, rng=rng), template=(B, SHAPE), n_noise=10, n_null=150,
                         seed=k).p_noise for k in range(40)]
    # with no signal, p_noise is roughly uniform: few "detections" at 5 %
    # (100 such spectra gave 5 % below 0.05 and 2 % below 0.01)
    assert np.mean(np.array(p_noise) < 0.05) < 0.15
    strong = intensity(spectrum(5.0, sigma=0.5, rng=rng), template=(B, SHAPE), n_noise=10, n_null=200)
    assert strong.p_noise < 0.01 and strong.detected() is True
    assert intensity(spectrum(5.0, sigma=0.5, rng=rng), "peak_to_peak").detected() is None
