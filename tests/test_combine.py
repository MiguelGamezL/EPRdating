"""Averaging repeated spectra of the same aliquot."""

import numpy as np
import pytest
from test_intensity_window import SHAPE, B, F

from eprdating.spectra import (
    DEFAULT_WINDOW,
    Intensity,
    Spectrum,
    combine_intensities,
    combine_spectra,
    combined_intensity,
    intensity,
    time_constant_filter,
)


def repeat(A, nscans, rng, sigma=1.0, shift=0.0, **kw):
    sh = np.interp(B, B + shift, SHAPE)
    rows = [A * sh + 0.5 + time_constant_filter(sigma * 2.0 * rng.standard_normal(B.size), 2)
            for _ in range(nscans)]
    return Spectrum(B=B, scans=np.array(rows), name=kw.pop("name", "s"), freq_GHz=kw.pop("freq_GHz", F),
                    power_mW=kw.pop("power_mW", 19.0), **kw)


def test_mean_not_sum_and_noise_weights():
    rng = np.random.default_rng(1)
    one, four = repeat(10.0, 1, rng, name="1scan"), repeat(10.0, 4, rng, name="4scan")
    c = combine_spectra([one, four])
    assert c.members == ["1scan#0", "4scan#0", "4scan#1", "4scan#2", "4scan#3"]
    assert c.weights.sum() == pytest.approx(1.0)
    assert np.all(np.abs(c.weights - 0.2) < 0.08)  # five scans of similar noise
    t = (B, SHAPE)
    combined = intensity(c.spectrum, template=t)
    assert abs(combined.value - 10.0) < 3 * combined.sigma  # the mean keeps the scale
    # noise of the average falls as 1/sqrt(5)
    out = (B < DEFAULT_WINDOW.bounds(F)[0] - 1)
    ratio = np.std(np.diff(c.spectrum.y[out])) / np.std(np.diff(one.y[out]))
    assert ratio == pytest.approx(1 / np.sqrt(5), rel=0.25)


def test_noisy_scans_weigh_less():
    rng = np.random.default_rng(2)
    c = combine_spectra([repeat(4.0, 1, rng, sigma=1.0, name="good"), repeat(4.0, 1, rng, sigma=3.0, name="bad")])
    assert c.weights[0] / c.weights[1] == pytest.approx(9.0, rel=0.3)
    assert np.all(combine_spectra([repeat(4, 1, rng), repeat(4, 1, rng, sigma=3)], weighting="equal").weights == 0.5)


def test_power_frequency_and_alignment():
    rng = np.random.default_rng(3)
    a = repeat(20.0, 4, rng, sigma=0.3)
    # same aliquot at a quarter of the power: sqrt(1/4) of the signal
    b = repeat(10.0, 4, rng, sigma=0.15, power_mW=4.75)
    c = combine_spectra([a, b])
    assert intensity(c.spectrum, template=(B, SHAPE)).value == pytest.approx(20.0, rel=0.03)
    # a file recorded at a higher frequency sits at a higher field
    hi = Spectrum(B=B, scans=np.array([20 * np.interp(B, B * 9.40 / F, SHAPE)]), freq_GHz=9.40, power_mW=19.0)
    lo = Spectrum(B=B, scans=np.array([20 * SHAPE]), freq_GHz=F, power_mW=19.0)
    assert intensity(combine_spectra([lo, hi]).spectrum, template=(B, SHAPE), max_shift=0).value == pytest.approx(
        20.0, rel=0.01)
    # a tuning shift without a frequency record: only the alignment removes it
    files = [repeat(30.0, 4, rng, sigma=0.3, shift=s) for s in (-0.3, 0.0, 0.3)]
    al = combine_spectra(files, align=True)
    assert al.shifts_mT[0] - al.shifts_mT[1] == pytest.approx(0.3, abs=0.03)
    assert al.shifts_mT[2] - al.shifts_mT[1] == pytest.approx(-0.3, abs=0.03)
    aligned = intensity(al.spectrum, template=(B, SHAPE)).value
    smeared = intensity(combine_spectra(files).spectrum, template=(B, SHAPE)).value
    assert aligned == pytest.approx(30.0, rel=0.04) and smeared < 0.88 * 30.0


def test_different_sweeps_are_not_averaged():
    rng = np.random.default_rng(4)
    a = repeat(4.0, 1, rng)
    wide = Spectrum(B=np.linspace(86, 586, 4096), scans=np.zeros((1, 4096)), freq_GHz=F, power_mW=19.0)
    with pytest.raises(ValueError, match="combine_intensities"):
        combine_spectra([a, wide])
    other_tc = repeat(4.0, 1, rng, time_constant_ms=40.0)
    with pytest.raises(ValueError, match="time constant"):
        combine_spectra([repeat(4.0, 1, rng, time_constant_ms=10.0), other_tc])


def test_rejects_a_scan_with_a_spike():
    rng = np.random.default_rng(5)
    s = repeat(4.0, 4, rng, name="x")
    s.scans[2, 255:258] += 40.0  # a spike inside the window
    assert combine_spectra([s]).rejected == []
    c = combine_spectra([s], reject_chi2=2.0)
    assert c.rejected == ["x#2"] and c.weights[2] == 0


def test_combined_intensity_gains_signal_to_noise():
    rng = np.random.default_rng(6)
    files = [repeat(3.0, 1, rng, name="a"), repeat(3.0, 4, rng, name="b"), repeat(3.0, 4, rng, name="c")]
    r = combined_intensity(files, template=(B, SHAPE), n_noise=100)
    assert r.n_repeats == 3 and len(r.repeats) == 3 and r.chi2_red is not None
    single = intensity(files[0], template=(B, SHAPE), n_noise=100)
    assert r.sigma < 0.55 * single.sigma  # nine scans against one
    assert r.value == pytest.approx(3.0, abs=3 * r.sigma)


def test_repeat_scatter_inflates_the_error():
    rng = np.random.default_rng(7)
    # strong signal, but the aliquot gives 20 % different amplitudes when repositioned
    files = [repeat(A, 4, rng, sigma=0.2, name=str(A)) for A in (40.0, 48.0, 36.0)]
    r = combined_intensity(files, template=(B, SHAPE), n_noise=60)
    assert r.chi2_red > 10
    plain = intensity(combine_spectra(files).spectrum, template=(B, SHAPE), n_noise=60)
    assert r.sigma == pytest.approx(plain.sigma * np.sqrt(r.chi2_red), rel=1e-6)


def test_combine_intensities():
    a = Intensity(10.0, 1.0, "template", (331.0, 341.0))
    b = Intensity(12.0, 1.0, "template", (331.0, 341.0))
    c = combine_intensities([a, b])
    assert c.value == pytest.approx(11.0)
    assert c.chi2_red == pytest.approx(2.0)
    assert c.sigma == pytest.approx(np.sqrt(0.5) * np.sqrt(2.0))
    assert combine_intensities([a]) is a
    with pytest.raises(ValueError):
        combine_intensities([a, Intensity(1.0, 0.1, "peak_to_peak", (331.0, 341.0))])
