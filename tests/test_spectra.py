import numpy as np
import pytest

from eprdating import fit_dose_response
from eprdating.spectra import (
    ComponentBasis,
    component_vs_dose,
    double_integral,
    field_for_g,
    g_for_field,
    peak_to_peak,
    t1_b2_amplitude,
)


def lorentz_derivative(B, B0, w):
    """First derivative of a Lorentzian absorption with peak-to-peak width w."""
    g = w * np.sqrt(3) / 2
    x = (B - B0) / g
    return -2 * x / (1 + x**2) ** 2


def test_g_field_roundtrip():
    B = field_for_g(2.0023, 9.5)
    assert 338 < B < 340
    assert g_for_field(B, 9.5) == pytest.approx(2.0023)


def test_peak_to_peak_and_t1b2():
    f = 9.5
    B = np.linspace(336, 342, 6001)
    b0 = field_for_g(2.0, f)
    y = lorentz_derivative(B, b0, 0.2)
    assert peak_to_peak(B, y) == pytest.approx(y.max() - y.min())
    # put a positive lobe at T1 and a negative lobe at B2 by construction
    y2 = 0.5 * lorentz_derivative(B, field_for_g(2.0018, f) + 0.05, 0.1) \
        + 0.5 * lorentz_derivative(B, field_for_g(1.9973, f) - 0.05, 0.1)
    assert t1_b2_amplitude(B, y2, f) > 0


def test_double_integral_scales_linearly():
    B = np.linspace(330, 350, 4001)
    y = lorentz_derivative(B, 340, 0.5)
    assert double_integral(B, 3 * y) == pytest.approx(3 * double_integral(B, y), rel=1e-6)


def test_deconvolution_recovers_amplitudes_and_De():
    B = np.linspace(336, 342, 2001)
    s1 = lorentz_derivative(B, 338.9, 0.25)  # "dating" component
    s2 = lorentz_derivative(B, 339.4, 0.6)  # interfering component
    basis = ComponentBasis(B, {"dating": s1, "other": s2}, baseline_order=1)
    rng = np.random.default_rng(0)
    doses = np.array([0, 100, 200, 400, 800, 1600, 3200])
    De_true, D0 = 120.0, 1500.0
    amp_true = 1000 * (1 - np.exp(-(doses + De_true) / D0))
    spectra = []
    for a in amp_true:
        y = a * s1 / np.ptp(s1) + 300 * s2 / np.ptp(s2) + 0.02 * (B - 339) + 2 * rng.standard_normal(B.size)
        spectra.append(y)
    amps, results = component_vs_dose(basis, spectra, "dating")
    assert np.allclose(amps, amp_true, rtol=0.01)
    assert all(r.r2 > 0.99 for r in results)
    fit = fit_dose_response(doses, amps, "SSE")
    assert fit.De == pytest.approx(De_true, rel=0.05)


def test_basis_interpolates_tuple_shapes():
    B = np.linspace(336, 342, 500)
    Bf = np.linspace(335, 343, 3000)
    basis = ComponentBasis(B, {"a": (Bf, lorentz_derivative(Bf, 339, 0.3))})
    assert basis.matrix.shape == (500, 1)
    assert np.ptp(basis.matrix[:, 0]) == pytest.approx(1.0, rel=1e-2)
