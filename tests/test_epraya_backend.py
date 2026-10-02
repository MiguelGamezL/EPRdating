"""Integration test with the real EPRAYA backend (skipped if unavailable)."""

import numpy as np
import pytest

pytestmark = pytest.mark.epraya


def test_simulated_basis_deconvolution():
    try:
        import epraya  # noqa: F401

        from eprdating.spectra.epraya_backend import Species, basis_from_species, simulate_species
    except ImportError as exc:
        pytest.skip(f"EPRAYA backend not available: {exc}")

    freq = 9.5
    B = np.linspace(336.5, 341.5, 1024)
    species = {
        "orthorhombic": Species(g=[2.0031, 1.9973, 2.0019], Hpp=[0, 0.25]),
        "isotropic": Species(g=2.0045, Hpp=[0, 0.6]),
    }
    basis = basis_from_species(B, species, freq, grid=40)
    Bo, yo = simulate_species(species["orthorhombic"], freq, (B[0], B[-1]), points=1024, grid=40)
    Bi, yi = simulate_species(species["isotropic"], freq, (B[0], B[-1]), points=1024, grid=40)
    y = 5 * np.interp(B, Bo, yo) / np.ptp(yo) + 2 * np.interp(B, Bi, yi) / np.ptp(yi)
    res = basis.fit(y)
    assert res.amplitudes["orthorhombic"] == pytest.approx(5, rel=0.02)
    assert res.amplitudes["isotropic"] == pytest.approx(2, rel=0.02)
