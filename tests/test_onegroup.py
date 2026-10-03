import math

import pytest

from eprdating.onegroup import (
    DENTINE,
    ELEMENTS,
    HYDROXYAPATITE,
    SILICA,
    Layer,
    ToothLayers,
    bethe_stopping,
    lewis_mus,
    mean_dose,
    solve_fluence,
)


def test_coefficients_match_prestwich_chan_2000():
    # Prestwich & Chan (2000), text after eq. 13: coefficients at 0.7 MeV (cm²/g)
    Z, A, I = ELEMENTS["H"]
    assert lewis_mus(Z, A, 0.7) == pytest.approx(1.790, abs=0.002)
    assert bethe_stopping(Z / A, I, 0.7) / 0.7 == pytest.approx(5.641, rel=0.005)
    Z, A, I = ELEMENTS["O"]
    assert lewis_mus(Z, A, 0.7) == pytest.approx(3.560, abs=0.002)
    assert bethe_stopping(Z / A, I, 0.7) / 0.7 == pytest.approx(2.401, rel=0.01)


def test_uniform_medium_gives_infinite_matrix_dose():
    layers = [Layer(SILICA, math.inf, 1.0), Layer(SILICA, 0.3, 1.0), Layer(SILICA, math.inf, 1.0)]
    sol = solve_fluence(layers, 0.8)
    assert mean_dose(sol, 1, 0.0, 0.3, 0.8) / 0.8 == pytest.approx(1.0, rel=1e-12)


@pytest.mark.parametrize("E", [0.05, 0.3, 0.8, 1.5])
def test_contributions_sum_to_one_in_a_single_material(E):
    geo = ToothLayers(1000, dentine_um=1e6, enamel=SILICA, dentine=SILICA, sediment=SILICA,
                      enamel_density=2.0, dentine_density=2.0, sediment_density=2.0)
    total = sum(geo.fraction(src, E) for src in ("enamel", "dentine", "sediment"))
    assert total == pytest.approx(1.0, abs=1e-9)


def test_materials_normalised():
    for m in (HYDROXYAPATITE, SILICA, DENTINE):
        assert sum(m.fractions.values()) == pytest.approx(1.0)


# Factors measured from ROSY 2.0 (tests/validation/ROSY_FINDINGS.md)
ROSY = {300: (0.391, 0.280, 0.303), 600: (0.563, 0.200, None), 1000: (0.692, 0.140, 0.152),
        1500: (0.781, 0.099, None), 3000: (0.888, 0.051, 0.056)}


@pytest.mark.parametrize("thick", sorted(ROSY))
def test_reproduces_rosy_beta_factors(thick):
    own, den, sed = ROSY[thick]
    geo = ToothLayers(thick)
    assert geo.chain_fraction("enamel", "U") == pytest.approx(own, rel=0.03)
    assert geo.chain_fraction("dentine", "U") == pytest.approx(den, rel=0.03)
    if sed is not None:
        assert geo.chain_fraction("sediment", "U") == pytest.approx(sed, rel=0.05)


def test_stripping_reduces_external_dose():
    a = ToothLayers(1000).chain_fraction("sediment", "U")
    b = ToothLayers(1000, strip_outer_um=50).chain_fraction("sediment", "U")
    c = ToothLayers(1000, strip_outer_um=100).chain_fraction("sediment", "U")
    assert a > b > c
    assert b == pytest.approx(0.136, rel=0.05) and c == pytest.approx(0.123, rel=0.05)  # ROSY


def test_segments_differ_in_self_absorption():
    geo = ToothLayers(1000)
    # 234mPa (U238 segment) is a harder emitter than the post-radon betas on average
    assert geo.chain_fraction("enamel", "U238") < geo.chain_fraction("enamel", "U234")


def test_geometry_with_uncertainties_and_validation():
    geo = ToothLayers(enamel_um=(1000, 100), strip_outer_um=(50, 20))
    assert geo.nominal_values()["enamel_um"] == 1000
    assert geo.at(enamel_um=900.0).chain_fraction("enamel", "U") < geo.at(enamel_um=1100.0).chain_fraction("enamel", "U")
    with pytest.raises(ValueError):
        geo.at(enamel_um=60.0, strip_outer_um=40.0, strip_inner_um=30.0)


def test_sediment_counts_both_sides():
    # a bare enamel fragment: sediment on both sides, equal contributions
    geo = ToothLayers(1000, dentine_um=0.0)
    outer, inner = geo.chain_fraction("sediment_outer", "U"), geo.chain_fraction("sediment_inner", "U")
    assert outer == pytest.approx(inner, rel=1e-6) and outer > 0
    assert geo.chain_fraction("sediment", "U") == pytest.approx(outer + inner, rel=1e-9)
    # thick dentine shields the inner side
    thick = ToothLayers(1000, dentine_um=5000.0)
    assert thick.chain_fraction("sediment_inner", "U") < 1e-4 * thick.chain_fraction("sediment_outer", "U")
