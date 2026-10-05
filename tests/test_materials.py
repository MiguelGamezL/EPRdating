import pytest

from eprdating import Material, ToothLayers, compound, dentine_material, mixture, sediment_material
from eprdating.alpha import alpha_range
from eprdating.onegroup import CALCITE, DENTINE, HYDROXYAPATITE, SILICA


def test_compound_and_mixture():
    calcite = compound("calcite", {"Ca": 1, "C": 1, "O": 3})
    assert calcite.fractions["Ca"] == pytest.approx(40.078 / 100.086, rel=1e-4)
    assert sum(calcite.fractions.values()) == pytest.approx(1.0)
    m = mixture("half", [(SILICA, 1.0), (CALCITE, 1.0)])
    assert m.fractions["Si"] == pytest.approx(0.5 * SILICA.fractions["Si"])
    with pytest.raises(ValueError):
        compound("x", {"U": 1})
    with pytest.raises(ValueError):
        mixture("x", [(SILICA, 0.0)])


def test_helpers():
    assert dentine_material().key == DENTINE.key
    assert sediment_material().key == SILICA.key
    s = sediment_material(quartz=0.6, calcite=0.3, kaolinite=0.1)
    assert set(s.fractions) == {"Si", "O", "Ca", "C", "Al", "H"}


def test_same_name_different_composition_do_not_share_caches():
    a = Material("custom", {"Si": 1, "O": 2})
    b = Material("custom", {"Ca": 1, "C": 1, "O": 3})
    assert a.coefficients(0.5) != b.coefficients(0.5)
    assert alpha_range(5.3, a) != alpha_range(5.3, b)


def test_sediment_composition_changes_the_beta_factor_slightly():
    kw = dict(enamel_um=1000, dentine_um=5000, strip_outer_um=50, strip_inner_um=50)  # noqa: C408
    quartz = ToothLayers(**kw).chain_fraction("sediment", "U")
    calc = ToothLayers(**kw, sediment=CALCITE).chain_fraction("sediment", "U")
    mix = ToothLayers(**kw, sediment=sediment_material(quartz=0.5, calcite=0.5)).chain_fraction("sediment", "U")
    assert calc < mix < quartz
    assert calc == pytest.approx(quartz, rel=0.04)


def test_materials_survive_the_monte_carlo_copy():
    geo = ToothLayers(enamel_um=(1000, 100), sediment=CALCITE, dentine=dentine_material(0.6, 0.3, 0.1))
    draw = geo.at(enamel_um=950.0)
    assert draw.sediment.key == CALCITE.key and draw.dentine.key != DENTINE.key
    assert draw.enamel.key == HYDROXYAPATITE.key
