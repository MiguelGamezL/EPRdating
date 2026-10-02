import math
from itertools import pairwise

import pytest

from eprdating import Sediment, available_factor_sets, conversion_factors, matrix_dose_rates, water_correction
from eprdating.dose_rate import cosmic_dose_rate, cosmic_dose_rate_sea_level, geomagnetic_latitude


def test_factor_sets_load():
    assert {"adamiec_aitken_1998", "guerin_2011", "liritzis_2013"} <= set(available_factor_sets())
    cf = conversion_factors("guerin_2011")
    assert cf.get("U", "alpha").value == pytest.approx(2.795)
    assert cf.get("K", "beta").value == pytest.approx(0.7982)
    assert cf.get("K", "alpha").value == 0.0
    with pytest.raises(KeyError):
        conversion_factors("does-not-exist")


def test_matrix_dose_rates_linear():
    cf = conversion_factors("adamiec_aitken_1998")
    r = matrix_dose_rates(U=1, Th=0, K=0, factors=cf)
    assert r["gamma"] == pytest.approx(0.113)
    r2 = matrix_dose_rates(U=2, Th=10, K=1, factors=cf)
    assert r2["beta"] == pytest.approx(2 * 0.146 + 10 * 0.0273 + 0.782)


def test_water_correction():
    assert water_correction(1.0, 0.0, "beta") == 1.0
    assert water_correction(1.0, 0.2, "gamma") == pytest.approx(1 / 1.228)
    with pytest.raises(ValueError):
        water_correction(1.0, -0.1, "alpha")


def test_sediment_helper():
    s = Sediment(U=(2.0, 0.1), Th=(8.0, 0.4), K=(1.2, 0.05), water=0.15)
    dry = matrix_dose_rates(2.0, 8.0, 1.2)["gamma"]
    assert s.dose_rate("gamma") == pytest.approx(dry / (1 + 1.14 * 0.15))


def test_cosmic_continuity_and_sanity():
    # two branches of the Prescott & Hutton parametrisation meet near 1.67 hg/cm2
    below = cosmic_dose_rate_sea_level(1.669, 1.0)
    above = cosmic_dose_rate_sea_level(1.671, 1.0)
    assert below == pytest.approx(above, rel=0.02)
    # monotonic attenuation with depth
    d = [cosmic_dose_rate_sea_level(x, 2.0) for x in (0.5, 1, 2, 5, 10)]
    assert all(a > b for a, b in pairwise(d))
    # surface value at mid latitude, sea level: ~0.2-0.3 Gy/ka
    v = cosmic_dose_rate(1.0, 1.8, 45.0, 10.0, 0.0)
    assert 0.1 < v.value < 0.3
    assert v.sigma == pytest.approx(0.1 * v.value)


def test_altitude_increases_cosmic():
    low = cosmic_dose_rate(2.0, 1.8, 4.6, -74.1, 0.0).value
    high = cosmic_dose_rate(2.0, 1.8, 4.6, -74.1, 2600.0).value  # Bogotá altitude
    assert 1.3 * low < high < 2.0 * low


def test_geomagnetic_latitude_pole():
    assert geomagnetic_latitude(78.3, 291.0) > 88.0  # dipole pole (coefficients rounded in the source)
    assert -90 <= geomagnetic_latitude(4.6, -74.1) <= 90
    assert math.isfinite(geomagnetic_latitude(0, 0))
