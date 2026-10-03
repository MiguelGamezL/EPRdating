import pytest

from eprdating import Sediment, ToothLayers, ToothSample, matrix_dose_rates
from eprdating.dose_rate import u_series_split


def test_split_fractions():
    for rad in ("alpha", "beta", "gamma"):
        pre, post = u_series_split(rad)
        assert pre + post == pytest.approx(1.0)
    assert u_series_split("gamma")[1] > 0.95  # gamma comes from 214Pb/214Bi
    assert 0.5 < u_series_split("beta")[1] < 0.65  # 234mPa carries ~40 % of the beta


def test_matrix_rates_reduce_to_equilibrium():
    eq = matrix_dose_rates(2.0, 8.0, 3.0)
    assert matrix_dose_rates(2.0, 8.0, 3.0, U_ra226=2.0) == pytest.approx(eq)
    no_ra = matrix_dose_rates(2.0, 0.0, 0.0, U_ra226=0.0)
    only_u = matrix_dose_rates(2.0, 0.0, 0.0)
    for rad in ("alpha", "beta", "gamma"):
        assert no_ra[rad] == pytest.approx(only_u[rad] * u_series_split(rad)[0])


def _tooth(sed):
    return ToothSample(De=5.0, enamel_U=0.0, dentine_U=0.0, sediment=sed, beta=ToothLayers(1000), cosmic=0.18)


def test_tooth_sediment_beta_with_disequilibrium():
    def comp(t):
        return {c.name: c.rate for c in t.components()}

    base = comp(_tooth(Sediment(U=2.0, Th=8.0, K=3.0, water=0.1)))
    same = comp(_tooth(Sediment(U=2.0, Th=8.0, K=3.0, water=0.1, U_ra226=2.0)))
    assert same == pytest.approx(base)
    low = comp(_tooth(Sediment(U=2.0, Th=8.0, K=3.0, water=0.1, U_ra226=1.0)))
    assert low["gamma"] < base["gamma"] and low["sediment beta"] < base["sediment beta"]
    # the drop of the U part follows the post-radium share
    u_only = comp(_tooth(Sediment(U=2.0, water=0.1)))
    u_half = comp(_tooth(Sediment(U=2.0, water=0.1, U_ra226=1.0)))
    assert u_half["gamma"] / u_only["gamma"] == pytest.approx(1 - 0.5 * u_series_split("gamma")[1])
    # through 1 mm of enamel 234mPa (the hardest beta, before 226Ra) is attenuated
    # least, so the post-radium share of the beta reaching the enamel is lower
    r = u_half["sediment beta"] / u_only["sediment beta"]
    assert 1 - 0.5 * u_series_split("beta")[1] < r < 0.85
    _tooth(Sediment(U=(2.0, 0.5), water=0.1, U_ra226=(1.0, 0.05))).age_mc(n=20, seed=1)
