import math

import pytest

from eprdating import (
    BetaGeometry,
    DoseRateComponent,
    EarlyUptake,
    LinearUptake,
    Sediment,
    ToothSample,
    USeries,
    USModel,
    solve_age,
)
from eprdating.series import activity_ratio_Th230_U238

# A synthetic partition (NOT published values) used only to test the algebra.
FAKE_PARTITION = {
    rad: {"U238": 0.2, "U234": 0.2, "Th230": 0.55, "U235": 0.05} for rad in ("alpha", "beta", "gamma")
}


def test_constant_dose_rate():
    r = solve_age(100.0, [DoseRateComponent("gamma", 0.6), DoseRateComponent("cosmic", 0.4)])
    assert r.age == pytest.approx(100.0)
    assert r.mean_dose_rate == pytest.approx(1.0)


@pytest.mark.parametrize("p,factor", [(-1, 1.0), (0, 0.5), (1, 1 / 3), (2.5, 1 / 4.5)])
def test_us_model_analytic(p, factor):
    # internal-only sample: De = rate * T / (p+2)
    comp = DoseRateComponent("internal", 2.0, USModel(p))
    T = solve_age(200.0, [comp]).age
    assert 2.0 * T * factor == pytest.approx(200.0, rel=1e-8)


def test_quadrature_matches_closed_form_in_equilibrium():
    G = lambda tau: tau
    for p in (-0.5, 0.0, 0.7, 3.0):
        assert USModel(p).accumulated(150.0, G) == pytest.approx(150.0 / (p + 2), rel=1e-6)


def test_ingrowth_limits():
    assert activity_ratio_Th230_U238(0.0) == 0.0
    assert activity_ratio_Th230_U238(5000.0) == pytest.approx(1.0, abs=1e-6)
    s = USeries(r0=1.0, partition=FAKE_PARTITION)
    g = s.G("alpha")
    # young: only U segments contribute (45 %); old: tends to equilibrium
    assert g(0.01) / 0.01 == pytest.approx(0.45, rel=1e-3)
    assert g(5000.0) / 5000.0 == pytest.approx(1.0, rel=0.05)


def test_ingrowth_makes_older_ages():
    eq = DoseRateComponent("int", 1.0, EarlyUptake())
    dis = DoseRateComponent("int", 1.0, EarlyUptake(), USeries(1.0, FAKE_PARTITION).G("alpha"))
    assert solve_age(50.0, [dis]).age > solve_age(50.0, [eq]).age


def test_partition_validation():
    with pytest.raises(ValueError):
        USeries(partition={r: {"U238": 0.5, "U234": 0.2, "Th230": 0.2, "U235": 0.0} for r in ("alpha", "beta", "gamma")})
    with pytest.raises(ValueError, match="not filled"):
        USeries()  # bundled table is still a template


def _tooth(**kw):
    base = dict(  # noqa: C408
        De=(200.0, 10.0),
        enamel_U=(0.5, 0.05),
        dentine_U=(10.0, 1.0),
        sediment=Sediment(U=(2.0, 0.1), Th=(6.0, 0.3), K=(1.0, 0.05), water=(0.1, 0.03)),
        beta=BetaGeometry(internal=(0.6, 0.03), dentine=(0.3, 0.02), external=(0.25, 0.02)),
        cosmic=(0.15, 0.015),
    )
    base.update(kw)
    return ToothSample(**base)


def test_tooth_sample_eu_lu_ordering():
    with pytest.warns(UserWarning, match="equilibrium"):
        eu = _tooth(uptake_enamel=EarlyUptake(), uptake_dentine=EarlyUptake()).age()
    with pytest.warns(UserWarning):
        lu = _tooth(uptake_enamel=LinearUptake(), uptake_dentine=LinearUptake()).age()
    assert lu.age > eu.age > 0
    assert sum(lu.accumulated.values()) == pytest.approx(200.0)


def test_tooth_mc():
    mc = _tooth().age_mc(n=400, seed=1)
    assert mc.samples.size > 390
    assert mc.mean == pytest.approx(mc.nominal.age, rel=0.05)
    lo, hi = mc.interval()
    assert lo < mc.nominal.age < hi
    assert math.isfinite(mc.std)
