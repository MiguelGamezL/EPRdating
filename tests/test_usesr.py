import math

import pytest

from eprdating import Sediment, ToothLayers, ToothSample, USModel
from eprdating.series import LAMBDA
from eprdating.usesr import (
    UseriesData,
    USESRSample,
    closed_system_age,
    incoming_ratio,
    predicted_ratios,
    solve_p,
    th230_u234,
)

BASE = dict(  # noqa: C408
    enamel_U=2.0, dentine_U=40.0, sediment=Sediment(U=2, Th=6, K=1, water=0.1),
    beta=ToothLayers(1000), cosmic=0.15,
)


def _synthetic(T, pe, pd, rin_e, rin_d):
    r48e, th_e = predicted_ratios(T, pe, rin_e)
    r48d, th_d = predicted_ratios(T, pd, rin_d)
    fwd = ToothSample(De=1.0, uptake_enamel=USModel(pe), uptake_dentine=USModel(pd),
                      u234_u238_enamel=rin_e, u234_u238_dentine=rin_d, u234_u238_is="initial", **BASE)
    De = sum(c.accumulated(T) for c in fwd.components())
    return De, (th_e, r48e), (th_d, r48d)


@pytest.mark.parametrize("T,pe,pd,rin_e,rin_d", [
    (300.0, 0.4, -0.5, 1.25, 1.35),
    (80.0, -1.0, 0.0, 1.10, 1.10),
    (600.0, 2.0, 0.5, 1.0, 1.4),
])
def test_round_trip_recovers_age_and_p(T, pe, pd, rin_e, rin_d):
    De, (th_e, r_e), (th_d, r_d) = _synthetic(T, pe, pd, rin_e, rin_d)
    us = USESRSample(ToothSample(De=De, **BASE), enamel=UseriesData(th_e, r_e), dentine=UseriesData(th_d, r_d))
    res = us.age()
    assert res.status == "ok"
    assert res.age == pytest.approx(T, rel=1e-6)
    assert res.p_enamel == pytest.approx(pe, abs=1e-5)
    assert res.p_dentine == pytest.approx(pd, abs=1e-5)
    assert res.r_in_dentine == pytest.approx(rin_d, rel=1e-6)


def test_closed_system_age_matches_textbook_equation():
    # 230Th/234U age equation with the measured (present-day) 234U/238U
    l0, l4 = LAMBDA["Th230"], LAMBDA["U234"]
    T, r = 150.0, 1.3
    d = r - 1
    th238 = 1 - math.exp(-l0 * T) + d * l0 / (l0 - l4) * (1 - math.exp(-(l0 - l4) * T))
    assert closed_system_age(th238 / r, r) == pytest.approx(T, rel=1e-6)


def test_p_monotonic_and_bounds():
    th = th230_u234(200.0, 0.0, 1.2)
    assert solve_p(200.0, th, 1.2) == pytest.approx(0.0, abs=1e-6)
    # later uptake -> less 230Th
    assert th230_u234(200.0, 2.0, 1.2) < th < th230_u234(200.0, -1.0, 1.2)
    # an age younger than the closed-system age cannot reproduce the ratio
    assert solve_p(50.0, th, 1.2) is None
    assert incoming_ratio(200.0, -1.0, 1.0) == pytest.approx(1.0)


def test_no_solution_when_esr_younger_than_useries():
    _, enamel, dentine = _synthetic(300.0, 0.0, 0.0, 1.2, 1.2)
    us = USESRSample(ToothSample(De=5.0, **BASE), enamel=UseriesData(*enamel), dentine=UseriesData(*dentine))
    res = us.age()
    assert res.status == "no_solution" and res.age is None
    assert "leaching" in res.summary()


def test_tissue_without_data_keeps_its_uptake_model():
    De, enamel, _ = _synthetic(250.0, 0.3, 0.0, 1.2, 1.0)
    us = USESRSample(ToothSample(De=De, uptake_dentine=USModel(0.0), **BASE), enamel=UseriesData(*enamel))
    res = us.age()
    assert res.p_dentine is None
    assert res.age == pytest.approx(250.0, rel=1e-4)


def test_monte_carlo():
    De, (th_e, r_e), (th_d, r_d) = _synthetic(300.0, 0.4, -0.5, 1.25, 1.35)
    us = USESRSample(
        ToothSample(De=(De, 0.05 * De), **BASE),
        enamel=UseriesData((th_e, 0.01), (r_e, 0.005)),
        dentine=UseriesData((th_d, 0.01), (r_d, 0.005)),
    )
    mc = us.age_mc(n=60, seed=1)
    assert mc.ages.size + mc.n_failed == 60
    assert 250 < mc.mean < 350 and mc.std > 0
