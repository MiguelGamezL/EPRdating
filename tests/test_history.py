"""Time-varying sediment water and burial depth (piecewise histories)."""

import pytest

from eprdating import (
    BetaGeometry,
    History,
    Sediment,
    ToothLayers,
    ToothSample,
    USModel,
    Value,
    cosmic_dose_rate,
    cosmic_history,
    water_correction,
)
from eprdating.dose_rate import matrix_dose_rates
from eprdating.history import integrate
from eprdating.usesr import UseriesData, USESRSample, predicted_ratios

SED = dict(U=2.0, Th=6.0, K=1.0)  # noqa: C408
GEO = BetaGeometry(internal=1.0, dentine=0.0, external=0.3)


def _external_only(water, De=50.0, **kw):
    """No U in the tooth: only sediment beta, gamma and cosmic."""
    return ToothSample(De=De, enamel_U=0.0, dentine_U=0.0, sediment=Sediment(**SED, water=water),
                       beta=GEO, **{"cosmic": 0.15, **kw})


def _rate(w, cosmic=0.15):
    dry = matrix_dose_rates(**SED)
    return 0.3 * water_correction(dry["beta"], w, "beta") + water_correction(dry["gamma"], w, "gamma") + cosmic


def test_history_validation_and_helpers():
    h = History([0.1, (0.3, 0.05), Value(0.2, 0.02)], breaks=[10, 60])
    assert h.nominal() == [0.1, 0.3, 0.2] and len(h) == 3
    assert h.at(5) == 0.1 and h.at(10) == 0.3 and h.at(500) == 0.2
    assert h.mean(100) == pytest.approx((0.1 * 10 + 0.3 * 50 + 0.2 * 40) / 100)
    for bad in ([[0.1, 0.2], []], [[0.1, 0.2], [0]], [[0.1, 0.2, 0.3], [20, 10]]):
        with pytest.raises(ValueError):
            History(*bad)
    assert integrate([10, 60], [1.0, 2.0, 3.0], 5) == 5
    assert integrate([10, 60], [1.0, 2.0, 3.0], 70) == 10 + 100 + 30


def test_constant_history_reproduces_constant_inputs():
    ref = _external_only(0.15).age().age
    assert _external_only(History([0.15])).age().age == pytest.approx(ref, rel=1e-12)
    assert _external_only(History([0.15, 0.15], [20])).age().age == pytest.approx(ref, rel=1e-12)
    assert _external_only(0.15, cosmic=History([0.15, 0.15], [5])).age().age == pytest.approx(ref, rel=1e-12)


def test_two_segment_water_history_matches_hand_calculation():
    De, t1 = 50.0, 20.0
    r1, r2 = _rate(0.10), _rate(0.35)
    assert r1 * t1 < De  # the age falls in the second segment
    expected = t1 + (De - r1 * t1) / r2
    res = _external_only(History([0.10, 0.35], [t1]), De=De).age()
    assert res.age == pytest.approx(expected, rel=1e-9)
    # present-day rates are those of the first segment; the dose is the integral
    assert res.components["gamma"] + res.components["sediment beta"] + 0.15 == pytest.approx(r1)
    assert sum(res.accumulated.values()) == pytest.approx(De)
    # a wetter past lowers the dose rate: older than with today's water throughout
    assert res.age > _external_only(0.10, De=De).age().age


def test_cosmic_history_from_burial_depth():
    site = dict(density=1.8, lat_deg=4.6, lon_deg=-74.1, altitude_m=2600.0)  # noqa: C408
    h = cosmic_history(History([(3.0, 0.3), 0.5], breaks=[30]), **site)
    deep = cosmic_dose_rate(3.0, **site)
    assert h.values[0].value == pytest.approx(deep.value)
    assert h.values[0].sigma > deep.sigma  # the depth uncertainty is propagated
    assert h.values[1].value > h.values[0].value  # shallower in the past
    assert h.values[1].sigma == pytest.approx(cosmic_dose_rate(0.5, **site).sigma)
    # in the age: T = t1 + (De - R1 t1) / R2
    De, t1 = 50.0, 30.0
    c1, c2 = h.nominal()
    expected = t1 + (De - _rate(0.1, c1) * t1) / _rate(0.1, c2)
    assert _external_only(0.1, De=De, cosmic=h).age().age == pytest.approx(expected, rel=1e-9)


def test_in_situ_gamma_follows_the_water_history():
    w = History([0.10, 0.30], [15])
    comps = {c.name: c for c in _external_only(w, gamma=0.8).components()}
    rates = comps["gamma"].profile[1]
    assert rates[0] == 0.8
    assert rates[1] == pytest.approx(0.8 * (1 + 1.14 * 0.10) / (1 + 1.14 * 0.30))
    # a gamma history is used as given
    g = {c.name: c for c in _external_only(w, gamma=History([0.8, 0.6], [40])).components()}["gamma"]
    assert g.profile == ((40.0,), [0.8, 0.6])


def test_onegroup_geometry_and_monte_carlo():
    geo = ToothLayers(enamel_um=1000, strip_outer_um=50, strip_inner_um=50)
    kw = dict(De=(150.0, 10.0), enamel_U=1.0, dentine_U=20.0, beta=geo, cosmic=0.1,  # noqa: C408
              uptake_enamel=USModel(0.0), uptake_dentine=USModel(0.0))
    dry = ToothSample(sediment=Sediment(**SED, water=0.1), **kw)
    wet_past = ToothSample(sediment=Sediment(**SED, water=History([0.1, (0.4, 0.1)], [30])), **kw)
    a, b = dry.age(), wet_past.age()
    assert b.age > a.age
    # the internal components keep the present-day water
    assert b.components["dentine beta"] == pytest.approx(a.components["dentine beta"])
    mc_a, mc_b = dry.age_mc(n=300, seed=1), wet_past.age_mc(n=300, seed=1)
    assert mc_b.n_failed == 0 and mc_b.std > mc_a.std


def test_us_esr_with_a_water_history():
    base = dict(enamel_U=2.0, dentine_U=40.0, beta=ToothLayers(1000), cosmic=0.15)  # noqa: C408
    sed = Sediment(**SED, water=History([0.1, 0.3], [50]))
    T, pe, pd, r = 200.0, 0.2, -0.3, 1.2
    (r_e, th_e), (r_d, th_d) = predicted_ratios(T, pe, r), predicted_ratios(T, pd, r)
    fwd = ToothSample(De=1.0, sediment=sed, uptake_enamel=USModel(pe), uptake_dentine=USModel(pd),
                      u234_u238_enamel=r, u234_u238_dentine=r, u234_u238_is="initial", **base)
    De = sum(c.accumulated(T) for c in fwd.components())
    res = USESRSample(ToothSample(De=De, sediment=sed, **base), enamel=UseriesData(th_e, r_e),
                      dentine=UseriesData(th_d, r_d)).age()
    assert res.status == "ok" and res.age == pytest.approx(T, rel=1e-6)
