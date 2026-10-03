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
from eprdating.series import LAMBDA, activity_ratio_Th230_U238

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
    s = USeries(ratio=1.0, partition=FAKE_PARTITION)
    g = s.G("alpha")
    # young: only U segments contribute (45 %); old: tends to equilibrium
    assert g(0.01) / 0.01 == pytest.approx(0.45, rel=1e-3)
    assert g(5000.0) / 5000.0 == pytest.approx(1.0, rel=0.05)


def test_ingrowth_makes_older_ages():
    eq = DoseRateComponent("int", 1.0, EarlyUptake())
    dis = DoseRateComponent("int", 1.0, EarlyUptake(), USeries(1.0, partition=FAKE_PARTITION).G("alpha"))
    assert solve_age(50.0, [dis]).age > solve_age(50.0, [eq]).age


def test_partition_validation():
    with pytest.raises(ValueError):
        USeries(partition={r: {"U238": 0.5, "U234": 0.2, "Th230": 0.2, "U235": 0.0} for r in ("alpha", "beta", "gamma")})
    with pytest.raises(ValueError, match="empty"):
        USeries(partition={r: {"U238": None, "Th230": 1.0} for r in ("alpha", "beta", "gamma")})


def test_bundled_partition_adamiec_aitken():
    s = USeries()
    for rad in ("alpha", "beta", "gamma"):
        assert sum(s.partition[rad].values()) == pytest.approx(1.0, abs=1e-6)
    # share of the equilibrium dose rate delivered by freshly incorporated U
    fresh = {rad: s.G(rad)(1e-3) / 1e-3 for rad in ("alpha", "beta", "gamma")}
    assert fresh["alpha"] == pytest.approx(0.203, abs=0.005)
    assert fresh["beta"] == pytest.approx(0.387, abs=0.005)
    assert fresh["gamma"] < 0.03
    # 234U excess raises the dose of young U
    assert USeries(ratio=1.5, ratio_is="initial").G("alpha")(10.0) > s.G("alpha")(10.0)


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
    eu = _tooth(uptake_enamel=EarlyUptake(), uptake_dentine=EarlyUptake()).age()
    lu = _tooth(uptake_enamel=LinearUptake(), uptake_dentine=LinearUptake()).age()
    assert lu.age > eu.age > 0
    assert sum(lu.accumulated.values()) == pytest.approx(200.0)


def test_equilibrium_opt_in_warns_and_is_younger():
    with pytest.warns(UserWarning, match="equilibrium"):
        eq = _tooth(ingrowth=False).age()
    assert _tooth().age().age > eq.age


def test_tooth_mc():
    mc = _tooth().age_mc(n=400, seed=1)
    assert mc.samples.size > 390
    assert mc.mean == pytest.approx(mc.nominal.age, rel=0.05)
    lo, hi = mc.interval()
    assert lo < mc.nominal.age < hi
    assert math.isfinite(mc.std)


# --- radon loss and 234U/238U ------------------------------------------------


def test_full_radon_loss_matches_pre_rn_factors():
    # Adamiec & Aitken (1998) Table 5, natural U: pre-Rn / full series
    expected = {"alpha": 1.26 / 2.78, "beta": 0.060 / 0.146, "gamma": 0.0044 / 0.113}
    s = USeries(radon_loss=1.0)
    for rad, ratio in expected.items():
        tau = 1e5  # long enough for full ingrowth
        assert s.G(rad)(tau) / tau == pytest.approx(ratio, rel=0.03), rad


def test_present_ratio_equals_back_corrected_initial():
    r_now, T = 1.3, 150.0
    present = USeries(ratio=r_now, ratio_is="present")
    r0 = present.initial_ratio(T)
    assert r0 == pytest.approx(1 + 0.3 * math.exp(LAMBDA["U234"] * T))
    initial = USeries(ratio=r0, ratio_is="initial")
    for rad in ("alpha", "beta", "gamma"):
        assert present.G(rad)(T) == pytest.approx(initial.G(rad)(T), rel=1e-12)


def test_useries_input_validation():
    with pytest.raises(ValueError):
        USeries(radon_loss=1.5)
    with pytest.raises(ValueError):
        USeries(ratio_is="today")
    with pytest.raises(ValueError, match="Rn222"):
        USeries(radon_loss=0.2, partition=FAKE_PARTITION)


def test_radon_loss_and_ratio_shift_tooth_age():
    base = _tooth().age().age
    assert _tooth(radon_loss_dentine=0.5, radon_loss_enamel=0.5).age().age > base
    assert _tooth(u234_u238_dentine=1.4).age().age < base
    # each tissue uses its own ratio
    assert _tooth(u234_u238_enamel=1.4).age().age != _tooth(u234_u238_dentine=1.4).age().age


def test_mc_with_uncertain_ratio_and_radon():
    mc = _tooth(u234_u238_dentine=(1.3, 0.05), radon_loss_dentine=(0.3, 0.2)).age_mc(n=200, seed=2)
    assert mc.samples.size == 200
    assert mc.interval()[0] < mc.nominal.age < mc.interval()[1]


def test_tooth_sample_with_onegroup_geometry():
    from eprdating import ToothLayers

    geo = ToothLayers(enamel_um=1000, strip_outer_um=50, strip_inner_um=50)
    s = _tooth(beta=geo)
    nominal = s.age()
    assert nominal.age > 0
    assert nominal.components["sediment beta"] > 0
    mc = s.age_mc(n=100, seed=3)
    assert mc.samples.size == 100


def test_geometry_uncertainty_widens_mc():
    from eprdating import ToothLayers

    geo = ToothLayers(enamel_um=(1000, 300), strip_outer_um=(50, 25), strip_inner_um=(50, 25))
    wide = _tooth(beta=geo, dentine_U=(30.0, 0.1)).age_mc(n=200, seed=5)
    narrow = _tooth(beta=geo, dentine_U=(30.0, 0.1), sample_geometry=False).age_mc(n=200, seed=5)
    assert wide.nominal.age == pytest.approx(narrow.nominal.age)
    assert wide.std > narrow.std


def test_alpha_efficiency_option():
    young_const = _tooth(De=(20.0, 1.0), enamel_U=(5.0, 0.1)).age().age
    young_energy = _tooth(De=(20.0, 1.0), enamel_U=(5.0, 0.1), alpha_efficiency="energy").age().age
    # young U is dominated by low-energy 238U/234U alphas -> lower k -> older age
    assert young_energy > young_const
    with pytest.raises(ValueError):
        _tooth(alpha_efficiency="whatever").age()


def test_beta_by_segment_option():
    from eprdating import LinearUptake, ToothLayers

    geo = ToothLayers(enamel_um=1000, strip_outer_um=50, strip_inner_um=50)
    # young tooth with U-rich dentine: before 226Ra grows in, the hard 234mPa
    # betas dominate, and they reach the enamel better than the chain average
    kw = dict(beta=geo, De=(10.0, 1.0), dentine_U=(50.0, 1.0), uptake_dentine=LinearUptake())  # noqa: C408
    seg = _tooth(**kw).age()
    chain = _tooth(beta_by_segment=False, **kw).age()
    assert seg.accumulated["dentine beta"] > 1.15 * chain.accumulated["dentine beta"]
    assert seg.age < chain.age
    # in equilibrium (old tooth, EU) the two are the same
    kw.update(De=(3000.0, 100.0), uptake_dentine=None)
    a, b = _tooth(**kw).age(), _tooth(beta_by_segment=False, **kw).age()
    assert a.age == pytest.approx(b.age, rel=0.01)
