"""Validation of eprdating.usesr against published US-ESR results.

See USESR_FINDINGS.md. The U-series part (p–T relation) is checked tightly;
dose rates and ages loosely, because the published calculations use beta
attenuation conventions (Marsh 1999; USESR program) that differ from the
one-group (ROSY) attenuation used here, and some geometry is not reported.
"""

import json
from pathlib import Path

import pytest

from eprdating import Sediment, ToothLayers, ToothSample, USModel
from eprdating.usesr import UseriesData, USESRSample, incoming_ratio, solve_p

pytestmark = pytest.mark.validation

REF = json.loads((Path(__file__).parent / "usesr_reference" / "published_usesr.json").read_text())
SHAO = REF["shao2015"]
DN = REF["denadale2026"]


def _t(x):
    return tuple(x[:2])


# --- U-series part: p at the published age -----------------------------------
@pytest.mark.parametrize("tissue", ["enamel", "dentine", "cementum"])
def test_shao2015_p_at_published_age(tissue):
    d = SHAO["inputs"][tissue]
    p = solve_p(SHAO["results"]["age_ka"][0], d["th230_u234"][0], d["u234_u238"][0])
    assert p == pytest.approx(SHAO["results"][f"p_{tissue}"][0], abs=0.015)


@pytest.mark.parametrize("sample", list(DN["samples"]))
@pytest.mark.parametrize("tissue", ["enamel", "dentine"])
def test_denadale_p_at_published_age(sample, tissue):
    s = DN["samples"][sample]
    d = s[tissue]
    p = solve_p(s["results"]["age_ka"][0], d["th230_u234"][0], d["u234_u238"][0])
    # published ages are rounded to 1 ka, which moves p by up to ~0.03
    assert p == pytest.approx(s["results"][f"p_{tissue}"][0], abs=0.03)


# --- dose rates and ages ------------------------------------------------------
def _denadale_tooth(s, De=None, **kw):
    c = DN["common"]["sediment"]
    w = c["water_pct_wet"][0] / (100 - c["water_pct_wet"][0])
    e = s["enamel"]
    return ToothSample(
        De=_t(s["De_Gy"]) if De is None else De, enamel_U=_t(e["U_ppm"]), dentine_U=_t(s["dentine"]["U_ppm"]),
        sediment=Sediment(U=_t(c["U_ppm"]), Th=_t(c["Th_ppm"]), K=_t(c["K_pct"]), water=(w, 0.25 * w)),
        beta=ToothLayers(enamel_um=_t(e["thickness_um"]), strip_outer_um=_t(e["removed_side1_um"]),
                         strip_inner_um=_t(e["removed_side2_um"])),
        gamma=(0.333, 0.014), cosmic=(0.084, 0.030), k_alpha=(0.13, 0.02),
        enamel_water=(0.03, 0.03), dentine_water=(0.05, 0.05), factors="adamiec_aitken_1998", **kw,
    )


@pytest.mark.parametrize("sample", ["CN-857", "CN-55"])
def test_denadale_internal_dose_at_published_age_and_p(sample):
    s = DN["samples"][sample]
    T, pe, pd = (s["results"][k][0] for k in ("age_ka", "p_enamel", "p_dentine"))
    tooth = _denadale_tooth(
        s, De=1.0, uptake_enamel=USModel(pe), uptake_dentine=USModel(pd), u234_u238_is="initial",
        u234_u238_enamel=incoming_ratio(T, pe, s["enamel"]["u234_u238"][0]),
        u234_u238_dentine=incoming_ratio(T, pd, s["dentine"]["u234_u238"][0]),
    )
    acc = {c.name: c.accumulated(T) / T * 1000 for c in tooth.components()}
    internal = acc["enamel alpha"] + acc["enamel beta"]
    assert internal == pytest.approx(s["results"]["internal_uGy_a"][0], rel=0.05)


@pytest.mark.parametrize("sample", ["CN-857", "CN-1322"])
def test_denadale_age_within_published_uncertainty(sample):
    s = DN["samples"][sample]
    us = USESRSample(_denadale_tooth(s),
                     enamel=UseriesData(_t(s["enamel"]["th230_u234"]), _t(s["enamel"]["u234_u238"])),
                     dentine=UseriesData(_t(s["dentine"]["th230_u234"]), _t(s["dentine"]["u234_u238"])))
    r = us.age()
    age, err = s["results"]["age_ka"]
    assert r.status == "ok"
    assert abs(r.age - age) < err


def test_denadale_cn55_has_no_solution():
    """CN-55 sits at the closed-system U-series bound (published p ≈ -0.9 in
    both tissues; the authors flag it as unreliable). With one-group beta
    doses, which are higher than USESR's, the dose is reached before that
    bound, so EPRdating reports no solution instead of an age."""
    s = DN["samples"]["CN-55"]
    us = USESRSample(_denadale_tooth(s),
                     enamel=UseriesData(_t(s["enamel"]["th230_u234"]), _t(s["enamel"]["u234_u238"])),
                     dentine=UseriesData(_t(s["dentine"]["th230_u234"]), _t(s["dentine"]["u234_u238"])))
    assert us.age().status == "no_solution"


def test_shao2015_age_within_published_uncertainty():
    i = SHAO["inputs"]
    w = i["sediment"]["water_pct"][0] / (100 - i["sediment"]["water_pct"][0])
    tooth = ToothSample(
        De=_t(i["De_Gy"]), enamel_U=_t(i["enamel"]["U_ppm"]), dentine_U=_t(i["dentine"]["U_ppm"]),
        cementum_U=_t(i["cementum"]["U_ppm"]),
        radon_loss_enamel=1 - i["enamel"]["pb210_th230"][0], radon_loss_dentine=1 - i["dentine"]["pb210_th230"][0],
        radon_loss_cementum=1 - i["cementum"]["pb210_th230"][0],
        dentine_water=0.07, cementum_water=0.07,
        sediment=Sediment(U=_t(i["sediment"]["U_ppm"]), Th=_t(i["sediment"]["Th_ppm"]), K=_t(i["sediment"]["K_pct"]),
                          water=w),
        beta=ToothLayers(enamel_um=_t(i["enamel"]["thickness_um"]), strip_outer_um=_t(i["enamel"]["removed_side1_um"]),
                         strip_inner_um=_t(i["enamel"]["removed_side2_um"]), cementum_um=1000.0),  # assumed
        gamma=0.832, cosmic=0.0, k_alpha=0.13, factors="guerin_2011",
    )
    us = USESRSample(tooth,
                     enamel=UseriesData(i["enamel"]["th230_u234"][0], i["enamel"]["u234_u238"][0]),
                     dentine=UseriesData(i["dentine"]["th230_u234"][0], i["dentine"]["u234_u238"][0]),
                     cementum=UseriesData(i["cementum"]["th230_u234"][0], i["cementum"]["u234_u238"][0]))
    r = us.age()
    age, up, down = SHAO["results"]["age_ka"]
    assert age + down < r.age < age + up
