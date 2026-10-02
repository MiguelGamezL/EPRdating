"""Regression tests against reference outputs of ROSY 2.0.

The reference data in ``rosy_reference/`` were produced with the original
ROSY program driven by ``tools/rosy_harness``. See ``ROSY_FINDINGS.md``.

What is checked here is the part of the physics that EPRdating already
implements in the same way as ROSY. Beta attenuation (one-group) is not
implemented yet, so beta factors measured from ROSY are used as inputs in the
end-to-end test.
"""

import json
from pathlib import Path

import pytest

from eprdating import BetaGeometry, Sediment, ToothLayers, ToothSample, USeries, USModel, conversion_factors
from eprdating.dose_rate import cosmic_dose_rate_sea_level, matrix_dose_rates, water_correction

pytestmark = pytest.mark.validation

REF = Path(__file__).parent / "rosy_reference"
CASES = json.loads((REF / "cases.json").read_text())["cases"]
RES = json.loads((REF / "results.json").read_text())
CF = conversion_factors("adamiec_aitken_1998")  # ROSY uses these factors (see findings)
CU = {r: CF.get("U", r).value for r in ("alpha", "beta", "gamma")}


def _v(x):
    """Central value of an input stored as a number or a [value, error] pair."""
    return x[0] if isinstance(x, list) else x


def _avg(rad, T, p=-1.0, ratio=1.0, radon_loss=0.0):
    us = USeries(ratio=ratio, ratio_is="initial", radon_loss=radon_loss)
    return USModel(p).accumulated(T, us.G(rad)) / T


ENV_CASES = [c for c in RES if CASES[c].get("env_option") == 2]


@pytest.mark.parametrize("cid", ENV_CASES)
def test_gamma_matches_rosy(cid):
    p, eu = CASES[cid], RES[cid]["EU"]
    m = matrix_dose_rates(_v(p["U_sed"]), _v(p["Th_sed"]), _v(p["K_sed"]), CF)
    ours = water_correction(m["gamma"], _v(p["water_sed"]) / 100, "gamma")
    assert eu["Total"][2] / 1000 == pytest.approx(ours, rel=0.005)


@pytest.mark.parametrize("cid", ENV_CASES)
def test_cosmic_matches_rosy(cid):
    p, eu = CASES[cid], RES[cid]["EU"]
    ours = cosmic_dose_rate_sea_level(_v(p["depth"]), _v(p["overburden_density"]))
    assert eu["Total"][3] / 1000 == pytest.approx(ours, rel=0.015)


ALPHA_CASES = [c for c in RES if c.startswith(("A_", "B_init", "D_"))]


@pytest.mark.parametrize("cid", ALPHA_CASES)
def test_enamel_alpha_ingrowth_consistent_with_rosy(cid):
    """ROSY's alpha efficiency varies with energy (k=0.15 at 5.3 MeV), so the
    effective k depends on which chain segments dominate. With our ingrowth
    model it stays within a few percent of 0.15 from 50 ka to 2 Ma."""
    p = CASES[cid]
    for mode, pp in (("EU", -1.0), ("LU", 0.0)):
        d = RES[cid][mode]
        T = d["age"] / 1000
        k_eff = (d["Enamel"][0] / 1000) / (p["U_en"] * CU["alpha"] * _avg("alpha", T, pp, p["ratio"]))
        assert k_eff == pytest.approx(0.15, rel=0.05)


E2E = [c for c in RES if c.startswith("H_")]


@pytest.mark.parametrize("cid", E2E)
def test_end_to_end_age_within_3_percent(cid):
    p = CASES[cid]
    for mode, pp in (("EU", -1.0), ("LU", 0.0)):
        s = ToothSample(
            De=p["De"], enamel_U=p["U_en"], dentine_U=p["U_den"],
            sediment=Sediment(U=p["U_sed"], Th=p["Th_sed"], K=p["K_sed"], water=p["water_sed"] / 100),
            beta=BetaGeometry(internal=0.69, dentine=0.14, external=0.148),  # measured from ROSY, 1000/2000 µm
            cosmic=cosmic_dose_rate_sea_level(p["depth"], p["overburden_density"]),
            k_alpha=0.15, uptake_enamel=USModel(pp), uptake_dentine=USModel(pp),
            factors="adamiec_aitken_1998",
        )
        assert s.age().age == pytest.approx(RES[cid][mode]["age"] / 1000, rel=0.03)


# --- full ages with one-group beta attenuation ------------------------------
AGE_CASES = [c for c in RES if c.startswith(("H_", "BR97"))]


@pytest.mark.parametrize("cid", AGE_CASES)
def test_ages_with_onegroup_within_2_percent_of_rosy(cid):
    """End-to-end ages, including the six teeth of Brennan et al. (1997),
    with beta attenuation computed by EPRdating's one-group solver."""
    p = CASES[cid]
    geo = ToothLayers(
        enamel_um=_v(p["thick_en"]), dentine_um=_v(p["thick_den"]),
        strip_outer_um=_v(p["strip_out"]), strip_inner_um=_v(p["strip_in"]),
        enamel_density=_v(p["density_en"]), dentine_density=_v(p["density_den"]),
        sediment_density=_v(p["density_sed"]),
    )
    if p["env_option"] == 2:
        env = {"gamma": None, "cosmic": cosmic_dose_rate_sea_level(_v(p["depth"]), _v(p["overburden_density"]))}
    else:
        env = {"gamma": _v(p["gamma_cosmic"]) / 1000, "cosmic": 0.0}
    cu = (-1.0 if p["uptake_en"] == 1 else 0.0, -1.0 if p["uptake_den"] == 1 else 0.0)
    for mode, (pe, pd) in (("EU", (-1.0, -1.0)), ("LU", (0.0, 0.0)), ("CU", cu)):
        s = ToothSample(
            De=_v(p["De"]), enamel_U=_v(p["U_en"]), dentine_U=_v(p["U_den"]),
            sediment=Sediment(U=_v(p["U_sed"]), Th=_v(p["Th_sed"]), K=_v(p["K_sed"]),
                              water=_v(p["water_sed"]) / 100),
            beta=geo, k_alpha=_v(p["alpha_eff"]), uptake_enamel=USModel(pe), uptake_dentine=USModel(pd),
            factors="adamiec_aitken_1998", **env,
        )
        assert s.age().age == pytest.approx(RES[cid][mode]["age"] / 1000, rel=0.02), mode
