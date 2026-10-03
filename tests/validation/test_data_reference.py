"""Regression tests against reference outputs of DATA (Grün 2009).

The reference data in ``data_reference/`` were produced with the original
DATA.EXE run as a black box in DOSBox-X by ``tools/data_harness``: 41 cases
around a base tooth, each giving EU and LU dose rates and ages. See
``DATA_FINDINGS.md``.

The inputs are translated with DATA's conventions:

- water contents are given as % of the wet mass;
- the 234U/238U entered is the ratio of the incoming U (no back-correction);
- radon loss applies to the dentine only;
- one beta attenuation factor per tissue for the whole U chain
  (``beta_by_segment=False``); the default per-segment attenuation is checked
  separately;
- the cosmic dose rate is computed for an overburden of 2 g/cm³ at sea level
  (DATA has no site coordinates);
- dentine on either side is thick (5 mm), and sediment replaces it when its
  U content is zero.
"""

import json
from dataclasses import replace
from pathlib import Path

import pytest

from eprdating import Sediment, ToothLayers, ToothSample, USModel
from eprdating.dose_rate import cosmic_dose_rate_sea_level

pytestmark = pytest.mark.validation

REF = Path(__file__).parent / "data_reference"
CASES = json.loads((REF / "cases.json").read_text())
RES = json.loads((REF / "results.json").read_text(encoding="utf-8"))
UPTAKE = {"EU": -1.0, "LU": 0.0}


def _wet_to_dry(w_percent):
    return w_percent / (100.0 - w_percent)


def tooth(c, up, beta_by_segment=False):
    """A ToothSample with the inputs of a DATA case, in DATA's conventions."""
    def v(k):
        return c[k][0]

    side1, side2 = v("U_d1") > 0, v("U_d2") > 0
    layers = ToothLayers(
        enamel_um=v("thick"), dentine_um=5000.0 if side1 else 0.0, cementum_um=5000.0 if side2 else 0.0,
        strip_inner_um=v("rem1"), strip_outer_um=v("rem2"), enamel_density=v("dens"),
        dentine_density=2.82, cementum_density=2.82,
    )
    p = UPTAKE[up]
    return ToothSample(
        De=v("De"), enamel_U=v("U_en"), dentine_U=v("U_d1"), cementum_U=v("U_d2"),
        sediment=Sediment(U=v("U_s"), Th=v("Th_s"), K=v("K_s"), water=_wet_to_dry(v("W_s"))),
        beta=layers,
        cosmic=cosmic_dose_rate_sea_level(v("depth"), 2.0) if v("depth") > 0 else 0.0,
        gamma=v("ext_gamma") / 1000.0 if c["beta_only"] == "Y" else None,
        k_alpha=v("k"), dentine_water=_wet_to_dry(v("W_d")), cementum_water=_wet_to_dry(v("W_d")),
        uptake_enamel=USModel(p), uptake_dentine=USModel(p), uptake_cementum=USModel(p),
        u234_u238_enamel=v("r48"), u234_u238_dentine=v("r48"), u234_u238_cementum=v("r48"),
        u234_u238_is="initial",
        radon_loss_enamel=0.0, radon_loss_dentine=v("Rn") / 100, radon_loss_cementum=v("Rn") / 100,
        factors="adamiec_aitken_1998", beta_by_segment=beta_by_segment,
    )


def components(c, up, T, beta_by_segment=False):
    """Mean dose rates (µGy/a) over T ka, grouped as DATA prints them."""
    acc = {x.name: x.accumulated(T) / T * 1000 for x in tooth(c, up, beta_by_segment).components()}
    return {
        "gamma_sed": acc["gamma"] + acc["cosmic"],
        "internal": acc["enamel alpha"] + acc["enamel beta"],
        "beta_de1": acc["dentine beta"],
        "beta_de2": acc["cementum beta"],
        "beta_sed": acc["sediment beta"],
    }


RUNS = [(name, up) for name in RES for up in ("EU", "LU")]
_cache = {}


def _ours(name, up, beta_by_segment=False):
    key = (name, up, beta_by_segment)
    if key not in _cache:
        _cache[key] = components(CASES[name], up, RES[name]["age"][up][0], beta_by_segment)
    return _cache[key]


def _pairs(item, min_rate=20.0, beta_by_segment=False):
    """(DATA, EPRdating) dose rates of one component, where DATA prints ≥ min_rate."""
    out = []
    for name, up in RUNS:
        d = RES[name].get(item, {}).get(up, [0.0])[0]
        if d >= min_rate:
            out.append((name, up, d, _ours(name, up, beta_by_segment)[item]))
    return out


def _dev(pairs):
    return {f"{n}/{u}": 100 * (e / d - 1) for n, u, d, e in pairs}


def test_reference_is_complete():
    assert set(RES) == set(CASES) and len(RES) == 41
    assert all("age" in r for r in RES.values())


def test_gamma_and_cosmic():
    # same conversion factors, wet-mass water, Prescott & Hutton at 2 g/cm3
    for name, up, d, e in _pairs("gamma_sed"):
        assert e == pytest.approx(d, rel=0.005, abs=1.0), name


def test_internal_dose_rate():
    dev = _dev(_pairs("internal"))
    assert max(abs(x) for x in dev.values()) < 3.5, dev


def test_dentine_beta_single_factor():
    dev = _dev(_pairs("beta_de1") + _pairs("beta_de2"))
    assert -2.5 < min(dev.values()) and max(dev.values()) < 7.0, dev


def test_dentine_beta_per_segment_is_higher_for_young_teeth():
    # DATA attenuates the dose with ingrowth by one chain factor; per segment
    # the hard 234mPa betas, which dominate before 226Ra grows in, reach the
    # enamel better. The difference fades as the chain approaches equilibrium.
    seg = {(n, u): e / d for n, u, d, e in _pairs("beta_de1", beta_by_segment=True)}
    assert min(seg.values()) > 0.99
    assert seg[("De5", "LU")] > 1.30
    assert seg[("De3000", "EU")] == pytest.approx(1.0, abs=0.02)
    assert seg[("De5", "LU")] > seg[("B0", "LU")] > seg[("De1000", "LU")]


def test_sediment_beta_is_systematically_higher():
    # one-group attenuation from the sediment (as in ROSY, see ROSY_FINDINGS)
    # gives 7-13 % more than DATA's factors
    dev = _dev(_pairs("beta_sed"))
    assert len(dev) >= 20
    assert 5.0 < min(dev.values()) and max(dev.values()) < 15.0, dev


@pytest.mark.parametrize("name,up", RUNS)
def test_age_single_factor(name, up):
    d = RES[name]["age"][up][0]
    ours = replace(tooth(CASES[name], up), beta_by_segment=False).age().age
    # DATA prints ages with 2-3 significant digits
    assert ours == pytest.approx(d, rel=0.03, abs=0.06)


def test_age_default_bounds():
    devs = []
    for name, up in RUNS:
        d = RES[name]["age"][up][0]
        devs.append(100 * (tooth(CASES[name], up, beta_by_segment=True).age().age / d - 1))
    assert -8.5 < min(devs) and max(devs) < 2.0


def test_beta_attenuation_factors():
    # DATA prints its beta factors (dose in enamel / infinite-matrix dose)
    dent, sed = [], {"U": [], "Th": [], "K": []}
    for name, r in RES.items():
        geo = tooth(CASES[name], "EU").beta
        if "DE-1" in r.get("beta_corr_U", {}) and CASES[name]["U_d1"][0] > 0:
            dent.append(geo.chain_fraction("dentine", "U") / r["beta_corr_U"]["DE-1"] - 1)
        for nuc, devs in sed.items():
            if "SED" in r.get(f"beta_corr_{nuc}", {}):
                devs.append(geo.chain_fraction("sediment", nuc) / r[f"beta_corr_{nuc}"]["SED"] - 1)
    # the one-group dentine factors are DATA's (ROSY one-group corrected by
    # Marsh 1999) within the rounding of the printout, except at the extremes
    assert len(dent) > 30 and max(abs(x) for x in dent) < 0.04
    assert sum(abs(x) < 0.01 for x in dent) > 0.8 * len(dent)
    # sediment: U and Th 0-8 % higher, 40K 5-16 % higher
    assert all(0.0 < x < 0.08 for x in sed["U"] + sed["Th"])
    assert all(0.03 < x < 0.17 for x in sed["K"])
