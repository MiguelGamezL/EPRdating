"""US-ESR and CSUS-ESR against the U-series/ESR part of DATA (Grün 2009).

The reference data in ``data_reference/useries_*`` were produced with the
original DATA.EXE (<F6> "read enamel into U-series") driven by
``tools/data_harness/run_useries.py``: 40 cases around a base tooth with
234U/238U and 230Th/234U for the enamel and each dentine side. See
``DATA_FINDINGS.md``.

Inputs follow DATA's conventions as in ``test_data_reference.py``; in
addition DATA ignores radon loss in its U-series/ESR calculations.
"""

import copy
import json
from pathlib import Path

import pytest
from test_data_reference import tooth  # tests/validation is on sys.path (rootdir import mode)

from eprdating.usesr import UseriesData, USESRSample, th230_u234

pytestmark = pytest.mark.validation

REF = Path(__file__).parent / "data_reference"
CASES = json.loads((REF / "useries_cases.json").read_text())
RES = json.loads((REF / "useries_results.json").read_text(encoding="utf-8"))

# DATA's dentine p does not reproduce the measured dentine ratio when the
# dentine took up its U much later than the enamel (see the test below)
DATA_UNCONVERGED = {"US_de_late", "US_d1d2", "US_dl05", "US_dl15"}
# DATA reports "U-SERIES TOO HIGH: NO RESULT" for solutions with p < ~-0.85
DATA_BOUND = {"US_cs_bound", "US_th05", "US_pd0.90", "US_pd0.95"}
NO_SOLUTION = {"US_De30", "US_nosol"}


def usesr(name, beta_by_segment=False):
    c = copy.deepcopy(CASES[name])
    c["Rn"] = [0.0, 0.0]  # not used by DATA's U-series/ESR part
    u = c["useries"]

    def data(t):
        (r48, _), (th, _) = u[t]
        return UseriesData(th, r48)

    return USESRSample(
        tooth(c, "EU", beta_by_segment),
        enamel=data("enamel"),
        dentine=data("de1") if c["U_d1"][0] > 0 else None,
        cementum=data("de2") if c["U_d2"][0] > 0 else None,
    )


_cache = {}


def _age(name, model="US"):
    if (name, model) not in _cache:
        _cache[(name, model)] = usesr(name).age(model=model)
    return _cache[(name, model)]


SOLVED = [n for n, r in RES.items() if "age" in r.get("usesr", {}) and n not in DATA_UNCONVERGED]
CSUS = [n for n, r in RES.items() if "csus" in r and n not in NO_SOLUTION]


def test_reference_is_complete():
    assert set(RES) == set(CASES) and len(RES) == 40
    assert len(SOLVED) == 29
    assert "crash" in RES["US_SS"]  # DATA's US-ESR fails without U in the dentine


@pytest.mark.parametrize("name", SOLVED)
def test_us_esr_age_and_p(name):
    d, ours = RES[name]["usesr"], _age(name)
    assert ours.status == "ok"
    # DATA prints the age with 2-4 digits
    assert ours.age == pytest.approx(d["age"][0], rel=0.015, abs=0.5)
    assert ours.p_enamel == pytest.approx(d["p_internal"][0], rel=0.02, abs=0.05)
    if ours.p_dentine is not None:
        assert ours.p_dentine == pytest.approx(d["p_beta_de1"][0], rel=0.02, abs=0.05)
    if ours.p_cementum is not None:
        assert ours.p_cementum == pytest.approx(d["p_beta_de2"][0], rel=0.02, abs=0.05)


@pytest.mark.parametrize("name", SOLVED)
def test_data_p_reproduces_the_measured_ratios(name):
    # the same U-series model: DATA's own (T, p) give back the measured ratios
    d, u = RES[name]["usesr"], CASES[name]["useries"]
    T = d["age"][0]
    for tissue, key in (("enamel", "p_internal"), ("de1", "p_beta_de1"), ("de2", "p_beta_de2")):
        if key in d:
            (r48, _), (th, _) = u[tissue]
            assert th230_u234(T, d[key][0], r48) == pytest.approx(th, abs=0.004)


@pytest.mark.parametrize("name", sorted(DATA_UNCONVERGED))
def test_data_dentine_p_unconverged(name):
    d, u = RES[name]["usesr"], CASES[name]["useries"]
    (r48, _), (th, _) = u["de1"]
    # DATA's dentine p predicts a 230Th/234U far from the one entered ...
    assert th230_u234(d["age"][0], d["p_beta_de1"][0], r48) - th > 0.02
    # ... EPRdating's solution reproduces it and gives a dentine p that is higher
    ours = _age(name)
    assert th230_u234(ours.age, ours.p_dentine, r48) == pytest.approx(th, abs=1e-6)
    assert ours.p_dentine > d["p_beta_de1"][0] + 0.5
    assert ours.age >= d["age"][0] - 0.5


@pytest.mark.parametrize("name", sorted(DATA_BOUND))
def test_near_the_bound_data_gives_no_result(name):
    assert "NO RESULT" in RES[name]["usesr"]["no_result"]
    ours = _age(name)
    assert ours.status == "ok"
    p = [x for x in (ours.p_enamel, ours.p_dentine, ours.p_cementum) if x is not None]
    assert min(p) < -0.85


@pytest.mark.parametrize("name", sorted(NO_SOLUTION))
def test_no_solution(name):
    assert "NO RESULT" in RES[name]["usesr"]["no_result"]
    assert _age(name).status == "no_solution"


def test_csus_younger_than_uptake_is_not_reported():
    # DATA prints a CS-US age (20 ka) younger than the U uptake it assumes
    # (31 and 38 ka); EPRdating reports that there is no solution
    assert RES["US_De30"]["csus"]["age"][0] == pytest.approx(20, abs=1)
    r = _age("US_De30", "CSUS")
    assert r.status == "no_solution" and r.min_age > 30


@pytest.mark.parametrize("name", CSUS)
def test_csus_age_and_u_dose(name):
    d, ours = RES[name]["csus"], _age(name, "CSUS")
    assert ours.status == "ok"
    # DATA prints CS-US ages rounded to 1 ka
    assert abs(ours.age - d["age"][0]) <= 0.5 + 0.015 * d["age"][0]
    # DATA's second CS-US column: dose from the tooth's own U since uptake (Gy)
    acc = ours.detail.accumulated
    u_dose = sum(acc[k] for k in ("enamel alpha", "enamel beta", "dentine beta", "cementum beta"))
    assert 0.94 < d["col2"][0] / u_dose < 1.0


def test_radon_loss_is_ignored_by_data():
    assert RES["US_Rn50"]["usesr"] == RES["US_B0"]["usesr"]
    assert RES["US_Rn50"]["csus"] == RES["US_B0"]["csus"]


@pytest.mark.parametrize("name", ["US_B0", "US_r20", "US_th80", "US_th95", "US_rmix"])
def test_eu_in_the_useries_screen_back_corrects_with_the_closed_system_age(name):
    # in the U-series screen DATA's EU/LU columns take the measured 234U/238U
    # back to the initial ratio over each tissue's closed-system U-series age
    import math
    from dataclasses import replace

    from eprdating.series import LAMBDA
    from eprdating.usesr import closed_system_age

    c = copy.deepcopy(CASES[name])
    c["Rn"] = [0.0, 0.0]
    u = c["useries"]

    def r0(t):
        (r, _), (th, _) = u[t]
        return 1 + (r - 1) * math.exp(LAMBDA["U234"] * closed_system_age(th, r))

    t = replace(tooth(c, "EU"), u234_u238_enamel=r0("enamel"), u234_u238_dentine=r0("de1"),
                u234_u238_cementum=r0("de2"), u234_u238_is="initial")
    T = RES[name]["age"]["EU"][0]
    acc = {x.name: x.accumulated(T) / T * 1000 for x in t.components()}
    assert acc["enamel alpha"] + acc["enamel beta"] == pytest.approx(RES[name]["internal"]["EU"][0], rel=0.02)
