"""Recompute published ESR tooth ages (benchmark in the spirit of DRAC's validation).

Inputs and published ages are in ``published/*.json`` (format in
``published/SCHEMA.md``), transcribed from the papers and checked
independently against the sources. ``published_cases.py`` translates them
with the conventions of the program each study used. See
``PUBLISHED_FINDINGS.md``.
"""

import math

import pytest
from published_cases import benchmark, build, load_studies, recompute, unconverged_dentine_p

pytestmark = pytest.mark.validation

STUDIES = load_studies()
ROWS = benchmark()
COMPARED = [r for r in ROWS if not r.get("excluded") and r.get("ours")]


def _dev(r):
    return 100 * (r["ours"] / r["published"] - 1)


def _published_sigma(r):
    s = next(x for x in STUDIES[r["study"]]["samples"] if x["sample"] == r["sample"])
    a = s["published"]["ages_ka"][r["model"]]
    up = a[1]
    down = a[2] if len(a) > 2 and a[2] is not None else a[1]
    f = 0.5 if r["study"].startswith("dirks") else 1.0  # Dirks et al. print 2σ
    return f * (up if r["ours"] > r["published"] else down)


def test_benchmark_size():
    assert len({r["study"] for r in COMPARED}) == 9
    assert len(COMPARED) == 58


def test_every_age_within_ten_percent():
    worst = max(COMPARED, key=lambda r: abs(_dev(r)))
    assert abs(_dev(worst)) < 10, worst


def test_most_ages_within_five_percent_and_one_sigma():
    assert sum(abs(_dev(r)) <= 5 for r in COMPARED) >= 0.75 * len(COMPARED)
    inside = sum(abs(r["ours"] - r["published"]) <= _published_sigma(r) for r in COMPARED)
    assert inside >= len(COMPARED) - 2


def test_systematic_offset():
    # EPRdating's beta doses are a few % higher than DATA's and USESR's
    mean = sum(_dev(r) for r in COMPARED) / len(COMPARED)
    assert -5 < mean < -1


@pytest.mark.parametrize("study,lo,hi", [
    ("richard2023_lovedale", -4, 0),          # DATA, EU
    ("duval2019_khok_sung", -6, 1),           # DATA, US (converged) and CS-US
    ("bahain2020_tourville_cenieh", -4.5, -2.5),
    ("bahain2020_tourville_mnhn", -7, 3.5),   # USESR, Rn loss per tissue
    ("dirks2017_rising_star", -8, -1),        # USESR, 234U/238U ~6
    ("rizal2020_ngandong", -7, -2.5),
    ("duval2021_olieboomspoort", -1.5, 0),
    ("falgueres2025_fumane", -10, 0.5),
    ("carvajal2011_aguazuque", 2, 5),         # ROSY, external dose only
])
def test_study_ranges(study, lo, hi):
    d = [_dev(r) for r in COMPARED if r["study"] == study]
    assert d and lo < min(d) and max(d) < hi, d


def test_khok_sung_dentine_p_not_converged():
    # the DATA failure found with the harness appears in a published study:
    # every dentine much later than its enamel has a p that does not give back
    # the measured 230Th/234U; the enamel p always do
    st = STUDIES["duval2019_khok_sung"]
    flagged = [s["sample"] for s in st["samples"] if unconverged_dentine_p(s)]
    assert len(flagged) == 9
    assert not any(x.startswith("3548") for x in flagged)  # dentine not late: converged
    late = [r for r in ROWS if r["study"] == "duval2019_khok_sung" and "dentine p" in (r.get("excluded") or "")]
    # where the dentine carries U (3545, 3547, 3549) the published ages are 5-13 % too young
    assert max(_dev(r) for r in late) > 10


def test_lovedale_dose_components():
    # DATA's U-series screen: 234U/238U back-corrected with the closed-system age
    st = STUDIES["richard2023_lovedale"]
    for s in st["samples"]:
        res = recompute(build(st, s))["EU"]
        pub = s["published"]["dose_rates_uGy_a"]
        dr = res["dose_rates"]
        assert dr["internal"] == pytest.approx(pub["internal"][0], rel=0.04)
        assert dr["beta_dentine"] == pytest.approx(pub["beta_dentine"][0], rel=0.06)
        assert dr["gamma_cosmic"] == pytest.approx(pub["gamma_cosmic"][0], abs=1)


def test_carvajal_external_dose():
    st = STUDIES["carvajal2011_aguazuque"]
    rc = build(st, st["samples"][0])
    rates = {c.name: c.rate * 1000 for c in rc.tooth.components()}
    assert rates["gamma"] == pytest.approx(347.31, rel=0.03)
    assert rates["cosmic"] == pytest.approx(251.0)


def test_excluded_cases_are_documented():
    for sid, st in STUDIES.items():
        if not st.get("benchmark"):
            assert st.get("benchmark_reason"), sid
    assert any(r.get("excluded", "").startswith("the printed dentine U") for r in ROWS)
    assert all(math.isfinite(r["ours"]) for r in COMPARED)
