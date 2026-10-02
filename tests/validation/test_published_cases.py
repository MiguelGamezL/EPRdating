"""Reproduce ages published with ROSY / DATA.

Each JSON file in ``cases/`` describes one published sample. A case runs
only when all its inputs and at least one expected age are filled in;
otherwise it is skipped with the list of missing fields, so this file doubles
as a to-do list for the validation campaign.

ROSY uses one-group beta attenuation; until that is implemented here, take
the beta geometry factors for each case from the original publication.
"""

import json
import warnings
from pathlib import Path

import pytest

from eprdating import BetaGeometry, Sediment, ToothSample, USeries, USModel

CASES = sorted((Path(__file__).parent / "cases").glob("*.json"))
UPTAKE = {"EU": -1.0, "LU": 0.0}


def _missing(d, prefix=""):
    out = []
    for k, v in d.items():
        if isinstance(v, dict):
            out += _missing(v, f"{prefix}{k}.")
        elif v is None and k != "u234_u238_initial":
            out.append(prefix + k)
    return out


def _v(x):
    return tuple(x) if isinstance(x, list) else x


@pytest.mark.validation
@pytest.mark.parametrize("path", CASES, ids=[p.stem for p in CASES])
def test_published_case(path):
    case = json.loads(path.read_text(encoding="utf-8"))
    inp = case["inputs"]
    missing = _missing(inp)
    expected = {k: v for k, v in case["expected_ages_ka"].items() if v is not None}
    if missing or not expected:
        pytest.skip(f"{case['id']}: missing {missing or 'expected ages'}")

    useries = USeries(r0=inp.get("u234_u238_initial") or 1.0)
    sed = inp["sediment"]
    for model, (age, _sigma) in expected.items():
        sample = ToothSample(
            De=_v(inp["De_Gy"]),
            enamel_U=_v(inp["enamel_U_ppm"]),
            dentine_U=_v(inp["dentine_U_ppm"]),
            sediment=Sediment(_v(sed["U_ppm"]), _v(sed["Th_ppm"]), _v(sed["K_pct"]), _v(sed["water"])),
            beta=BetaGeometry(**{k: _v(v) for k, v in inp["beta_geometry"].items()}),
            gamma=_v(inp["gamma_Gy_ka"]),
            cosmic=_v(inp["cosmic_Gy_ka"]),
            k_alpha=_v(inp["k_alpha"]),
            uptake_enamel=USModel(UPTAKE[model]),
            uptake_dentine=USModel(UPTAKE[model]),
            useries=useries,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            got = sample.age().age
        assert got == pytest.approx(age, rel=case["tolerance_rel"]), f"{model}: {got:.1f} vs {age} ka"
