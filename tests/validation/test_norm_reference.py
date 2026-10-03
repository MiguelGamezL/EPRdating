"""Gamma spectrometry of the corte 0 sediment (UNAL, HPGe 40 %, 24 h,
August 2024) against an independent analysis of the same spectra.

The spectra are not in the repository. Point ``EPRDATING_NORM_DATA`` to the
folder with ``U_24h.txt``, ``Th_24h.txt``, ``K_24h.txt``, ``corte0_24h.txt``,
``fondo_24h.txt`` and ``calibracion-57Co-22Na-137Cs-88Y_300s.txt``.

Reference (separate script, same spectra, 500 g for sample and references):
K = 3.32 ± 0.04 %, U = 1.74 ± 0.03 µg/g, Th = 7.92 ± 0.16 µg/g.
"""

import os
from pathlib import Path

import pytest

from eprdating.gamma import (
    CALIBRATION_LINES,
    IAEA_RGK_1,
    IAEA_RGTH_1,
    IAEA_RGU_1,
    Calibration,
    analyse,
    calibrate_natural,
    read_spectrum_txt,
)

pytestmark = pytest.mark.validation
DATA = os.environ.get("EPRDATING_NORM_DATA")


@pytest.fixture(scope="module")
def result():
    if not DATA:
        pytest.skip("set EPRDATING_NORM_DATA to the folder with the NORM spectra")
    d = Path(DATA)
    src = read_spectrum_txt(d / "calibracion-57Co-22Na-137Cs-88Y_300s.txt")
    lines = [e for k in ("57Co", "22Na", "137Cs", "88Y") for e in CALIBRATION_LINES[k]]
    cal = src.calibrate(lines, guess=Calibration((1.86, 0.3517), (0.5, 0.0005)))
    S = {k: read_spectrum_txt(d / f"{k}_24h.txt") for k in ("U", "Th", "K", "corte0", "fondo")}
    for s in S.values():
        calibrate_natural(s, cal, min_significance=5)
    refs = {"U": (S["U"], IAEA_RGU_1), "Th": (S["Th"], IAEA_RGTH_1), "K": (S["K"], IAEA_RGK_1)}
    return analyse(S["corte0"], 500.0, refs, background=S["fondo"])


def test_corte0_contents_match_independent_analysis(result):
    assert result.K.value == pytest.approx(3.32, rel=0.01)
    assert result.U.value == pytest.approx(1.74, rel=0.02)
    assert result.Th.value == pytest.approx(7.92, rel=0.01)


def test_corte0_lines_agree_within_group(result):
    assert result.groups["Ra226"].chi2_red < 2
    assert result.groups["Th232"].chi2_red < 3
