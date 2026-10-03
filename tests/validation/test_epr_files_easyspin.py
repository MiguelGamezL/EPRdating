"""Bruker readers on EasySpin's collection of real files.

Not copied into this repository::

    git clone --depth 1 https://github.com/StollLab/EasySpin
    EPRDATING_EASYSPIN_FILES=EasySpin/tests/eprfiles pytest -m validation

The same measurements stored in two formats must give the same spectrum.
"""

import os
from pathlib import Path

import numpy as np
import pytest

from eprdating.spectra import read_epr

pytestmark = pytest.mark.validation
D = os.environ.get("EPRDATING_EASYSPIN_FILES")

PAIRS = [
    ("strong1.dsc", "strong1esp.par"),  # BES3T vs ESP (int32, other scale)
    ("esp.par", "win.par"),  # ESP vs WinEPR
    ("be3tintlit.dsc", "bes3tint.dsc"),  # little vs big endian
]


def _p(name):
    if not D:
        pytest.skip("set EPRDATING_EASYSPIN_FILES to EasySpin's tests/eprfiles")
    return Path(D) / name


@pytest.mark.parametrize("a,b", PAIRS)
def test_same_measurement_in_two_formats(a, b):
    x, y = read_epr(_p(a)), read_epr(_p(b))
    np.testing.assert_allclose(x.B, y.B)
    assert np.corrcoef(x.y, y.y)[0, 1] > 0.9999


@pytest.mark.parametrize("name,freq", [("E580_Xepr26b6_cwX_WillMyers.DSC", 9.487074), ("EMX_field1d.par", 9.702),
                                       ("99090211.dsc", 9.80963), ("frem_gly.par", 9.766096)])
def test_parameters(name, freq):
    s = read_epr(_p(name))
    assert s.freq_GHz == pytest.approx(freq)
    assert np.isfinite(s.y).all() and np.ptp(s.y) > 0
    assert 300 < s.B.mean() < 400  # X-band field sweeps, in mT
