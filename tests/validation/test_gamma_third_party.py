"""Readers and automatic calibration on spectra from other laboratories.

Uses the sample spectra of the becquerel project (LBNL), which are not
copied into this repository::

    git clone --depth 1 https://github.com/lbl-anp/becquerel
    EPRDATING_BQ_SAMPLES=becquerel/tests/samples pytest -m validation

CNF and SPC files are read through becquerel itself (optional dependency).
For each HPGe spectrum the automatic calibration, done without a first
guess, must place the 351.9, 1460.8 and 2614.5 keV peaks within 0.5 keV and
agree with the calibration stored by the laboratory's software at 1460.8 keV.
(The stored calibrations are not always good at high energy: in Alcatraz14
the stored one puts 2614.5 keV 18 channels away from the actual peak.)
"""

import os
from pathlib import Path

import pytest

from eprdating.gamma import auto_calibrate, read_gamma
from eprdating.gamma.spectrum import fit_single_peak

pytestmark = pytest.mark.validation
D = os.environ.get("EPRDATING_BQ_SAMPLES")

HPGE = [
    "1110C NAA cave background May 2017.spe",  # ORTEC PopTop, lead cave
    "Mendocino_07-10-13_Acq-10-10-13.Spe",  # ORTEC PopTop, Marinelli beaker
    "01122014152731-GT01122014182338-GA37.4963000N-GO122.4633000W.cnf",  # Canberra Falcon 5000, field
    "Alcatraz14.Spc",  # ORTEC Trans-SPEC, field
]


def _path(name):
    if not D:
        pytest.skip("set EPRDATING_BQ_SAMPLES to becquerel's tests/samples folder")
    p = Path(D) / name
    if p.suffix.lower() in (".cnf", ".spc"):
        pytest.importorskip("becquerel")
    return p


@pytest.mark.parametrize("name", HPGE)
def test_auto_calibration_matches_stored_calibration(name):
    s = read_gamma(_path(name))
    stored = s.calibration
    cal = auto_calibrate(s)
    for e in (351.932, 1460.820, 2614.511):
        p = fit_single_peak(s, cal, e)
        assert p is not None and p["area"] > 5 * p["area_sigma"]
        assert cal.energy(p["centroid"]) == pytest.approx(e, abs=0.5)
    assert cal.energy(stored.channel(1460.82)) == pytest.approx(1460.82, abs=1.0)
    assert 2.3548 * cal.sigma_keV(1332) < 3.0  # HPGe resolution


@pytest.mark.parametrize("name", ["nai_detector.spe", "sim_spec.spe"])
def test_scintillator_spectra_are_rejected(name):
    s = read_gamma(_path(name))
    with pytest.raises(RuntimeError, match="HPGe"):
        auto_calibrate(s)
