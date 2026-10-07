"""A real dose series: the Calio P4 enamel powder (Hakim et al. 2025, Nature).

Public MS5000 spectra (tests/data/calio_p4, from Zenodo 15515771, CC-BY 4.0):
natural and nine gamma doses, each aliquot in two or three rotations.
Published De (T1-B2 peak-to-peak, SSE): 2267 ± 99 Gy. See docs/validation.md.
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent / "validation"))
from calio_ms5000 import DATA, PUBLISHED, WINDOW, load

from eprdating import fit_dose_response
from eprdating.spectra import combined_intensity

pytestmark = pytest.mark.skipif(not DATA.is_dir(), reason="Calio spectra not available")


@pytest.fixture(scope="module")
def series():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        groups = load(DATA)
        doses = np.array(sorted(groups))
        I = np.array([combined_intensity(groups[d], "peak_to_peak", None, WINDOW, n_noise=50).value
                      for d in doses])
    return groups, doses, I


def test_files(series):
    groups, doses, _ = series
    assert list(doses) == [0, 50, 100, 250, 600, 1200, 2400, 4000, 8000, 15000]
    assert sum(len(v) for v in groups.values()) == 25
    s = groups[0][0]
    assert s.freq_GHz == pytest.approx(9.3895, abs=1e-3) and s.power_mW == 2.0 and s.mod_amp_mT == 0.1
    assert groups[100][0].n_scans == 10  # a rotation saved only as its ten runs


def test_published_De_with_the_8000_and_15000_Gy_labels_exchanged(series):
    _, doses, I = series
    # as labelled, the 8000 Gy aliquot lies far above the 15000 Gy one
    assert I[doses == 8000][0] > 1.2 * I[doses == 15000][0]
    swapped = doses.copy()
    swapped[doses == 8000], swapped[doses == 15000] = 15000, 8000
    f = fit_dose_response(swapped, I, "SSE", weighting="1/I^2")
    assert f.De == pytest.approx(PUBLISHED[0], abs=PUBLISHED[1])
    assert f.De_sigma < 1.5 * PUBLISHED[1]
    assert np.sqrt(f.chi2_red) < 0.04  # relative scatter about the curve
    labelled = fit_dose_response(doses, I, "SSE", weighting="1/I^2")
    assert np.sqrt(labelled.chi2_red) > 0.07 and labelled.De < 0.85 * PUBLISHED[0]
