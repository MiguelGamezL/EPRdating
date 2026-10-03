"""Gamma-ray spectrometry of sediments: U, Th and K by the comparative method.

* :mod:`.spectrum`     reading ASCII spectra, energy and resolution calibration.
* :mod:`.comparative`  peak areas, comparison with reference materials,
                       226Ra/238U equilibrium check, output as a
                       :class:`~eprdating.Sediment`.
"""

from .comparative import (
    BQ_PER_KG,
    IAEA_RGK_1,
    IAEA_RGTH_1,
    IAEA_RGU_1,
    LINES,
    GammaResult,
    Line,
    Reference,
    analyse,
    calibrate_natural,
    line_area,
    water_content,
)
from .spectrum import CALIBRATION_LINES, NATURAL_LINES, Calibration, GammaSpectrum, read_spectrum_txt

__all__ = [
    "BQ_PER_KG",
    "CALIBRATION_LINES",
    "IAEA_RGK_1",
    "IAEA_RGTH_1",
    "IAEA_RGU_1",
    "LINES",
    "NATURAL_LINES",
    "Calibration",
    "GammaResult",
    "GammaSpectrum",
    "Line",
    "Reference",
    "analyse",
    "calibrate_natural",
    "line_area",
    "read_spectrum_txt",
    "water_content",
]
