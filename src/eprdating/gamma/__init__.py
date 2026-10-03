"""Gamma-ray spectrometry of sediments: U, Th and K by the comparative method.

* :mod:`.readers`      reading ORTEC .Spe/.Chn, N42, ASCII and column files
                       (CNF, SPC, IEC through the optional package becquerel).
* :mod:`.spectrum`     energy/resolution calibration, automatic calibration.
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
    analyse_files,
    calibrate_natural,
    line_area,
    water_content,
)
from .readers import read_chn, read_columns, read_gamma, read_n42, read_spe
from .spectrum import (
    CALIBRATION_LINES,
    NATURAL_LINES,
    Calibration,
    GammaSpectrum,
    auto_calibrate,
    find_peaks_channels,
    read_spectrum_txt,
)

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
    "analyse_files",
    "auto_calibrate",
    "calibrate_natural",
    "find_peaks_channels",
    "line_area",
    "read_chn",
    "read_columns",
    "read_gamma",
    "read_n42",
    "read_spe",
    "read_spectrum_txt",
    "water_content",
]
