"""Spectral processing: from cw-ESR spectra to intensities.

* :mod:`.io`, :mod:`.bruker`  reading spectra: Bruker BES3T and ESP/WinEPR,
                         ``.dat``/``.par`` ASCII pairs, plain columns.
* :mod:`.preprocess`     baseline, normalisation, alignment, modulation and
                         time-constant broadening of simulated shapes.
* :mod:`.intensity`      peak-to-peak, T1–B2 and double-integral intensities.
* :mod:`.deconvolution`  non-negative decomposition into component shapes.
* :mod:`.epraya_backend` (optional) component shapes simulated with EPRAYA.
"""

from .bruker import read_bes3t, read_esp
from .deconvolution import ComponentBasis, DeconvolutionResult, component_vs_dose
from .intensity import double_integral, field_for_g, g_for_field, peak_to_peak, t1_b2_amplitude
from .io import Spectrum, read_columns, read_dat, read_epr, read_par, read_series
from .preprocess import (
    aligned_average,
    noise_sigma,
    normalise,
    pseudo_modulation,
    subtract_baseline,
    time_constant_filter,
)

__all__ = [
    "ComponentBasis",
    "DeconvolutionResult",
    "Spectrum",
    "aligned_average",
    "component_vs_dose",
    "double_integral",
    "field_for_g",
    "g_for_field",
    "noise_sigma",
    "normalise",
    "peak_to_peak",
    "pseudo_modulation",
    "read_bes3t",
    "read_columns",
    "read_dat",
    "read_epr",
    "read_esp",
    "read_par",
    "read_series",
    "subtract_baseline",
    "t1_b2_amplitude",
    "time_constant_filter",
]
