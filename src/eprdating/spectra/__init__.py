"""Spectral processing: from cw-ESR spectra to intensities.

* :mod:`.io`, :mod:`.bruker`, :mod:`.freiberg`  reading spectra: Bruker
                         BES3T and ESP/WinEPR, Freiberg MS5000 (ESRStudio),
                         ``.dat``/``.par`` ASCII pairs, plain columns.
* :mod:`.preprocess`     baseline, normalisation, alignment, modulation and
                         time-constant broadening of simulated shapes.
* :mod:`.intensity`      intensity window; peak-to-peak, T1–B2 and
                         double-integral intensities.
* :mod:`.measure`        intensity of a spectrum inside the window, with errors.
* :mod:`.combine`        weighted average of repeated spectra of an aliquot.
* :mod:`.fragments`      (optional) enamel fragments measured at several
                         angles: merged spectrum, angular profile, X-ray
                         calibration.
* :mod:`.deconvolution`  non-negative decomposition into component shapes.
* :mod:`.epraya_backend` (optional) component shapes simulated with EPRAYA.
"""

from .bruker import read_bes3t, read_esp
from .combine import Combination, combine_spectra, combined_intensity
from .deconvolution import ComponentBasis, DeconvolutionResult, component_vs_dose
from .fragments import (
    AngularSpectrum,
    XrayCalibration,
    angular_profile,
    angular_set,
    fragment_intensity,
    merge_angular,
    parse_angular_name,
)
from .freiberg import read_ms5000
from .intensity import (
    DEFAULT_WINDOW,
    IntensityWindow,
    double_integral,
    field_for_g,
    g_for_field,
    peak_to_peak,
    t1_b2_amplitude,
)
from .io import Spectrum, read_columns, read_dat, read_epr, read_par, read_series
from .measure import Intensity, combine_intensities, empirical_template, intensity
from .preprocess import (
    aligned_average,
    noise_sigma,
    normalise,
    pseudo_modulation,
    subtract_baseline,
    time_constant_filter,
)

__all__ = [
    "DEFAULT_WINDOW",
    "AngularSpectrum",
    "Combination",
    "ComponentBasis",
    "DeconvolutionResult",
    "Intensity",
    "IntensityWindow",
    "Spectrum",
    "XrayCalibration",
    "aligned_average",
    "angular_profile",
    "angular_set",
    "combine_intensities",
    "combine_spectra",
    "combined_intensity",
    "component_vs_dose",
    "double_integral",
    "empirical_template",
    "field_for_g",
    "fragment_intensity",
    "g_for_field",
    "intensity",
    "merge_angular",
    "noise_sigma",
    "normalise",
    "parse_angular_name",
    "peak_to_peak",
    "pseudo_modulation",
    "read_bes3t",
    "read_columns",
    "read_dat",
    "read_epr",
    "read_esp",
    "read_ms5000",
    "read_par",
    "read_series",
    "subtract_baseline",
    "t1_b2_amplitude",
    "time_constant_filter",
]
