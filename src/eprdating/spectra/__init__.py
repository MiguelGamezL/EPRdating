"""Spectral processing: from cw-ESR spectra to intensities.

* :mod:`.intensity`      peak-to-peak, T1–B2 and double-integral intensities.
* :mod:`.deconvolution`  non-negative decomposition into component shapes.
* :mod:`.epraya_backend` (optional) component shapes simulated with EPRAYA.
"""

from .deconvolution import ComponentBasis, DeconvolutionResult, component_vs_dose
from .intensity import double_integral, field_for_g, g_for_field, peak_to_peak, t1_b2_amplitude

__all__ = [
    "ComponentBasis",
    "DeconvolutionResult",
    "component_vs_dose",
    "double_integral",
    "field_for_g",
    "g_for_field",
    "peak_to_peak",
    "t1_b2_amplitude",
]
