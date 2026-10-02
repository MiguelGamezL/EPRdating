"""EPRdating — ESR/EPR dating in Python.

Pipeline: spectra → intensity → dose-response → De → dose rate → age.
"""

from ._types import Value
from .age import AgeMC, AgeResult, DoseRateComponent, ToothSample, solve_age
from .beta import BetaGeometry
from .dose_rate import (
    Sediment,
    available_factor_sets,
    conversion_factors,
    cosmic_dose_rate,
    matrix_dose_rates,
    water_correction,
)
from .doseresponse import DoseResponseResult, bootstrap_De, fit_dose_response
from .onegroup import ToothLayers
from .series import USeries
from .uptake import EarlyUptake, LinearUptake, USModel

__version__ = "0.1.0.dev0"

__all__ = [
    "AgeMC",
    "AgeResult",
    "BetaGeometry",
    "DoseRateComponent",
    "DoseResponseResult",
    "EarlyUptake",
    "LinearUptake",
    "Sediment",
    "ToothLayers",
    "ToothSample",
    "USModel",
    "USeries",
    "Value",
    "available_factor_sets",
    "bootstrap_De",
    "conversion_factors",
    "cosmic_dose_rate",
    "fit_dose_response",
    "matrix_dose_rates",
    "solve_age",
    "water_correction",
]
