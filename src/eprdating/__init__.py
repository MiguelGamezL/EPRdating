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
    cosmic_history,
    matrix_dose_rates,
    water_correction,
)
from .doseresponse import DoseResponseResult, bootstrap_De, fit_dose_response
from .history import History
from .onegroup import Material, ToothLayers, compound, dentine_material, mixture, sediment_material
from .series import USeries
from .uptake import DelayedUptake, EarlyUptake, LinearUptake, USModel
from .usesr import UseriesData, USESRSample

__version__ = "0.2.0"

__all__ = [
    "AgeMC",
    "AgeResult",
    "BetaGeometry",
    "DelayedUptake",
    "DoseRateComponent",
    "DoseResponseResult",
    "EarlyUptake",
    "History",
    "LinearUptake",
    "Material",
    "Sediment",
    "ToothLayers",
    "ToothSample",
    "USESRSample",
    "USModel",
    "USeries",
    "UseriesData",
    "Value",
    "available_factor_sets",
    "bootstrap_De",
    "compound",
    "conversion_factors",
    "cosmic_dose_rate",
    "cosmic_history",
    "dentine_material",
    "fit_dose_response",
    "matrix_dose_rates",
    "mixture",
    "sediment_material",
    "solve_age",
    "water_correction",
]
