"""Age calculation: solve ∫_0^T Ḋ(t) dt = De for the age T.

Generic layer
-------------
:class:`DoseRateComponent` is one contribution to the dose rate with its
present-day (equilibrium) value, an optional uptake model and an optional
U-series ingrowth function. :func:`solve_age` combines any list of them.

Tooth-enamel layer
------------------
:class:`ToothSample` assembles the usual components of ESR dating of
enamel (internal alpha and beta, beta from dentine and sediment, gamma,
cosmic) and provides a Monte Carlo age with :meth:`ToothSample.age_mc`.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from ._types import ValueLike, as_value
from .beta import BetaGeometry
from .dose_rate import (
    DEFAULT_FACTORS,
    K_ENAMEL,
    ConversionFactors,
    Sediment,
    conversion_factors,
    matrix_dose_rates,
    water_correction,
)
from .series import USeries
from .uptake import USModel


@dataclass
class DoseRateComponent:
    """One contribution to the total dose rate (Gy/ka).

    rate   : present-day dose rate the source would deliver in secular
             equilibrium with its present U content.
    uptake : None for a constant source, else a :class:`USModel`.
    G      : time-integrated activity ratio of incorporated U
             (``tau -> ∫ D/D_eq``); None means secular equilibrium.
    """

    name: str
    rate: float
    uptake: USModel | None = None
    G: Callable[[float], float] | None = None

    def accumulated(self, T: float) -> float:
        """Dose (Gy) delivered by this component over ``T`` ka."""
        if self.rate == 0:
            return 0.0
        if self.uptake is None:
            if self.G is None:
                return self.rate * T
            return self.rate * self.G(T)
        return self.rate * self.uptake.accumulated(T, self.G)


@dataclass
class AgeResult:
    age: float  # ka
    De: float  # Gy
    components: dict[str, float]  # present-day dose rates, Gy/ka
    accumulated: dict[str, float]  # dose delivered by each component, Gy

    @property
    def mean_dose_rate(self) -> float:
        """Time-averaged dose rate De / T (Gy/ka)."""
        return self.De / self.age

    def summary(self) -> str:
        lines = [f"Age = {self.age:.4g} ka   (De = {self.De:.4g} Gy, <Ḋ> = {self.mean_dose_rate:.4g} Gy/ka)"]
        for k, r in self.components.items():
            share = 100 * self.accumulated[k] / self.De
            lines.append(f"  {k:16s} Ḋ_now = {r:8.4f} Gy/ka   contributes {share:5.1f} % of De")
        return "\n".join(lines)


def solve_age(De: float, components: Sequence[DoseRateComponent], t_max: float = 1e5) -> AgeResult:
    """Find the age T (ka) such that the accumulated dose equals ``De`` (Gy)."""
    if De <= 0:
        raise ValueError("De must be positive")
    comps = list(components)

    def f(T: float) -> float:
        return sum(c.accumulated(T) for c in comps) - De

    lo, hi = 1e-9, 1.0
    while f(hi) < 0:
        hi *= 4.0
        if hi > t_max:
            raise RuntimeError(f"no solution below {t_max} ka: dose rate too low for this De")
    T = brentq(f, lo, hi, xtol=1e-10, rtol=1e-10)
    return AgeResult(
        age=T,
        De=De,
        components={c.name: c.rate for c in comps},
        accumulated={c.name: c.accumulated(T) for c in comps},
    )


# --------------------------------------------------------------------------
# Tooth enamel
# --------------------------------------------------------------------------


@dataclass
class AgeMC:
    """Monte Carlo age distribution."""

    samples: np.ndarray
    nominal: AgeResult
    n_failed: int = 0

    @property
    def mean(self) -> float:
        return float(np.mean(self.samples))

    @property
    def std(self) -> float:
        return float(np.std(self.samples, ddof=1))

    def interval(self, level: float = 0.68) -> tuple:
        a = (1 - level) / 2
        return tuple(np.quantile(self.samples, [a, 1 - a]))

    def summary(self) -> str:
        lo, hi = self.interval(0.68)
        lo95, hi95 = self.interval(0.95)
        return (
            f"Age (nominal) = {self.nominal.age:.4g} ka\n"
            f"MC: mean = {self.mean:.4g} ± {self.std:.2g} ka "
            f"(68 %: {lo:.4g}–{hi:.4g}; 95 %: {lo95:.4g}–{hi95:.4g}; n = {self.samples.size}"
            + (f", {self.n_failed} failed" if self.n_failed else "")
            + ")"
        )


@dataclass
class ToothSample:
    """ESR dating of tooth enamel with the classical component model.

    Concentrations: U in ppm. Dose in Gy, dose rates in Gy/ka, ages in ka.

    Parameters
    ----------
    De : equivalent dose of the enamel.
    enamel_U, dentine_U : present-day U in each tissue.
    sediment : U, Th, K and water content of the surrounding sediment.
    beta : geometry factors (see :class:`eprdating.beta.BetaGeometry`).
    gamma : external gamma dose rate. If None, computed from ``sediment``
        as an infinite matrix (use in-situ measurements when available).
    cosmic : cosmic dose rate (see :func:`eprdating.dose_rate.cosmic_dose_rate`).
    k_alpha : alpha efficiency of enamel.
    dentine_water : water content used for the dentine beta contribution.
    uptake_enamel, uptake_dentine : uptake models (EU, LU, US with p).
    useries : :class:`eprdating.series.USeries` to account for 230Th
        ingrowth; None assumes secular equilibrium (warns).
    factors : name of the conversion-factor set.
    """

    De: ValueLike
    enamel_U: ValueLike
    dentine_U: ValueLike
    sediment: Sediment
    beta: BetaGeometry
    cosmic: ValueLike
    gamma: ValueLike | None = None
    k_alpha: ValueLike = K_ENAMEL
    dentine_water: ValueLike = 0.0
    uptake_enamel: USModel = field(default_factory=lambda: USModel(-1.0))
    uptake_dentine: USModel = field(default_factory=lambda: USModel(-1.0))
    useries: USeries | None = None
    factors: str = DEFAULT_FACTORS

    # ---- parameter handling -------------------------------------------
    def _nominal(self) -> dict[str, float]:
        v = {
            "De": as_value(self.De).value,
            "enamel_U": as_value(self.enamel_U).value,
            "dentine_U": as_value(self.dentine_U).value,
            "cosmic": as_value(self.cosmic).value,
            "k_alpha": as_value(self.k_alpha).value,
            "dentine_water": as_value(self.dentine_water).value,
        }
        if self.gamma is not None:
            v["gamma"] = as_value(self.gamma).value
        for k in ("U", "Th", "K", "water"):
            v["sed_" + k] = as_value(getattr(self.sediment, k)).value
        for k, val in self.beta.values().items():
            v["beta_" + k] = val.value
        return v

    def _sample(self, rng: np.random.Generator, n: int) -> dict[str, np.ndarray]:
        s = {k: as_value(getattr(self, k)).sample(rng, n)
             for k in ("De", "enamel_U", "dentine_U", "cosmic", "k_alpha", "dentine_water")}
        if self.gamma is not None:
            s["gamma"] = as_value(self.gamma).sample(rng, n)
        for k in ("U", "Th", "K", "water"):
            s["sed_" + k] = as_value(getattr(self.sediment, k)).sample(rng, n)
        for k, val in self.beta.values().items():
            s["beta_" + k] = val.sample(rng, n)
        # physical constraints: no negative concentrations / fractions
        for k, arr in s.items():
            np.clip(arr, 0.0, None, out=arr)
        return s

    # ---- model -----------------------------------------------------------
    def components(self, v: dict[str, float] | None = None) -> list[DoseRateComponent]:
        v = v or self._nominal()
        cf: ConversionFactors = conversion_factors(self.factors)
        cU = {r: cf.get("U", r).value for r in ("alpha", "beta", "gamma")}
        Ga = self.useries.G("alpha") if self.useries else None
        Gb = self.useries.G("beta") if self.useries else None
        sed = {"U": v["sed_U"], "Th": v["sed_Th"], "K": v["sed_K"], "water": v["sed_water"]}
        sed_beta = water_correction(matrix_dose_rates(sed["U"], sed["Th"], sed["K"], cf)["beta"], sed["water"], "beta")
        if "gamma" in v:
            gamma = v["gamma"]
        else:
            gamma = water_correction(matrix_dose_rates(sed["U"], sed["Th"], sed["K"], cf)["gamma"], sed["water"], "gamma")
        dentine_beta = water_correction(v["dentine_U"] * cU["beta"], v["dentine_water"], "beta")
        return [
            DoseRateComponent("enamel alpha", v["k_alpha"] * v["enamel_U"] * cU["alpha"], self.uptake_enamel, Ga),
            DoseRateComponent("enamel beta", v["beta_internal"] * v["enamel_U"] * cU["beta"], self.uptake_enamel, Gb),
            DoseRateComponent("dentine beta", v["beta_dentine"] * dentine_beta, self.uptake_dentine, Gb),
            DoseRateComponent("sediment beta", v["beta_external"] * sed_beta),
            DoseRateComponent("gamma", gamma),
            DoseRateComponent("cosmic", v["cosmic"]),
        ]

    def age(self) -> AgeResult:
        """Nominal age with the central value of every input."""
        if self.useries is None and (as_value(self.enamel_U).value > 0 or as_value(self.dentine_U).value > 0):
            warnings.warn(
                "Secular equilibrium assumed for U in the dental tissues (useries=None). "
                "This ignores 230Th ingrowth and overestimates the internal dose rate, "
                "especially for samples younger than a few hundred ka.",
                stacklevel=2,
            )
        v = self._nominal()
        return solve_age(v["De"], self.components(v))

    def age_mc(self, n: int = 2000, seed: int | None = None) -> AgeMC:
        """Monte Carlo age: every input is sampled from its Gaussian uncertainty."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nominal = self.age()
        rng = np.random.default_rng(seed)
        s = self._sample(rng, n)
        ages = np.full(n, np.nan)
        for i in range(n):
            v = {k: float(arr[i]) for k, arr in s.items()}
            if v["De"] <= 0:
                continue
            try:
                ages[i] = solve_age(v["De"], self.components(v)).age
            except (RuntimeError, ValueError):
                pass
        ok = np.isfinite(ages)
        return AgeMC(samples=ages[ok], nominal=nominal, n_failed=int((~ok).sum()))
