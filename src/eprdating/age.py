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
from .onegroup import ToothLayers
from .series import USeries, load_partition
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
    beta : either fixed geometry factors (:class:`eprdating.beta.BetaGeometry`)
        or a layered geometry (:class:`eprdating.onegroup.ToothLayers`), in
        which case beta attenuation is computed with one-group theory for each
        emitter and U-series segment (as ROSY does). Uncertainties on the
        layer thicknesses, stripping and densities are sampled in
        :meth:`age_mc` together with the water contents.
    gamma : external gamma dose rate. If None, computed from ``sediment``
        as an infinite matrix (use in-situ measurements when available).
    cosmic : cosmic dose rate (see :func:`eprdating.dose_rate.cosmic_dose_rate`).
    k_alpha : alpha efficiency of enamel.
    dentine_water : water content used for the dentine beta contribution.
    uptake_enamel, uptake_dentine : uptake models (EU, LU, US with p).
    u234_u238_enamel, u234_u238_dentine : present-day (measured) 234U/238U
        activity ratio of each tissue.
    radon_loss_enamel, radon_loss_dentine : fraction of 222Rn escaping each
        tissue (0 to 1).
    ingrowth : model U-series daughter ingrowth after uptake (default). With
        ``False`` the tissues are taken in secular equilibrium (warns).
    partition : optional custom U-series segment table
        (see :mod:`eprdating.series`).
    factors : name of the conversion-factor set.
    sample_geometry : with a :class:`~eprdating.onegroup.ToothLayers`
        geometry, recompute the one-group factors for every Monte Carlo draw
        (default). With ``False`` they are kept at their nominal values.
    """

    De: ValueLike
    enamel_U: ValueLike
    dentine_U: ValueLike
    sediment: Sediment
    beta: BetaGeometry | ToothLayers
    cosmic: ValueLike
    gamma: ValueLike | None = None
    k_alpha: ValueLike = K_ENAMEL
    dentine_water: ValueLike = 0.0
    uptake_enamel: USModel = field(default_factory=lambda: USModel(-1.0))
    uptake_dentine: USModel = field(default_factory=lambda: USModel(-1.0))
    u234_u238_enamel: ValueLike = 1.0
    u234_u238_dentine: ValueLike = 1.0
    radon_loss_enamel: ValueLike = 0.0
    radon_loss_dentine: ValueLike = 0.0
    ingrowth: bool = True
    partition: dict | None = None
    factors: str = DEFAULT_FACTORS
    sample_geometry: bool = True

    _SCALARS = (
        "De", "enamel_U", "dentine_U", "cosmic", "k_alpha", "dentine_water",
        "u234_u238_enamel", "u234_u238_dentine", "radon_loss_enamel", "radon_loss_dentine",
    )

    # ---- parameter handling -------------------------------------------
    def _nominal(self) -> dict[str, float]:
        v = {k: as_value(getattr(self, k)).value for k in self._SCALARS}
        if self.gamma is not None:
            v["gamma"] = as_value(self.gamma).value
        for k in ("U", "Th", "K", "water"):
            v["sed_" + k] = as_value(getattr(self.sediment, k)).value
        if isinstance(self.beta, BetaGeometry):
            for k, val in self.beta.values().items():
                v["beta_" + k] = val.value
        else:
            for k, val in self.beta.nominal_values().items():
                v["geo_" + k] = val
        return v

    def _sample(self, rng: np.random.Generator, n: int) -> dict[str, np.ndarray]:
        s = {k: as_value(getattr(self, k)).sample(rng, n) for k in self._SCALARS}
        if self.gamma is not None:
            s["gamma"] = as_value(self.gamma).sample(rng, n)
        for k in ("U", "Th", "K", "water"):
            s["sed_" + k] = as_value(getattr(self.sediment, k)).sample(rng, n)
        if isinstance(self.beta, BetaGeometry):
            for k, val in self.beta.values().items():
                s["beta_" + k] = val.sample(rng, n)
        elif self.sample_geometry:
            for k in self.beta.GEOMETRY:
                s["geo_" + k] = as_value(getattr(self.beta, k)).sample(rng, n)
        # physical constraints: no negative concentrations / fractions
        for k, arr in s.items():
            np.clip(arr, 0.0, None, out=arr)
        for k in ("radon_loss_enamel", "radon_loss_dentine"):
            np.clip(s[k], 0.0, 1.0, out=s[k])
        return s

    def _useries(self, tissue: str, v: dict[str, float]) -> USeries | None:
        if not self.ingrowth:
            return None
        return USeries(
            ratio=v[f"u234_u238_{tissue}"],
            ratio_is="present",
            radon_loss=v[f"radon_loss_{tissue}"],
            partition=self._partition(),
        )

    def _partition(self) -> dict:
        if self.partition is None:
            self.partition = load_partition()
        return self.partition

    # ---- one-group beta factors ------------------------------------------
    def _onegroup(self, v: dict[str, float]) -> dict:
        """Per-source, per-chain/segment one-group factors for the geometry and
        water contents in ``v`` (nominal geometry if ``v`` has none)."""
        geo_vals = {k: v.get("geo_" + k, as_value(getattr(self.beta, k)).value) for k in self.beta.GEOMETRY}
        if not self.sample_geometry:
            geo_vals = self.beta.nominal_values()
            water = (as_value(self.sediment.water).value, as_value(self.dentine_water).value)
        else:
            water = (v["sed_water"], v["dentine_water"])
        key = (tuple(geo_vals.values()), water)
        cache = self.__dict__.setdefault("_og_cache", {})
        if key not in cache:
            geo = self.beta.at(**geo_vals, sediment_water=water[0], dentine_water=water[1])
            segs = ("U238", "U234", "Th230", "Rn222", "U235", "Pa231")
            if len(cache) > 8:  # keep the nominal entries, not every MC draw
                cache.clear()
            cache[key] = {
                "enamel": {s: geo.chain_fraction("enamel", s) for s in segs},
                "dentine": {s: geo.chain_fraction("dentine", s) for s in segs},
                "enamel_U": geo.chain_fraction("enamel", "U"),
                "dentine_U": geo.chain_fraction("dentine", "U"),
                "sediment": {c: geo.chain_fraction("sediment", c) for c in ("U", "Th", "K")},
            }
        return cache[key]

    # ---- model -----------------------------------------------------------
    def components(self, v: dict[str, float] | None = None) -> list[DoseRateComponent]:
        v = v or self._nominal()
        cf: ConversionFactors = conversion_factors(self.factors)
        cU = {r: cf.get("U", r).value for r in ("alpha", "beta", "gamma")}
        use, usd = self._useries("enamel", v), self._useries("dentine", v)
        sed = {"U": v["sed_U"], "Th": v["sed_Th"], "K": v["sed_K"], "water": v["sed_water"]}
        dry = matrix_dose_rates(sed["U"], sed["Th"], sed["K"], cf)
        if "gamma" in v:
            gamma = v["gamma"]
        else:
            gamma = water_correction(dry["gamma"], sed["water"], "gamma")
        Ga_e = use.G("alpha") if use else None

        if isinstance(self.beta, BetaGeometry):
            Gb_e = use.G("beta") if use else None
            Gb_d = usd.G("beta") if usd else None
            en_beta = v["beta_internal"] * v["enamel_U"] * cU["beta"]
            den_beta = v["beta_dentine"] * water_correction(v["dentine_U"] * cU["beta"], v["dentine_water"], "beta")
            sed_beta = v["beta_external"] * water_correction(dry["beta"], sed["water"], "beta")
        else:
            og = self._onegroup(v)
            # one-group factors are relative to the wet medium's own infinite-matrix
            # dose, i.e. the dry dose rate diluted by (1 + water)
            # rates are the attenuated equilibrium values; the per-segment factors
            # enter G relative to the whole-chain factor
            w_e = {s: f / og["enamel_U"] for s, f in og["enamel"].items()}
            w_d = {s: f / og["dentine_U"] for s, f in og["dentine"].items()}
            Gb_e = use.G("beta", w_e) if use else None
            Gb_d = usd.G("beta", w_d) if usd else None
            en_beta = v["enamel_U"] * cU["beta"] * og["enamel_U"]
            den_beta = v["dentine_U"] * cU["beta"] / (1.0 + v["dentine_water"]) * og["dentine_U"]
            sed_beta = sum(
                sed[nuc] * cf.get(nuc, "beta").value * og["sediment"][nuc] for nuc in ("U", "Th", "K")
            ) / (1.0 + sed["water"])
        return [
            DoseRateComponent("enamel alpha", v["k_alpha"] * v["enamel_U"] * cU["alpha"], self.uptake_enamel, Ga_e),
            DoseRateComponent("enamel beta", en_beta, self.uptake_enamel, Gb_e),
            DoseRateComponent("dentine beta", den_beta, self.uptake_dentine, Gb_d),
            DoseRateComponent("sediment beta", sed_beta),
            DoseRateComponent("gamma", gamma),
            DoseRateComponent("cosmic", v["cosmic"]),
        ]

    def age(self) -> AgeResult:
        """Nominal age with the central value of every input."""
        if not self.ingrowth and (as_value(self.enamel_U).value > 0 or as_value(self.dentine_U).value > 0):
            warnings.warn(
                "Secular equilibrium assumed for U in the dental tissues (ingrowth=False). "
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
