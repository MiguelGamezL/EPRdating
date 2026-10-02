"""Combined U-series / ESR (US-ESR) dating of teeth.

Grün, Schwarcz & Chadam (1988): the uptake parameter ``p`` of each dental
tissue is not assumed (EU, LU) but solved together with the age from the
tissue's measured U-series activity ratios.

Model
-----
Uranium accumulates as ``U(t) = U_m (t/T)^(p+1)`` (:class:`~eprdating.uptake.USModel`).
Every parcel of uranium arrives with the same 234U/238U ratio ``r_in`` (that of
the groundwater) and without 230Th. A parcel that has spent a time ``τ`` in the
tissue shows today

    234U/238U = 1 + (r_in - 1) e^{-λ234 τ}
    230Th/238U = 1 - e^{-λ230 τ} + (r_in - 1) λ230/(λ230 - λ234) (e^{-λ234 τ} - e^{-λ230 τ})

and the tissue ratios are averages over parcels weighted by ``dU``. For given
``T`` and ``p`` the measured 234U/238U fixes ``r_in`` (linear relation), and
the measured 230Th/234U then fixes ``p``. The age is the ``T`` at which the
ESR dose equation is also satisfied, with each tissue's ``p(T)`` and
``r_in(T)``.

The U-series data alone set a lower bound on the age: the closed-system
(early uptake) U-series age. When the ESR dose is already exceeded at that
bound there is no US-ESR solution (typically uranium leaching); this is
reported, not forced.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from ._types import ValueLike, as_value
from .age import AgeResult, ToothSample, solve_age
from .series import LAMBDA, activity_ratio_Th230_U238
from .uptake import USModel

P_MAX = 1000.0


def _expect(f, T: float, p: float) -> float:
    """Average of ``f(τ)`` over the uranium parcels of a tissue of age T."""
    if p == -1.0:
        return f(T)
    val, _ = quad(lambda x: f(T * (1.0 - x)), 0.0, 1.0, weight="alg", wvar=(p, 0.0), limit=200)
    return (p + 1.0) * val


def incoming_ratio(T: float, p: float, u234_u238: float) -> float:
    """234U/238U of the incoming uranium that reproduces the measured ratio."""
    m = _expect(lambda tau: math.exp(-LAMBDA["U234"] * tau), T, p)
    return 1.0 + (u234_u238 - 1.0) / m


def predicted_ratios(T: float, p: float, r_in: float) -> tuple[float, float]:
    """Present-day (234U/238U, 230Th/234U) of a tissue for age T, uptake p and
    incoming ratio r_in."""
    l4 = LAMBDA["U234"]
    r48 = 1.0 + (r_in - 1.0) * _expect(lambda tau: math.exp(-l4 * tau), T, p)
    r08 = _expect(lambda tau: activity_ratio_Th230_U238(tau, r_in), T, p)
    return r48, r08 / r48


def th230_u234(T: float, p: float, u234_u238: float) -> float:
    """230Th/234U predicted for a tissue with the measured 234U/238U."""
    return predicted_ratios(T, p, incoming_ratio(T, p, u234_u238))[1]


def closed_system_age(th230_u234_meas: float, u234_u238: float, t_max: float = 5000.0) -> float:
    """Closed-system (early uptake) U-series age in ka; ``inf`` if beyond range."""
    if th230_u234_meas <= 0:
        return 0.0
    f = lambda T: th230_u234(T, -1.0, u234_u238) - th230_u234_meas
    if f(t_max) < 0:
        return math.inf
    return brentq(f, 1e-6, t_max, xtol=1e-9)


def solve_p(T: float, th230_u234_meas: float, u234_u238: float) -> float | None:
    """Uptake parameter p of a tissue for age T, or None if T is younger than
    the tissue's closed-system U-series age (no p >= -1 fits)."""
    g = lambda p: th230_u234(T, p, u234_u238) - th230_u234_meas
    g_eu = g(-1.0)
    if g_eu < 0:
        return None
    if abs(g_eu) < 1e-12:
        return -1.0
    if g(P_MAX) > 0:
        return P_MAX
    return brentq(g, -1.0, P_MAX, xtol=1e-10)


@dataclass
class UseriesData:
    """Measured U-series activity ratios of one dental tissue."""

    th230_u234: ValueLike
    u234_u238: ValueLike


@dataclass
class USESRResult:
    age: float | None  # ka
    p_enamel: float | None
    p_dentine: float | None
    r_in_enamel: float | None
    r_in_dentine: float | None
    min_age: float  # closed-system U-series bound, ka
    status: str  # "ok" or "no_solution"
    detail: AgeResult | None = None

    def summary(self) -> str:
        if self.status != "ok":
            return (f"No US-ESR solution: the ESR dose is reached before the closed-system "
                    f"U-series age ({self.min_age:.4g} ka); uranium leaching is likely.")
        fmt = lambda x: "—" if x is None else f"{x:.3g}"
        return (f"US-ESR age = {self.age:.4g} ka   p(enamel) = {fmt(self.p_enamel)}, "
                f"p(dentine) = {fmt(self.p_dentine)}   (U-series lower bound {self.min_age:.4g} ka)")


@dataclass
class USESRMC:
    ages: np.ndarray
    p_enamel: np.ndarray
    p_dentine: np.ndarray
    nominal: USESRResult
    n_failed: int = 0

    @property
    def mean(self) -> float:
        return float(np.mean(self.ages))

    @property
    def std(self) -> float:
        return float(np.std(self.ages, ddof=1))

    def summary(self) -> str:
        lo, hi = np.quantile(self.ages, [0.16, 0.84])
        return (f"US-ESR MC: {self.mean:.4g} ± {self.std:.2g} ka (68 %: {lo:.4g}–{hi:.4g}; "
                f"n = {self.ages.size}, failed {self.n_failed})")


@dataclass
class USESRSample:
    """A tooth for US-ESR dating.

    ``tooth`` carries everything except the uptake history (its uptake models
    and 234U/238U fields are ignored for tissues with U-series data).
    """

    tooth: ToothSample
    enamel: UseriesData | None = None
    dentine: UseriesData | None = None

    _tissues = ("enamel", "dentine")

    # ---- inputs --------------------------------------------------------------
    def _has(self, tissue: str, v: dict) -> bool:
        return getattr(self, tissue) is not None and v[f"{tissue}_U"] > 0

    def _nominal(self) -> dict:
        v = self.tooth._nominal()
        for t in self._tissues:
            d = getattr(self, t)
            if d is not None:
                v[f"us_{t}_th"] = as_value(d.th230_u234).value
                v[f"us_{t}_r"] = as_value(d.u234_u238).value
        return v

    def _sample(self, rng: np.random.Generator, n: int) -> dict:
        s = self.tooth._sample(rng, n)
        for t in self._tissues:
            d = getattr(self, t)
            if d is not None:
                s[f"us_{t}_th"] = np.clip(as_value(d.th230_u234).sample(rng, n), 0.0, None)
                s[f"us_{t}_r"] = np.clip(as_value(d.u234_u238).sample(rng, n), 0.0, None)
        return s

    # ---- solver --------------------------------------------------------------
    def _state(self, T: float, v: dict):
        """Uptake models and incoming ratios of both tissues at age T."""
        ups, rins = {}, {}
        for t in self._tissues:
            if self._has(t, v):
                p = solve_p(T, v[f"us_{t}_th"], v[f"us_{t}_r"])
                if p is None:
                    return None
                ups[t] = USModel(p)
                rins[t] = incoming_ratio(T, p, v[f"us_{t}_r"])
            else:
                ups[t] = getattr(self.tooth, f"uptake_{t}")
                rins[t] = None
        return ups, rins

    def _dose(self, T: float, v: dict) -> float:
        st = self._state(T, v)
        if st is None:
            raise ValueError("age below the closed-system U-series age")
        ups, rins = st
        vv = dict(v)
        for t in self._tissues:
            if rins[t] is not None:
                vv[f"u234_u238_{t}"] = rins[t]
        comps = self._components(vv, ups)
        return sum(c.accumulated(T) for c in comps)

    def _components(self, v: dict, ups: dict):
        tooth = self.tooth
        saved = tooth.u234_u238_is
        tooth.u234_u238_is = "initial"
        try:
            return tooth.components(v, uptake_enamel=ups["enamel"], uptake_dentine=ups["dentine"])
        finally:
            tooth.u234_u238_is = saved

    def _solve(self, v: dict) -> USESRResult:
        bounds = [closed_system_age(v[f"us_{t}_th"], v[f"us_{t}_r"]) for t in self._tissues if self._has(t, v)]
        t_min = max(bounds) if bounds else 1e-6
        if not math.isfinite(t_min):
            return USESRResult(None, None, None, None, None, t_min, "no_solution")
        t_min = max(t_min * (1 + 1e-9), 1e-6)
        De = v["De"]
        d_min = self._dose(t_min, v)
        if d_min > De * (1 + 1e-6):
            return USESRResult(None, None, None, None, None, t_min, "no_solution")
        if d_min >= De:  # solution at the bound: early uptake in the limiting tissue
            T = t_min
        else:
            hi = max(2 * t_min, 1.0)
            while self._dose(hi, v) < De:
                hi *= 2
                if hi > 1e5:
                    raise RuntimeError("no age below 100 Ma reproduces De")
            T = brentq(lambda x: self._dose(x, v) - De, t_min, hi, xtol=1e-8, rtol=1e-10)
        ups, rins = self._state(T, v)
        vv = dict(v)
        for t in self._tissues:
            if rins[t] is not None:
                vv[f"u234_u238_{t}"] = rins[t]
        detail = solve_age(De, self._components(vv, ups))
        pe = ups["enamel"].p if self._has("enamel", v) else None
        pd = ups["dentine"].p if self._has("dentine", v) else None
        return USESRResult(T, pe, pd, rins["enamel"], rins["dentine"], t_min, "ok", detail)

    def age(self) -> USESRResult:
        """Nominal US-ESR age and uptake parameters."""
        return self._solve(self._nominal())

    def age_mc(self, n: int = 1000, seed: int | None = None) -> USESRMC:
        """Monte Carlo over all inputs, U-series ratios included."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nominal = self.age()
        rng = np.random.default_rng(seed)
        s = self._sample(rng, n)
        ages, pe, pd = [], [], []
        failed = 0
        for i in range(n):
            v = {k: float(a[i]) for k, a in s.items()}
            try:
                r = self._solve(v) if v["De"] > 0 else None
            except (RuntimeError, ValueError):
                r = None
            if r is None or r.status != "ok":
                failed += 1
                continue
            ages.append(r.age)
            pe.append(np.nan if r.p_enamel is None else r.p_enamel)
            pd.append(np.nan if r.p_dentine is None else r.p_dentine)
        return USESRMC(np.array(ages), np.array(pe), np.array(pd), nominal, failed)


__all__ = [
    "USESRMC",
    "USESRResult",
    "USESRSample",
    "UseriesData",
    "closed_system_age",
    "incoming_ratio",
    "predicted_ratios",
    "solve_p",
    "th230_u234",
]

