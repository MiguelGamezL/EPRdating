"""Uranium-uptake models for dental tissues.

The US model (Grün et al. 1988) describes the U content of a tissue at time
``t`` after burial (``0 <= t <= T``, with ``T`` the age) as

    U(t) = U_m * (t / T) ** (p + 1),       p >= -1

where ``U_m`` is the U concentration measured today. Special cases:

* ``p = -1``  early uptake (EU): all U present since burial.
* ``p =  0``  linear uptake (LU).
* ``p >  0``  increasingly recent uptake.

Uptake only affects the *U-derived* dose-rate contributions of the tissue
(alpha and beta of enamel, beta of dentine); sediment gamma and cosmic
contributions are taken as constant.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad


@dataclass(frozen=True)
class USModel:
    """Uptake following ``(t/T)**(p+1)``; ``p=-1`` is EU, ``p=0`` is LU."""

    p: float = 0.0

    def __post_init__(self) -> None:
        if self.p < -1:
            raise ValueError("p must be >= -1")

    @property
    def is_early(self) -> bool:
        return self.p == -1

    def fraction(self, t, T: float):
        """U(t)/U_m at time ``t`` after burial for a sample of age ``T``."""
        x = np.clip(np.asarray(t, float) / T, 0.0, 1.0)
        return x ** (self.p + 1)

    def accumulated(self, T: float, G=None) -> float:
        """Dose accumulated over ``T`` per unit present-day dose rate.

        ``G(tau)`` is the time-integrated activity ratio of U that has been in
        the tissue for a time ``tau`` (``G(tau) = tau`` in secular
        equilibrium, see :mod:`eprdating.series`). The result is

            A(T) = ∫ G(T - t') dU(t')/U_m

        which reduces to ``T / (p + 2)`` in equilibrium.
        """
        if G is None:
            return T / (self.p + 2.0)
        if self.is_early:
            return float(G(T))
        # dU/U_m = (p+1) x**p dx with x = t'/T; algebraic weight handles x**p
        val, _ = quad(lambda x: G(T * (1.0 - x)), 0.0, 1.0, weight="alg", wvar=(self.p, 0.0), limit=200)
        return (self.p + 1.0) * val


def EarlyUptake() -> USModel:
    return USModel(-1.0)


def LinearUptake() -> USModel:
    return USModel(0.0)


@dataclass(frozen=True)
class DelayedUptake:
    """All U taken up at once ``t_uptake`` ka before present (none before).

    Used by the CSUS-ESR model (Grün 2000), where ``t_uptake`` is the
    closed-system U-series age of the tissue. For a sample older than
    ``t_uptake`` the tissue delivers dose only during its last ``t_uptake``.
    """

    t_uptake: float

    def __post_init__(self) -> None:
        if self.t_uptake < 0:
            raise ValueError("t_uptake must be >= 0")

    def fraction(self, t, T: float):
        """U(t)/U_m at time ``t`` after burial for a sample of age ``T``."""
        return (np.asarray(t, float) >= T - self.t_uptake).astype(float)

    def accumulated(self, T: float, G=None) -> float:
        tau = min(self.t_uptake, T)
        return float(tau if G is None else G(tau))
