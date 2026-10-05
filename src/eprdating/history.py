"""Piecewise-constant histories of environmental parameters.

The burial environment of a tooth need not stay as it is today: a wetter
period raises the sediment water and lowers the external dose rate, and
erosion or deposition changes the overburden and so the cosmic dose rate.
A :class:`History` gives one parameter as a sequence of constant values,
from the present backwards in time::

    # 10 % water today and back to 12 ka, 25 % from 12 to 60 ka, 15 % before
    water = History([(0.10, 0.03), (0.25, 0.05), (0.15, 0.05)], breaks=[12, 60])

Times are in ka before present and the last value holds back to any age.
In the Monte Carlo every segment value is sampled independently; the break
times are fixed.

Histories are accepted by :class:`eprdating.age.ToothSample` for the
sediment water (``Sediment(water=History(...))``), the cosmic dose rate and
the external gamma dose rate; :func:`eprdating.dose_rate.cosmic_history`
turns a burial-depth history into a cosmic dose-rate history. The internal
components (enamel, dentine and cementum) keep the present-day values.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise

import numpy as np

from ._types import Value, ValueLike, as_value


@dataclass(frozen=True, init=False)
class History:
    """Piecewise-constant parameter history, in ka before present.

    ``values[0]`` holds from today back to ``breaks[0]``, ``values[i]``
    between ``breaks[i-1]`` and ``breaks[i]``, and the last value before
    the last break.
    """

    values: tuple[Value, ...]
    breaks: tuple[float, ...]

    def __init__(self, values: Sequence[ValueLike], breaks: Sequence[float] = ()) -> None:
        vals = tuple(as_value(x) for x in values)
        brk = tuple(float(b) for b in breaks)
        if len(vals) != len(brk) + 1:
            raise ValueError(f"a history with {len(brk)} breaks needs {len(brk) + 1} values, got {len(vals)}")
        if any(b <= 0 or not math.isfinite(b) for b in brk) or any(b2 <= b1 for b1, b2 in pairwise(brk)):
            raise ValueError("breaks must be positive and strictly increasing (ka before present)")
        object.__setattr__(self, "values", vals)
        object.__setattr__(self, "breaks", brk)

    @classmethod
    def constant(cls, value: ValueLike) -> History:
        return cls([value])

    def __len__(self) -> int:
        return len(self.values)

    def nominal(self) -> list[float]:
        return [v.value for v in self.values]

    def sample(self, rng: np.random.Generator, n: int) -> list[np.ndarray]:
        """``n`` draws of every segment value, sampled independently."""
        return [v.sample(rng, n) for v in self.values]

    def at(self, t: float) -> float:
        """Nominal value at ``t`` ka before present."""
        return self.nominal()[int(np.searchsorted(self.breaks, t, side="right"))]

    def mean(self, T: float) -> float:
        """Time-weighted nominal mean over the last ``T`` ka."""
        return integrate(self.breaks, self.nominal(), T) / T


def integrate(breaks: Sequence[float], rates: Sequence[float], T: float) -> float:
    """∫_0^T of a piecewise-constant function (``rates`` between ``breaks``)."""
    total, lo = 0.0, 0.0
    for r, hi in zip(rates, (*breaks, math.inf), strict=True):
        if T <= lo:
            break
        total += r * (min(T, hi) - lo)
        lo = hi
    return total
