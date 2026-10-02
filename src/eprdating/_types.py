"""Small shared types: values with uncertainty."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np


@dataclass(frozen=True)
class Value:
    """A measured quantity with a 1-sigma (Gaussian) uncertainty."""

    value: float
    sigma: float = 0.0

    def __post_init__(self) -> None:
        if self.sigma < 0:
            raise ValueError(f"sigma must be >= 0, got {self.sigma}")

    def sample(self, rng: np.random.Generator, n: int) -> np.ndarray:
        """Draw ``n`` Gaussian samples (constant array if ``sigma == 0``)."""
        if self.sigma == 0:
            return np.full(n, float(self.value))
        return rng.normal(self.value, self.sigma, n)

    @property
    def rel(self) -> float:
        """Relative uncertainty sigma/|value| (``inf`` if value is 0)."""
        return self.sigma / abs(self.value) if self.value else float("inf")

    def __float__(self) -> float:  # allows float(Value(...))
        return float(self.value)

    def __repr__(self) -> str:
        return f"{self.value:.6g} ± {self.sigma:.2g}"


ValueLike = Union[Value, float, int, tuple]  # noqa: UP007 - runtime alias, py3.10 compatible


def as_value(x: ValueLike) -> Value:
    """Coerce ``x`` (a :class:`Value`, a number or a ``(value, sigma)`` tuple)."""
    if isinstance(x, Value):
        return x
    if isinstance(x, tuple):
        if len(x) != 2:
            raise ValueError("tuple must be (value, sigma)")
        return Value(float(x[0]), float(x[1]))
    return Value(float(x), 0.0)
