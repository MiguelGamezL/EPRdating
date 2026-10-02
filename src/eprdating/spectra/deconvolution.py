"""Spectral deconvolution into known components.

The measured derivative spectrum is modelled as a non-negative combination
of component line shapes (e.g. the axial and orthorhombic CO2- species, the
native organic signal…) plus a free polynomial baseline:

    y(B) ≈ Σ_i a_i · s_i(B) + Σ_k c_k · B^k ,   a_i >= 0

The component shapes ``s_i`` can come from any source: reference spectra,
analytic line shapes, or spin-Hamiltonian simulations with EPRAYA
(:mod:`eprdating.spectra.epraya_backend`). Each shape is normalised to unit
peak-to-peak so the amplitudes ``a_i`` are peak-to-peak intensities.

Fitting the amplitudes is a bounded *linear* problem, so it is fast and has
a unique solution; the (slow) simulation of the shapes is done once.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Union

import numpy as np
from scipy.optimize import lsq_linear

Shape = Union[np.ndarray, tuple[np.ndarray, np.ndarray]]  # noqa: UP007 - runtime alias


def _on_grid(B: np.ndarray, shape: Shape) -> np.ndarray:
    if isinstance(shape, tuple):
        Bs, ys = (np.asarray(a, float) for a in shape)
        return np.interp(B, Bs, ys, left=0.0, right=0.0)
    ys = np.asarray(shape, float)
    if ys.shape != B.shape:
        raise ValueError("component array must match the field grid; pass (B, y) to interpolate")
    return ys


def normalise_pp(y: np.ndarray) -> np.ndarray:
    pp = y.max() - y.min()
    if pp == 0:
        raise ValueError("component shape is flat")
    return y / pp


@dataclass
class DeconvolutionResult:
    amplitudes: dict[str, float]
    baseline: np.ndarray  # polynomial coefficients, highest power first
    fitted: np.ndarray
    residuals: np.ndarray
    r2: float

    def component(self, name: str) -> float:
        return self.amplitudes[name]


class ComponentBasis:
    """A set of normalised component shapes on a fixed field grid."""

    def __init__(self, B, components: Mapping[str, Shape], baseline_order: int = 1):
        self.B = np.asarray(B, float)
        self.names = list(components)
        self.matrix = np.column_stack([normalise_pp(_on_grid(self.B, s)) for s in components.values()])
        self.baseline_order = baseline_order
        # baseline columns on a scaled field to keep the system well conditioned
        self._x = (self.B - self.B.mean()) / (np.ptp(self.B) or 1.0)
        self._poly = np.column_stack([self._x**k for k in range(baseline_order, -1, -1)]) if baseline_order >= 0 else None

    def fit(self, spectrum) -> DeconvolutionResult:
        y = np.asarray(spectrum, float)
        if y.shape != self.B.shape:
            raise ValueError("spectrum must be sampled on the basis field grid")
        A = self.matrix if self._poly is None else np.hstack([self.matrix, self._poly])
        nc = self.matrix.shape[1]
        lb = np.r_[np.zeros(nc), np.full(A.shape[1] - nc, -np.inf)]
        ub = np.full(A.shape[1], np.inf)
        sol = lsq_linear(A, y, bounds=(lb, ub), method="bvls")
        fitted = A @ sol.x
        resid = y - fitted
        ss = np.sum((y - y.mean()) ** 2)
        return DeconvolutionResult(
            amplitudes=dict(zip(self.names, map(float, sol.x[:nc]))),
            baseline=sol.x[nc:],
            fitted=fitted,
            residuals=resid,
            r2=float(1 - np.sum(resid**2) / ss) if ss > 0 else float("nan"),
        )


def component_vs_dose(
    basis: ComponentBasis,
    spectra: Sequence[np.ndarray],
    component: str,
) -> tuple[np.ndarray, list]:
    """Amplitude of ``component`` in each spectrum of a dose series.

    Returns ``(amplitudes, results)``; feed the amplitudes with their doses
    to :func:`eprdating.doseresponse.fit_dose_response`.
    """
    results = [basis.fit(s) for s in spectra]
    return np.array([r.amplitudes[component] for r in results]), results
