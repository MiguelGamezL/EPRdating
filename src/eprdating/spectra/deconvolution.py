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
from dataclasses import dataclass, field
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
    errors: dict[str, float] = field(default_factory=dict)
    shift: float = 0.0  # field shift applied to all components (mT)
    noise_corr_length: int = 1  # points over which the residual noise is correlated

    def component(self, name: str) -> float:
        return self.amplitudes[name]


def _noise_autocov(resid: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Autocovariance of the residual noise, kept up to its first zero
    crossing (filtered cw-EPR noise is correlated over a few points by the
    time constant)."""
    r = resid - resid.mean()
    n = r.size
    max_lag = max_lag or max(1, n // 8)
    acov = np.array([np.dot(r[: n - k], r[k:]) / n for k in range(max_lag + 1)])
    neg = np.where(acov <= 0)[0]
    return acov[: int(neg[0])] if neg.size else acov


def _param_covariance(W: np.ndarray, acov: np.ndarray) -> np.ndarray:
    """``W C Wᵀ`` for a stationary noise covariance C given by ``acov``."""
    cov = acov[0] * W @ W.T
    for k in range(1, acov.size):
        t = W[:, :-k] @ W[:, k:].T
        cov += acov[k] * (t + t.T)
    return cov


class ComponentBasis:
    """A set of normalised component shapes on a fixed field grid."""

    def __init__(self, B, components: Mapping[str, Shape], baseline_order: int = 1):
        self.B = np.asarray(B, float)
        self.names = list(components)
        self._shapes = {}
        for k, s in components.items():
            if isinstance(s, tuple):
                self._shapes[k] = tuple(np.asarray(a, float) for a in s)
            else:
                self._shapes[k] = (self.B, np.asarray(s, float))
        self.matrix = self._matrix(0.0)
        self.baseline_order = baseline_order
        # baseline columns on a scaled field to keep the system well conditioned
        self._x = (self.B - self.B.mean()) / (np.ptp(self.B) or 1.0)
        self._poly = np.column_stack([self._x**k for k in range(baseline_order, -1, -1)]) if baseline_order >= 0 else None

    def _matrix(self, shift: float) -> np.ndarray:
        cols = []
        for (Bs, ys) in self._shapes.values():
            if shift == 0.0 and Bs is self.B:
                cols.append(normalise_pp(_on_grid(self.B, ys)))
            else:
                cols.append(normalise_pp(np.interp(self.B, Bs + shift, ys, left=0.0, right=0.0)))
        return np.column_stack(cols)

    def _solve(self, M: np.ndarray, y: np.ndarray, nonnegative: bool):
        A = M if self._poly is None else np.hstack([M, self._poly])
        nc = M.shape[1]
        if nonnegative:
            lb = np.r_[np.zeros(nc), np.full(A.shape[1] - nc, -np.inf)]
            x = lsq_linear(A, y, bounds=(lb, np.full(A.shape[1], np.inf)), method="bvls").x
        else:
            x, *_ = np.linalg.lstsq(A, y, rcond=None)
        return A, x

    def _best(self, y: np.ndarray, shifts: np.ndarray, nonnegative: bool):
        best = None
        for sh in shifts:
            A, x = self._solve(self._matrix(float(sh)), y, nonnegative)
            chi = float(np.sum((y - A @ x) ** 2))
            if best is None or chi < best[0]:
                best = (chi, float(sh), A, x)
        return best[1:]

    def fit(self, spectrum, nonnegative: bool = True, max_shift: float = 0.0,
            shift_step: float | None = None, noise=None, n_noise: int = 300,
            seed: int | None = 0) -> DeconvolutionResult:
        """Fit amplitudes (peak-to-peak units) and the polynomial baseline.

        nonnegative : constrain amplitudes to be >= 0. For dose-response work
            on weak signals use ``False``: the constraint biases noisy
            amplitudes upwards.
        max_shift : if > 0, a common field shift of all components within
            ``±max_shift`` mT is searched on a grid (step ``shift_step``,
            default one field step) and the best one kept.
        noise : optional signal-free stretch of the same spectrum (baseline
            corrected, same field step), at least as long as the fit window.
            The errors are then obtained by noise injection: ``n_noise``
            blocks of it are added to the fitted model and the fit, shift
            search included, is repeated. This captures the correlation of
            the noise and the extra scatter of the shift search on weak
            signals.

        Without ``noise`` the ``errors`` come from the linear least-squares
        weights and the noise covariance estimated from the residual
        autocovariance, at fixed shift. On short windows this underestimates
        the error by up to ~25 %, and more for weak signals when the shift
        is searched.
        """
        y = np.asarray(spectrum, float)
        if y.shape != self.B.shape:
            raise ValueError("spectrum must be sampled on the basis field grid")
        if max_shift > 0:
            step = shift_step or float(np.median(np.diff(self.B)))
            n = int(np.floor(max_shift / step + 1e-9))
            shifts = step * np.arange(-n, n + 1)  # symmetric, and includes zero
        else:
            shifts = np.array([0.0])
        shift, A, x = self._best(y, shifts, nonnegative)
        nc = len(self.names)
        fitted = A @ x
        resid = y - fitted
        ss = np.sum((y - y.mean()) ** 2)
        acov = _noise_autocov(resid)
        if noise is None:
            err = np.sqrt(np.diag(_param_covariance(np.linalg.pinv(A), acov)))[:nc]
        else:
            noise = np.asarray(noise, float)
            if noise.size < y.size:
                raise ValueError("noise stretch must be at least as long as the spectrum")
            rng = np.random.default_rng(seed)
            draws = np.empty((n_noise, nc))
            for k in range(n_noise):
                i = rng.integers(0, noise.size - y.size + 1)
                blk = noise[i:i + y.size]
                draws[k] = self._best(fitted + (blk - blk.mean()), shifts, nonnegative)[2][:nc]
            err = draws.std(axis=0, ddof=1)
        return DeconvolutionResult(
            amplitudes=dict(zip(self.names, map(float, x[:nc]))),
            baseline=x[nc:],
            fitted=fitted,
            residuals=resid,
            r2=float(1 - np.sum(resid**2) / ss) if ss > 0 else float("nan"),
            errors=dict(zip(self.names, map(float, err))),
            shift=shift,
            noise_corr_length=int(acov.size),
        )


def component_vs_dose(
    basis: ComponentBasis,
    spectra: Sequence[np.ndarray],
    component: str,
    **fit_kw,
) -> tuple[np.ndarray, np.ndarray, list]:
    """Amplitude of ``component`` in each spectrum of a dose series.

    Returns ``(amplitudes, errors, results)``; feed amplitudes and errors
    with their doses to :func:`eprdating.doseresponse.fit_dose_response`.
    ``fit_kw`` go to :meth:`ComponentBasis.fit`.
    """
    results = [basis.fit(s, **fit_kw) for s in spectra]
    amps = np.array([r.amplitudes[component] for r in results])
    errs = np.array([r.errors[component] for r in results])
    return amps, errs, results
