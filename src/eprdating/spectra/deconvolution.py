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

    def _best_many(self, Y: np.ndarray, shifts: np.ndarray, nonnegative: bool) -> np.ndarray:
        """Coefficients of the best shift for each row of ``Y``. Unconstrained
        fits are done for all rows at once, one pseudo-inverse per shift."""
        if nonnegative:
            return np.array([self._best(y, shifts, True)[2] for y in Y])
        best_chi = np.full(len(Y), np.inf)
        best_x = None
        for sh in shifts:
            M = self._matrix(float(sh))
            A = M if self._poly is None else np.hstack([M, self._poly])
            X = np.linalg.pinv(A) @ Y.T  # (n_params, n_rows)
            chi = np.sum((Y.T - A @ X) ** 2, axis=0)
            if best_x is None:
                best_x = X.T.copy()
                best_chi = chi
            else:
                better = chi < best_chi
                best_x[better] = X.T[better]
                best_chi = np.where(better, chi, best_chi)
        return best_x

    def _shift_grid(self, max_shift: float, shift_step: float | None) -> np.ndarray:
        if max_shift <= 0:
            return np.array([0.0])
        step = shift_step or float(np.median(np.diff(self.B)))
        n = int(np.floor(max_shift / step + 1e-9))
        return step * np.arange(-n, n + 1)  # symmetric, and includes zero

    def null_amplitudes(self, noise, n_noise: int = 300, max_shift: float = 0.0,
                        shift_step: float | None = None, nonnegative: bool = False,
                        seed: int | None = 0) -> np.ndarray:
        """Amplitudes fitted to signal-free noise alone, ``(n_noise, n_components)``.

        ``noise`` is either a signal-free stretch, from which ``n_noise``
        blocks are drawn (as in :meth:`fit`), or a 2-D array of noise
        realisations on this grid (one per row), each fitted with the same
        shift search. Their distribution is what the fit returns when there is no
        signal: with a shift search it is not centred on zero for one sign,
        because the search finds the place where the noise looks most like
        the component. Compare a measured amplitude with it to decide whether
        a signal is detected.
        """
        noise = np.asarray(noise, float)
        shifts = self._shift_grid(max_shift, shift_step)
        nc = len(self.names)
        if noise.ndim == 2:
            if noise.shape[1] != self.B.size:
                raise ValueError("noise realisations must be sampled on the basis grid")
            blocks = noise
        else:
            if noise.size < self.B.size:
                raise ValueError("noise stretch must be at least as long as the spectrum")
            rng = np.random.default_rng(seed)
            starts = rng.integers(0, noise.size - self.B.size + 1, n_noise)
            blocks = np.array([noise[i:i + self.B.size] for i in starts])
        blocks = blocks - blocks.mean(axis=1, keepdims=True)
        return self._best_many(blocks, shifts, nonnegative)[:, :nc]

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
            corrected, same field step), at least as long as the fit window,
            or a 2-D array of noise realisations on the basis grid.
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
        shifts = self._shift_grid(max_shift, shift_step)
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
            if noise.ndim == 2:  # ready-made noise realisations
                if noise.shape[1] != y.size:
                    raise ValueError("noise realisations must be sampled on the basis grid")
                blocks = noise
            else:
                if noise.size < y.size:
                    raise ValueError("noise stretch must be at least as long as the spectrum")
                rng = np.random.default_rng(seed)
                starts = rng.integers(0, noise.size - y.size + 1, n_noise)
                blocks = np.array([noise[i:i + y.size] for i in starts])
            blocks = blocks - blocks.mean(axis=1, keepdims=True)
            draws = self._best_many(fitted + blocks, shifts, nonnegative)[:, :nc]
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
