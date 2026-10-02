"""Optional EPRAYA backend: component line shapes from spin-Hamiltonian
simulations of powder samples.

Install with ``pip install eprdating[spectra]``. EPRAYA is imported lazily, so
the rest of the library works without it (and without JAX).

Only two EPRAYA entry points are used (``Start`` and ``Powder``), which keeps
the coupling between both projects small.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .deconvolution import ComponentBasis


def _import_epraya():
    try:
        import epraya
    except ModuleNotFoundError as exc:
        if exc.name == "epraya":
            raise ImportError(
                "EPRAYA is not installed. Install the optional extra: pip install 'eprdating[spectra]'"
            ) from exc
        if exc.name and exc.name.startswith("tkinter"):
            raise ImportError(
                "EPRAYA imports tkinter at import time; install your Python's Tk support "
                "(e.g. `apt install python3-tk`) to use the EPRAYA backend."
            ) from exc
        raise
    return epraya


@dataclass
class Species:
    """Spin-Hamiltonian description of one paramagnetic species (S=1/2 default).

    g    : isotropic float or [gx, gy, gz] principal values.
    Hpp  : [Gaussian, Lorentzian] peak-to-peak linewidths in mT, as in EPRAYA.
    eta  : Gaussian weight of the Voigt profile (EPRAYA convention).
    extra: any other ``Ham`` attribute understood by EPRAYA (A, I, Nucl, ...).
    """

    g: Sequence[float]
    Hpp: Sequence[float] = (0.0, 0.3)
    S: float = 0.5
    eta: float = 0.5
    extra: dict[str, object] = field(default_factory=dict)


def simulate_species(
    species: Species,
    freq_GHz: float,
    field_range_mT: tuple[float, float],
    points: int = 2048,
    grid: int = 70,
) -> tuple[np.ndarray, np.ndarray]:
    """Simulate the first-derivative powder spectrum of one species.

    Returns ``(B_mT, intensity)``. ``grid`` is EPRAYA's Delaunay ``M``.
    """
    epr = _import_epraya()
    Ham, Exp, _ = epr.Start()
    Ham.S = species.S
    Ham.g = list(np.atleast_1d(species.g).astype(float)) if np.ndim(species.g) else float(species.g)
    Ham.Hpp = list(species.Hpp)
    Ham.eta = species.eta
    for k, v in species.extra.items():
        setattr(Ham, k, v)
    Exp.Freq = float(freq_GHz)
    Exp.Points = int(points)
    Exp.Frange = [float(field_range_mT[0]), float(field_range_mT[1])]
    B, y = epr.Powder(Ham, Exp, M=grid, graph=False)
    return np.asarray(B, float), np.asarray(y, float)


def basis_from_species(
    B,
    species: dict[str, Species],
    freq_GHz: float,
    baseline_order: int = 1,
    points: int | None = None,
    grid: int = 70,
) -> ComponentBasis:
    """Build a :class:`ComponentBasis` on the field grid ``B`` (mT) by
    simulating every species with EPRAYA."""
    B = np.asarray(B, float)
    rng = (float(B.min()), float(B.max()))
    shapes = {
        name: simulate_species(sp, freq_GHz, rng, points or max(len(B), 1024), grid)
        for name, sp in species.items()
    }
    return ComponentBasis(B, shapes, baseline_order=baseline_order)
