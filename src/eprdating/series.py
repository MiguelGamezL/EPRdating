"""U-series disequilibrium after uranium uptake.

Uranium taken up by a tooth arrives without its daughters: ``230Th`` (and
everything below it) grows in over ~10^5 years, and ``234U`` may be in excess
of ``238U``. The dose rate of freshly incorporated U is therefore much lower
than the secular-equilibrium value given by the conversion factors.

We split the U-238 chain into segments whose activity relative to ``238U``
evolves as (``tau`` = time since that U was incorporated, ``r0`` = initial
234U/238U activity ratio, 230Th initially absent):

* ``U238``  : 238U → 234Th → 234Pa                 ratio 1
* ``U234``  : 234U                                  1 + (r0-1) e^{-λ4 τ}
* ``Th230`` : 230Th and its daughters to 206Pb      Bateman ingrowth
  (226Ra and below are assumed in equilibrium with 230Th; radon loss is
  not yet modelled)
* ``U235``  : 235U chain, taken as constant (231Pa ingrowth ignored)

The fraction of the equilibrium dose rate carried by each segment, for each
radiation type, must come from a published table (e.g. Adamiec & Aitken
1998; Guérin et al. 2011). Those numbers are **not** bundled yet: fill
``data/u_series_partition.json`` (the loader checks that each radiation's
fractions sum to 1) or pass a dict to :class:`USeries`.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path

#: Half-lives in ka (Cheng et al. 2013 for 234U and 230Th).
HALF_LIFE_KA = {"U238": 4.4683e6, "U234": 245.620, "Th230": 75.584}
LAMBDA = {k: math.log(2) / v for k, v in HALF_LIFE_KA.items()}

SEGMENTS = ("U238", "U234", "Th230", "U235")


def load_partition(path: str | None = None) -> dict[str, dict[str, float]]:
    """Load segment fractions ``{radiation: {segment: fraction}}``."""
    src = resources.files("eprdating.data").joinpath("u_series_partition.json") if path is None else Path(path)
    with src.open(encoding="utf-8") as fh:
        raw = json.load(fh)
    part = {rad: raw[rad] for rad in ("alpha", "beta", "gamma")}
    _check_partition(part)
    return part


def _check_partition(part: dict[str, dict[str, float]]) -> None:
    for rad, seg in part.items():
        if any(v is None for v in seg.values()):
            raise ValueError(
                "U-series partition table is not filled in yet "
                "(eprdating/data/u_series_partition.json). Add the per-segment "
                "fractions from Adamiec & Aitken (1998) or Guérin et al. (2011), "
                "or pass `partition=` explicitly."
            )
        unknown = set(seg) - set(SEGMENTS)
        if unknown:
            raise ValueError(f"unknown segments {unknown}; expected {SEGMENTS}")
        total = sum(seg.values())
        if abs(total - 1.0) > 1e-3:
            raise ValueError(f"{rad} fractions sum to {total:.4f}, expected 1")


def _G_U238(tau: float, r0: float) -> float:
    return tau


def _G_U234(tau: float, r0: float) -> float:
    l4 = LAMBDA["U234"]
    return tau + (r0 - 1.0) * (1.0 - math.exp(-l4 * tau)) / l4


def _G_Th230(tau: float, r0: float) -> float:
    l4, l0 = LAMBDA["U234"], LAMBDA["Th230"]
    e4 = (1.0 - math.exp(-l4 * tau)) / l4
    e0 = (1.0 - math.exp(-l0 * tau)) / l0
    return tau - e0 + (r0 - 1.0) * l0 / (l0 - l4) * (e4 - e0)


_G = {"U238": _G_U238, "U234": _G_U234, "Th230": _G_Th230, "U235": _G_U238}


def activity_ratio_Th230_U238(tau: float, r0: float = 1.0) -> float:
    """(230Th/238U) activity ratio after a time ``tau`` (ka) of a closed system
    that started with 234U/238U = ``r0`` and no 230Th."""
    l4, l0 = LAMBDA["U234"], LAMBDA["Th230"]
    return 1.0 - math.exp(-l0 * tau) + (r0 - 1.0) * l0 / (l0 - l4) * (math.exp(-l4 * tau) - math.exp(-l0 * tau))


@dataclass
class USeries:
    """Time-integrated dose of U incorporated in a tissue, per unit of its
    equilibrium dose rate, accounting for daughter ingrowth.

    ``G(radiation)(tau)`` returns ∫_0^tau (D(s)/D_eq) ds for that radiation.
    """

    r0: float = 1.0
    partition: dict[str, dict[str, float]] | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        if self.partition is None:
            self.partition = load_partition()
        else:
            _check_partition(self.partition)

    def G(self, radiation: str):
        frac = self.partition[radiation]

        def g(tau: float) -> float:
            tau = max(float(tau), 0.0)
            return sum(f * _G[seg](tau, self.r0) for seg, f in frac.items())

        return g
