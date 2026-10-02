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
* ``Th230`` : 230Th → 226Ra                         Bateman ingrowth
* ``Rn222`` : 222Rn and its daughters to 210Po      as Th230, times (1 - f_Rn)
* ``U235``  : 235U → 231Th                          ratio 1
* ``Pa231`` : 231Pa and its daughters               1 - e^{-λPa τ}

226Ra is assumed in equilibrium with 230Th (half-life 1.6 ka) and ``f_Rn``
is the fraction of 222Rn that escapes the tissue (radon loss).

The 234U/238U ratio can be given as the **present-day** (measured) value
(default) or as the initial value. For a present-day ratio, every U parcel is
assumed to show today the measured value, so a parcel incorporated a time
``tau`` ago started with ``r0 = 1 + (r_now - 1) e^{λ4 τ}``. This is exact for
early uptake and a consistent approximation for continuous uptake. For old
samples with a ratio far from 1 the implied ``r0`` grows quickly; check it
with :meth:`USeries.initial_ratio`.

The fraction of the equilibrium dose rate of natural U carried by each
segment, per radiation type, is bundled in ``data/u_series_partition.json``.
It is derived from the energies per disintegration of Adamiec & Aitken
(1998, Tables 2, 3 and 5) by ``tools/derive_u_series_partition.py``. A
different table can be passed to :class:`USeries` (fractions must sum to 1).

For freshly incorporated U only ~20 % (alpha), ~39 % (beta) and ~2 %
(gamma) of the equilibrium dose rate is delivered.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path

#: Half-lives in ka (Cheng et al. 2013 for 234U and 230Th; 32.76 ka for 231Pa).
HALF_LIFE_KA = {"U238": 4.4683e6, "U234": 245.620, "Th230": 75.584, "Pa231": 32.76}
LAMBDA = {k: math.log(2) / v for k, v in HALF_LIFE_KA.items()}

SEGMENTS = ("U238", "U234", "Th230", "Rn222", "U235", "Pa231")


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
            raise ValueError(f"U-series partition table has empty entries for {rad}")
        unknown = set(seg) - set(SEGMENTS)
        if unknown:
            raise ValueError(f"unknown segments {unknown}; expected {SEGMENTS}")
        total = sum(seg.values())
        if abs(total - 1.0) > 1e-3:
            raise ValueError(f"{rad} fractions sum to {total:.4f}, expected 1")


def _G_unit(tau: float, k: float) -> float:
    return tau


def _G_U234(tau: float, k: float) -> float:
    l4 = LAMBDA["U234"]
    return tau + k * (1.0 - math.exp(-l4 * tau)) / l4


def _G_Th230(tau: float, k: float) -> float:
    l4, l0 = LAMBDA["U234"], LAMBDA["Th230"]
    e4 = (1.0 - math.exp(-l4 * tau)) / l4
    e0 = (1.0 - math.exp(-l0 * tau)) / l0
    return tau - e0 + k * l0 / (l0 - l4) * (e4 - e0)


def _G_Pa231(tau: float, k: float) -> float:
    lp = LAMBDA["Pa231"]
    return tau - (1.0 - math.exp(-lp * tau)) / lp


# k = r0 - 1, the initial 234U excess of the parcel
_G = {
    "U238": _G_unit,
    "U234": _G_U234,
    "Th230": _G_Th230,
    "Rn222": _G_Th230,
    "U235": _G_unit,
    "Pa231": _G_Pa231,
}


def activity_ratio_Th230_U238(tau: float, r0: float = 1.0) -> float:
    """(230Th/238U) activity ratio after a time ``tau`` (ka) of a closed system
    that started with 234U/238U = ``r0`` and no 230Th."""
    l4, l0 = LAMBDA["U234"], LAMBDA["Th230"]
    return 1.0 - math.exp(-l0 * tau) + (r0 - 1.0) * l0 / (l0 - l4) * (math.exp(-l4 * tau) - math.exp(-l0 * tau))


@dataclass
class USeries:
    """Time-integrated dose of U incorporated in a tissue, per unit of its
    equilibrium dose rate, accounting for daughter ingrowth.

    Parameters
    ----------
    ratio : 234U/238U activity ratio.
    ratio_is : ``"present"`` (measured today, default) or ``"initial"``.
    radon_loss : fraction of 222Rn escaping the tissue (0 to 1).
    partition : segment fractions; default is the bundled table.

    ``G(radiation)(tau)`` returns ∫_0^tau (D(s)/D_eq) ds for U that has been
    in the tissue for ``tau`` ka.
    """

    ratio: float = 1.0
    ratio_is: str = "present"
    radon_loss: float = 0.0
    partition: dict[str, dict[str, float]] | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        if self.ratio_is not in ("present", "initial"):
            raise ValueError("ratio_is must be 'present' or 'initial'")
        if not 0.0 <= self.radon_loss <= 1.0:
            raise ValueError("radon_loss must be between 0 and 1")
        if self.ratio < 0:
            raise ValueError("234U/238U ratio must be >= 0")
        if self.partition is None:
            self.partition = load_partition()
        else:
            _check_partition(self.partition)
        if self.radon_loss > 0 and not all("Rn222" in seg for seg in self.partition.values()):
            raise ValueError("radon_loss needs a partition table with a separate 'Rn222' segment")

    def initial_ratio(self, tau: float) -> float:
        """Initial 234U/238U of a U parcel incorporated ``tau`` ka ago."""
        return 1.0 + self._k(tau)

    def _k(self, tau: float) -> float:
        if self.ratio_is == "initial":
            return self.ratio - 1.0
        return (self.ratio - 1.0) * math.exp(LAMBDA["U234"] * tau)

    def G(self, radiation: str):
        frac = self.partition[radiation]
        keep_rn = 1.0 - self.radon_loss

        def g(tau: float) -> float:
            tau = max(float(tau), 0.0)
            k = self._k(tau)
            total = 0.0
            for seg, f in frac.items():
                w = f * keep_rn if seg == "Rn222" else f
                total += w * _G[seg](tau, k)
            return total

        return g
