"""Beta-dose geometry for layered (planar) tooth samples.

ROSY (Brennan et al. 1997, 1999) computes the beta dose deposited in the
enamel layer by each neighbouring medium (sediment or cement | enamel |
dentine) with "one-group" beta-transport theory, which accounts for the
finite thickness of each layer and for the removed outer enamel.

**Status (v0.1):** the one-group solver is not implemented. Supply the
attenuation/geometry factors yourself (fraction of the infinite-matrix beta
dose rate of each source that the dated enamel receives), e.g. from the
tables in Brennan et al. (1997) or from Monte Carlo (DosiVox, Geant4)
simulations. :class:`BetaGeometry` documents the expected inputs so the
solver can be dropped in later without changing the age API.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._types import ValueLike, as_value


@dataclass
class BetaGeometry:
    """Fractions of infinite-matrix beta dose rate reaching the dated enamel.

    internal : enamel self-dose (U inside the enamel layer).
    dentine  : from U in the adjacent dentine.
    external : from the sediment (or cement) on the outer side.
    """

    internal: ValueLike
    dentine: ValueLike
    external: ValueLike

    def values(self) -> dict:
        return {k: as_value(getattr(self, k)) for k in ("internal", "dentine", "external")}


@dataclass
class Layers:
    """Thicknesses (µm) and densities (g/cm³) of a planar tooth geometry.

    Intended input for the future one-group solver. The default densities
    are indicative placeholders; set the values you use explicitly.
    """

    enamel_thickness_um: float
    enamel_removed_outer_um: float = 0.0
    enamel_removed_inner_um: float = 0.0
    dentine_thickness_um: float = 2000.0
    cement_thickness_um: float = 0.0
    enamel_density: float = 2.95
    dentine_density: float = 2.85


def one_group_attenuation(layers: Layers) -> BetaGeometry:
    """Beta geometry factors from one-group theory (Brennan et al. 1997).

    Not implemented yet; planned for v0.3.
    """
    raise NotImplementedError(
        "One-group beta attenuation is planned for v0.3. Pass a BetaGeometry "
        "with factors from published tables or Monte Carlo simulations."
    )
