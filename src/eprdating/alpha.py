"""Energy-dependent alpha efficiency.

The luminescence/ESR signal induced by an alpha particle is, to a good
approximation, proportional to the length of its track rather than to the
energy it deposits (Zimmerman 1971; Aitken & Bowman 1975; see Adamiec & Aitken
1998, p. 38). An alpha of initial energy E deposits E over its range R(E), so
its efficiency per unit dose scales as R(E)/E. Relative to a reference energy
(ROSY uses 5.3 MeV)::

    k(E) = k_ref · [R(E)/E] / [R(E_ref)/E_ref]

R(E) is the CSDA range in the medium from the Bethe stopping power for alpha
particles (below 1 MeV the range is extrapolated with R ∝ E^1.5). Alpha
energies per decay are those of Adamiec & Aitken (1998), Tables 2 and 3.
"""

from __future__ import annotations

import math
from functools import cache

from .onegroup import ELEMENTS, HYDROXYAPATITE, Material

M_ALPHA = 3727.379  # MeV
MC2_E = 0.51099895  # MeV
K_BETHE = 0.307075  # MeV cm²/mol

#: (alpha energy of the emitter in MeV, alpha energy released per decay of the
#: chain parent in MeV) per U-series segment, Adamiec & Aitken (1998).
ALPHA_EMITTERS = {
    "U238": [(4.19, 4.19)],
    "U234": [(4.68, 4.68)],
    "Th230": [(4.58, 4.58), (4.77, 4.77)],  # 230Th, 226Ra
    "Rn222": [(5.49, 5.49), (6.00, 6.00), (7.68, 7.68), (5.31, 5.31)],  # 222Rn, 218Po, 214Po, 210Po
    "U235": [(4.27, 4.27)],
    "Pa231": [(4.84, 4.84), (4.95, 0.07), (5.70 / 0.986, 5.70), (5.67, 5.67),  # 231Pa, 227Ac, 227Th, 223Ra
              (6.63, 6.63), (7.39, 7.39), (6.54, 6.54)],  # 219Rn, 215Po, 211Bi
}

#: weights of the U-238 and U-235 chains in the alpha dose of natural U per MeV
#: released (abundance × dose per ppm / energy per decay; A&A 1998 Table 5).
_CHAIN_WEIGHT = {"U238": 0.9929 * 2.685 / 42.7, "U235": 0.0071 * 16.6 / 41.1}
_U235_SEGMENTS = ("U235", "Pa231")


def stopping_power(E: float, material: Material = HYDROXYAPATITE) -> float:
    """Bethe electronic stopping power for alpha particles, MeV cm²/g."""
    g = 1 + E / M_ALPHA
    b2 = 1 - 1 / g**2
    za = sum(w * ELEMENTS[el][0] / ELEMENTS[el][1] for el, w in material.fractions.items())
    lnI = sum(w * ELEMENTS[el][0] / ELEMENTS[el][1] * math.log(ELEMENTS[el][2])
              for el, w in material.fractions.items()) / za
    i_mev = math.exp(lnI) * 1e-6
    return K_BETHE * 4 * za / b2 * (math.log(2 * MC2_E * b2 * g * g / i_mev) - b2)


@cache
def _range_cached(E: float, key: tuple, n: int = 2000) -> float:
    material = _MATERIALS[key]
    e0 = 1.0
    r0 = e0 / (1.5 * stopping_power(e0, material))  # ∫0^e0 dE/S with S ∝ E^-1/2
    if E <= e0:
        return r0 * (E / e0) ** 1.5
    h = (E - e0) / n
    r = sum(h / stopping_power(e0 + (k + 0.5) * h, material) for k in range(n))
    return r + r0


_MATERIALS = {HYDROXYAPATITE.key: HYDROXYAPATITE}  # keyed by composition


def alpha_range(E: float, material: Material = HYDROXYAPATITE) -> float:
    """CSDA range of an alpha particle of energy E (MeV), in g/cm²."""
    _MATERIALS.setdefault(material.key, material)
    return _range_cached(round(E, 6), material.key)


def k_ratio(E: float, e_ref: float = 5.3, material: Material = HYDROXYAPATITE) -> float:
    """Alpha efficiency at energy E relative to that at ``e_ref``."""
    return (alpha_range(E, material) / E) / (alpha_range(e_ref, material) / e_ref)


@cache
def _segment_ratios(e_ref: float, key: tuple) -> tuple:
    material = _MATERIALS[key]
    out = {}
    for seg, emitters in ALPHA_EMITTERS.items():
        tot = sum(w for _, w in emitters)
        out[seg] = sum(w * k_ratio(E, e_ref, material) for E, w in emitters) / tot
    return tuple(out.items())


def segment_k_ratios(e_ref: float = 5.3, material: Material = HYDROXYAPATITE) -> dict[str, float]:
    """Dose-weighted alpha-efficiency ratio k/k_ref for each U-series segment."""
    _MATERIALS.setdefault(material.key, material)
    return dict(_segment_ratios(e_ref, material.key))


def natural_u_k_ratio(e_ref: float = 5.3, material: Material = HYDROXYAPATITE) -> float:
    """k/k_ref for natural U in secular equilibrium."""
    seg = segment_k_ratios(e_ref, material)
    num = den = 0.0
    for s, emitters in ALPHA_EMITTERS.items():
        w = _CHAIN_WEIGHT["U235" if s in _U235_SEGMENTS else "U238"] * sum(x for _, x in emitters)
        num += w * seg[s]
        den += w
    return num / den
