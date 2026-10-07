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


# --------------------------------------------------------------------------
# Alpha escape at layer surfaces
# --------------------------------------------------------------------------
#
# An alpha particle born within its range R of a surface may leave the layer,
# and alphas from the neighbouring medium may enter it. With straight tracks
# and an efficiency proportional to track length (see above), a point at
# depth x (u = x/R <= 1) loses the fraction
#
#     e(u) = [(1 - u) + u ln u] / 2
#
# of the track length of its own alphas through that surface (isotropic
# emission), and receives the same fraction e(u) of the infinite-matrix alpha
# dose of the medium on the other side, scaled by the ratio of the alpha
# ranges per unit mass in that medium and in the layer (an alpha crossing
# a medium of shorter mass range has less of its track left for the layer:
# 0.89 for dentine, 0.96 for silica against hydroxyapatite). Averaged over a
# whole layer of thickness T >= R the loss is R/(8T) per surface. Each
# emitter has its own range.

#: alpha emitters of the 232Th chain: (energy in MeV, energy released per
#: decay of the chain parent in MeV), Adamiec & Aitken (1998); 212Bi decays
#: by alpha in 36 % of cases, 212Po follows the other 64 %.
TH232_EMITTERS = [(4.01, 4.01), (5.42, 5.42), (5.69, 5.69), (6.29, 6.29), (6.78, 6.78),
                  (6.05, 0.36 * 6.05), (8.78, 0.64 * 8.78)]


def _escape_integral(u: float) -> float:
    """∫_0^u e(t) dt, with e(t) = [(1 - t) + t ln t]/2 for t <= 1 and 0 beyond."""
    u = min(max(u, 0.0), 1.0)
    if u == 0.0:
        return 0.0
    return 0.5 * (u - 0.75 * u * u + 0.5 * u * u * math.log(u))


def surface_fraction(a_um: float, b_um: float, range_um: float) -> float:
    """Mean of e(x/R) over distances x from a surface between ``a_um`` and
    ``b_um``: the fraction of the alpha track length lost through that
    surface (or received from beyond it), averaged over the slice."""
    if b_um <= a_um:
        x = max(a_um, 0.0)
        u = x / range_um
        return 0.0 if u >= 1 else 0.5 * ((1 - u) + (u * math.log(u) if u > 0 else 0.0))
    return range_um * (_escape_integral(b_um / range_um) - _escape_integral(a_um / range_um)) / (b_um - a_um)


def slab_fractions(thickness_um: float, strip_outer_um: float, strip_inner_um: float,
                   range_um: float) -> tuple[float, float, float]:
    """``(own, from_outer, from_inner)`` for the measured part of a layer.

    own        : fraction of the layer's own infinite-matrix alpha dose kept;
    from_outer : fraction of the outer medium's infinite-matrix alpha dose
                 received (likewise ``from_inner``).
    The measured part is what remains after stripping ``strip_outer_um`` and
    ``strip_inner_um`` from the two faces.
    """
    a, b = strip_outer_um, thickness_um - strip_inner_um
    if b < a:
        raise ValueError("the stripping removes the whole layer")
    outer = surface_fraction(a, b, range_um)
    inner = surface_fraction(strip_inner_um, thickness_um - strip_outer_um, range_um)
    return 1.0 - outer - inner, outer, inner


def _weighted(emitters, fn, energy: bool, e_ref: float, material: Material) -> float:
    w = [x * (k_ratio(E, e_ref, material) if energy else 1.0) for E, x in emitters]
    return sum(wi * fn(E) for wi, (E, _) in zip(w, emitters, strict=True)) / sum(w)


def escape_fractions(thickness_um: float, strip_outer_um: float = 0.0, strip_inner_um: float = 0.0,
                     density: float = 3.0, energy: bool = False, e_ref: float = 5.3,
                     material: Material = HYDROXYAPATITE, outer_material: Material | None = None,
                     inner_material: Material | None = None) -> dict:
    """Alpha escape for a layer, per U-series segment and for the 232Th chain.

    Returns ``{"own": {segment: f}, "outer": {...}, "inner": {...}}`` with the
    U-series segments of :data:`ALPHA_EMITTERS` and ``"Th232"``; each value is
    averaged over the emitters of the segment weighted by their alpha dose
    (times the relative efficiency k(E)/k_ref when ``energy``). Ranges are
    those in ``material`` at ``density`` (g/cm³); the incoming fractions are
    scaled by the mass-range ratio of ``outer_material`` / ``inner_material``
    (the neighbouring media) to ``material`` when they are given.
    """
    out: dict = {"own": {}, "outer": {}, "inner": {}}
    chains = {**ALPHA_EMITTERS, "Th232": TH232_EMITTERS}
    fractions = {}
    for emitters in chains.values():
        for E, _ in emitters:
            if E not in fractions:
                r_mass = alpha_range(E, material)
                own, outer, inner = slab_fractions(thickness_um, strip_outer_um, strip_inner_um,
                                                   r_mass / density * 1e4)  # range in µm
                if outer_material is not None:
                    outer *= alpha_range(E, outer_material) / r_mass
                if inner_material is not None:
                    inner *= alpha_range(E, inner_material) / r_mass
                fractions[E] = (own, outer, inner)
    for seg, emitters in chains.items():
        for i, key in enumerate(("own", "outer", "inner")):
            out[key][seg] = _weighted(emitters, lambda E, i=i: fractions[E][i], energy, e_ref, material)
    return out


def th232_k_ratio(e_ref: float = 5.3, material: Material = HYDROXYAPATITE) -> float:
    """k/k_ref for the 232Th chain in secular equilibrium."""
    tot = sum(x for _, x in TH232_EMITTERS)
    return sum(x * k_ratio(E, e_ref, material) for E, x in TH232_EMITTERS) / tot
