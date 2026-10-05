"""One-group (double-P0) beta transport in planar layers.

Theory: O'Brien, Samson, Sanna & McLaughlin (1964) Nucl. Sci. Eng. 18, 90;
coefficients as in Prestwich & Chan (2000) Radiat. Phys. Chem. 59, 221.
Applied to tooth enamel by Brennan et al. (1997) Radiat. Meas. 27, 307 (ROSY).

For each beta emitter, the spectrum is replaced by a single group of
electrons at the mean energy E. The forward/backward fluences obey

    ± dΦ±/dz + (2 μa + μs) Φ± = Y(z) + μs Φ∓

with z the mass depth (g/cm²), Y the decays per gram, μa = S(E)/E (collision
stopping power over energy) and μs the Lewis transport cross section of a
screened Coulomb potential. The dose is D = μa E (Φ+ + Φ−); in an infinite
homogeneous medium D = E Y. Every layer is homogeneous, so the solution in a
layer is a constant plus a rising and a decaying exponential, with
attenuation ν = 2 sqrt(μa (μa + μs)). Fluences are continuous at interfaces.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np

from ._types import ValueLike, as_value

N_A = 6.02214076e23
R_E = 2.8179403262e-13  # cm
MC2 = 0.51099895  # MeV
ALPHA = 1 / 137.035999084
K_BETHE = 0.1535375  # 2π N_A r_e² m c² / 2, MeV cm²/mol (with Z/A in mol/g)

#: Z, atomic mass, mean excitation energy (eV, ICRU 37)
ELEMENTS = {
    "H": (1, 1.008, 19.2),
    "C": (6, 12.011, 78.0),
    "N": (7, 14.007, 82.0),
    "O": (8, 15.999, 95.0),
    "Na": (11, 22.990, 149.0),
    "Mg": (12, 24.305, 156.0),
    "Al": (13, 26.982, 166.0),
    "Si": (14, 28.085, 173.0),
    "P": (15, 30.974, 173.0),
    "K": (19, 39.098, 190.0),
    "Ca": (20, 40.078, 191.0),
    "Fe": (26, 55.845, 286.0),
}

#: Below this energy the one-group equations are not meaningful; such
#: emitters are treated at 25 keV, i.e. absorbed almost locally (as in
#: Brennan et al. 1997 for Auger and conversion electrons).
E_MIN = 0.025

_COEF_CACHE: dict = {}


def _kinematics(E):
    tau = E / MC2
    beta2 = tau * (tau + 2) / (tau + 1) ** 2
    p2 = tau * (tau + 2)  # (p / mc)^2
    return tau, beta2, p2


def lewis_mus(Z, A, E):
    """Lewis transport (scattering) cross section, cm²/g (Prestwich & Chan eq. 12-13)."""
    _, b2, p2 = _kinematics(E)
    psi = 0.319 * ALPHA**2 * Z ** (2 / 3) / p2 * (1.13 + 3.76 * ALPHA**2 * Z**2 / b2)
    return 2 * math.pi * N_A * R_E**2 * Z * (Z + 1) / (A * p2 * b2) * (math.log(1 / psi) - 1)


def bethe_stopping(z_over_a, I_eV, E):
    """Collision stopping power for electrons, MeV cm²/g (no density effect)."""
    tau, b2, _ = _kinematics(E)
    i = I_eV * 1e-6 / MC2
    F = 1 - b2 + (tau**2 / 8 - (2 * tau + 1) * math.log(2)) / (tau + 1) ** 2
    return K_BETHE * z_over_a / b2 * (math.log(tau**2 * (tau + 2) / (2 * i**2)) + F)


@dataclass(frozen=True)
class Material:
    """Mass fractions of elements (normalised on construction)."""

    name: str
    fractions: dict

    def __post_init__(self):
        tot = sum(self.fractions.values())
        object.__setattr__(self, "fractions", {k: v / tot for k, v in self.fractions.items()})

    @property
    def key(self) -> tuple:
        """Identity by composition (two materials with the same name but
        different compositions never share cached coefficients)."""
        return tuple(sorted((el, round(w, 12)) for el, w in self.fractions.items()))

    def with_water(self, water: float) -> Material:
        """Add ``water`` grams of water per gram of dry material."""
        if water <= 0:
            return self
        f = {k: v for k, v in self.fractions.items()}
        f["H"] = f.get("H", 0) + water * 2 * 1.008 / 18.015
        f["O"] = f.get("O", 0) + water * 15.999 / 18.015
        return Material(f"{self.name}+{water!r}water", f)

    def coefficients(self, E: float, scatter_factor: float = 1.0) -> tuple[float, float]:
        """(μa, μs) in cm²/g at energy E (MeV)."""
        key = (self.key, round(E, 9), scatter_factor)
        hit = _COEF_CACHE.get(key)
        if hit is not None:
            return hit
        E = max(E, E_MIN)
        za = sum(w * ELEMENTS[el][0] / ELEMENTS[el][1] for el, w in self.fractions.items())
        lnI = sum(w * ELEMENTS[el][0] / ELEMENTS[el][1] * math.log(ELEMENTS[el][2])
                  for el, w in self.fractions.items()) / za
        mua = bethe_stopping(za, math.exp(lnI), E) / E
        mus = sum(w * lewis_mus(ELEMENTS[el][0], ELEMENTS[el][1], E) for el, w in self.fractions.items())
        _COEF_CACHE[key] = (mua, scatter_factor * mus)
        return _COEF_CACHE[key]


def compound(name: str, formula: dict) -> Material:
    """A material from its chemical formula, e.g. ``compound("calcite", {"Ca": 1, "C": 1, "O": 3})``.

    Elements available: ``ELEMENTS`` (H, C, N, O, Na, Mg, Al, Si, P, K, Ca, Fe).
    """
    unknown = set(formula) - set(ELEMENTS)
    if unknown:
        raise ValueError(f"no data for element(s) {sorted(unknown)}; available: {sorted(ELEMENTS)}")
    return Material(name, {el: n * ELEMENTS[el][1] for el, n in formula.items()})


_compound = compound  # backwards compatibility

HYDROXYAPATITE = compound("hydroxyapatite", {"Ca": 10, "P": 6, "O": 26, "H": 2})
SILICA = compound("silica", {"Si": 1, "O": 2})
WATER = compound("water", {"H": 2, "O": 1})
COLLAGEN = Material("collagen", {"C": 0.50, "H": 0.07, "N": 0.17, "O": 0.26})
#: Common sediment minerals (ideal formulas).
CALCITE = compound("calcite", {"Ca": 1, "C": 1, "O": 3})
DOLOMITE = compound("dolomite", {"Ca": 1, "Mg": 1, "C": 2, "O": 6})
KAOLINITE = compound("kaolinite", {"Al": 2, "Si": 2, "O": 9, "H": 4})
ILLITE = compound("illite", {"K": 0.65, "Al": 2.65, "Si": 3.35, "O": 12, "H": 2})
ORTHOCLASE = compound("orthoclase", {"K": 1, "Al": 1, "Si": 3, "O": 8})
HEMATITE = compound("hematite", {"Fe": 2, "O": 3})


def mixture(name: str, parts: list[tuple[Material, float]]) -> Material:
    """Mix materials by mass fraction: ``[(material, mass_fraction), ...]``
    (fractions are normalised)."""
    if any(w < 0 for _, w in parts) or not sum(w for _, w in parts) > 0:
        raise ValueError("mass fractions must be non-negative and not all zero")
    f: dict = {}
    for mat, w in parts:
        for el, x in mat.fractions.items():
            f[el] = f.get(el, 0.0) + w * x
    return Material(name, f)


def sediment_material(quartz: float = 1.0, calcite: float = 0.0, dolomite: float = 0.0, kaolinite: float = 0.0,
                      illite: float = 0.0, feldspar: float = 0.0, iron_oxide: float = 0.0,
                      name: str = "sediment") -> Material:
    """Dry sediment from its mineral mass fractions (normalised), e.g.
    ``sediment_material(quartz=0.6, calcite=0.3, kaolinite=0.1)``; water is
    added separately (``ToothLayers.sediment_water``)."""
    parts = [(SILICA, quartz), (CALCITE, calcite), (DOLOMITE, dolomite), (KAOLINITE, kaolinite), (ILLITE, illite),
             (ORTHOCLASE, feldspar), (HEMATITE, iron_oxide)]
    return mixture(name, [(m, w) for m, w in parts if w > 0])


def dentine_material(mineral: float = 0.70, collagen: float = 0.20, water: float = 0.10,
                     name: str = "dentine") -> Material:
    """Dentine (or cementum) as hydroxyapatite + collagen + water by mass."""
    return mixture(name, [(HYDROXYAPATITE, mineral), (COLLAGEN, collagen), (WATER, water)])


#: Default dentine: 70 % mineral, 20 % collagen, 10 % water by mass (indicative).
DENTINE = dentine_material()


@dataclass
class Layer:
    material: Material
    thickness: float  # g/cm²; math.inf for a half-space
    source: float = 0.0  # decays per gram (relative units)


def solve_fluence(layers: list[Layer], E: float, scatter_factor: float = 1.0):
    """Solve the double-P0 equations for one emitter of energy E.

    ``layers`` runs from left to right; the first and last must be half-spaces
    (``thickness = inf``). Returns a list of per-layer solutions
    ``(mua, mus, nu, P, A, B, t)`` with, in local depth ζ ∈ [0, t]::

        Φ+ = P + A e^{ν(ζ - t)} + B e^{-ν ζ}
        Φ- = P + A g_A e^{ν(ζ - t)} + B g_B e^{-ν ζ}
    """
    n = len(layers)
    if n < 2 or not (math.isinf(layers[0].thickness) and math.isinf(layers[-1].thickness)):
        raise ValueError("first and last layers must be half-spaces (thickness=inf)")
    par = []
    for L in layers:
        mua, mus = L.material.coefficients(E, scatter_factor)
        nu = 2 * math.sqrt(mua * (mua + mus))
        gA = (nu + 2 * mua + mus) / mus
        gB = (2 * mua + mus - nu) / mus
        P = L.source / (2 * mua)
        par.append((mua, mus, nu, gA, gB, P))
    # unknowns: left half-space A0 (rising toward its right edge), finite layers A,B, right half-space B
    idx = {}
    k = 0
    for i, L in enumerate(layers):
        if i == 0:
            idx[(i, "A")] = k; k += 1
        elif i == n - 1:
            idx[(i, "B")] = k; k += 1
        else:
            idx[(i, "A")] = k; idx[(i, "B")] = k + 1; k += 2
    M = np.zeros((k, k))
    rhs = np.zeros(k)
    row = 0
    for i in range(n - 1):
        # interface between layer i (at its right edge) and i+1 (at its left edge)
        for comp in (0, 1):  # 0: Φ+, 1: Φ-
            _, _, nu_l, gA_l, gB_l, P_l = par[i]
            _, _, nu_r, gA_r, gB_r, P_r = par[i + 1]
            t_l = layers[i].thickness
            t_r = layers[i + 1].thickness
            # left side value at ζ = t_l: A e^0 + B e^{-ν t}
            if (i, "A") in idx:
                M[row, idx[(i, "A")]] += 1.0 if comp == 0 else gA_l
            if (i, "B") in idx:
                M[row, idx[(i, "B")]] += math.exp(-nu_l * t_l) * (1.0 if comp == 0 else gB_l)
            # right side value at ζ = 0: A e^{-ν t} + B
            if (i + 1, "A") in idx:
                M[row, idx[(i + 1, "A")]] -= math.exp(-nu_r * t_r) * (1.0 if comp == 0 else gA_r)
            if (i + 1, "B") in idx:
                M[row, idx[(i + 1, "B")]] -= 1.0 if comp == 0 else gB_r
            rhs[row] = P_r - P_l
            row += 1
    x = np.linalg.solve(M, rhs)
    out = []
    for i, L in enumerate(layers):
        mua, mus, nu, gA, gB, P = par[i]
        A = x[idx[(i, "A")]] if (i, "A") in idx else 0.0
        B = x[idx[(i, "B")]] if (i, "B") in idx else 0.0
        out.append((mua, mus, nu, gA, gB, P, A, B, L.thickness))
    return out


def mean_dose(sol, i: int, z0: float, z1: float, E: float) -> float:
    """Mean dose (E × decays/g units) over local depths [z0, z1] of layer i."""
    mua, _, nu, gA, gB, P, A, B, t = sol[i]
    w = z1 - z0
    if w <= 0:
        raise ValueError("empty averaging interval")
    ia = (math.exp(nu * (z1 - t)) - math.exp(nu * (z0 - t))) / nu if A else 0.0
    ib = (math.exp(-nu * z0) - math.exp(-nu * z1)) / nu if B else 0.0
    phi = 2 * P * w + A * (1 + gA) * ia + B * (1 + gB) * ib
    return mua * E * phi / w


@dataclass
class Emitter:
    name: str
    energy: float  # mean energy of the emitted electrons, MeV
    weight: float  # beta energy released per decay of the chain parent, MeV


# Beta emitters with energy per disintegration from Adamiec & Aitken (1998),
# Tables 1-3 and 6. Group energy = energy per decay / beta branching.
U238_EMITTERS = {
    "U238": [Emitter("U-238", 0.007, 0.007), Emitter("Th-234", 0.060, 0.060),
             Emitter("Pa-234m", 0.818, 0.818), Emitter("Pa-234", 0.001 / 0.0016, 0.001)],
    "U234": [Emitter("U-234", 0.012, 0.012)],
    "Th230": [Emitter("Th-230", 0.013, 0.013), Emitter("Ra-226", 0.0038, 0.0038)],
    "Rn222": [Emitter("Pb-214", 0.294, 0.294), Emitter("Bi-214", 0.652, 0.652),
              Emitter("Pb-210", 0.033, 0.033), Emitter("Bi-210", 0.389, 0.389)],
}
U235_EMITTERS = {
    "U235": [Emitter("U-235", 0.037, 0.037), Emitter("Th-231", 0.1506, 0.1506)],
    "Pa231": [Emitter("Pa-231", 0.032, 0.032), Emitter("Th-227", 0.028 / 0.986, 0.028),
              Emitter("Fr-223", 0.006 / 0.014, 0.006), Emitter("Ra-223", 0.066, 0.066),
              Emitter("Rn-219", 0.007, 0.007), Emitter("Pb-211", 0.455, 0.455),
              Emitter("Bi-211", 0.001, 0.001), Emitter("Tl-207", 0.494 / 0.997, 0.494)],
}
TH232_EMITTERS = [Emitter("Ra-228", 0.01, 0.01), Emitter("Ac-228", 0.413, 0.413),
                  Emitter("Th-228", 0.019, 0.019), Emitter("Ra-224", 0.002, 0.002),
                  Emitter("Pb-212", 0.173, 0.173), Emitter("Bi-212", 0.502 / 0.641, 0.502),
                  Emitter("Tl-208", 0.209 / 0.359, 0.209)]
K40_EMITTERS = [Emitter("K-40", 0.501 / 0.893, 0.501)]


def weighted_fraction(emitters, frac_fn) -> float:
    """Energy-weighted mean of ``frac_fn(E)`` over a list of emitters."""
    tot = sum(e.weight for e in emitters)
    return sum(e.weight * frac_fn(e.energy) for e in emitters) / tot


@dataclass
class ToothLayers:
    """Planar tooth geometry for one-group beta attenuation.

    Thicknesses in µm, densities in g/cm³, water in g per g of dry material.
    Order: outer sediment | cementum | enamel | dentine | inner sediment.
    Geometric inputs accept a number, a ``(value, sigma)`` tuple or a
    :class:`~eprdating.Value`; the uncertainties are sampled by
    :meth:`eprdating.ToothSample.age_mc`.
    """

    enamel_um: ValueLike
    dentine_um: ValueLike = 2000.0
    cementum_um: ValueLike = 0.0
    strip_outer_um: ValueLike = 0.0
    strip_inner_um: ValueLike = 0.0
    enamel_density: ValueLike = 3.0
    dentine_density: ValueLike = 2.82
    cementum_density: ValueLike = 2.54
    sediment_density: ValueLike = 2.0
    sediment_water: float = 0.0
    dentine_water: float = 0.0
    cementum_water: float = 0.0
    enamel: Material = HYDROXYAPATITE
    dentine: Material = DENTINE
    cementum: Material = DENTINE
    sediment: Material = SILICA
    scatter_factor: float = 1.0
    _cache: dict = field(default_factory=dict, repr=False)

    #: geometric inputs that may carry an uncertainty
    GEOMETRY = (
        "enamel_um", "dentine_um", "cementum_um", "strip_outer_um", "strip_inner_um",
        "enamel_density", "dentine_density", "cementum_density", "sediment_density",
    )

    def _x(self, name: str) -> float:
        return as_value(getattr(self, name)).value

    def nominal_values(self) -> dict:
        return {k: self._x(k) for k in self.GEOMETRY}

    def at(self, **values) -> ToothLayers:
        """Copy with the given fields set (floats), e.g. one Monte Carlo draw."""
        geo = replace(self, _cache={}, **values)
        en = geo._x("enamel_um")
        if en <= 0:
            raise ValueError("enamel thickness must be positive")
        if geo._x("strip_outer_um") + geo._x("strip_inner_um") >= en:
            raise ValueError("enamel stripped from both sides exceeds the enamel thickness")
        return geo

    def _stack(self, source: str) -> tuple[list[Layer], int]:
        x = self._x
        sed = self.sediment.with_water(self.sediment_water)
        # "sediment" is the sediment on both sides (the inner one matters when
        # there is little or no dentine); "sediment_outer"/"sediment_inner" split it
        layers = [Layer(sed, math.inf, 1.0 if source in ("sediment", "sediment_outer") else 0.0)]
        if x("cementum_um") > 0:
            layers.append(Layer(self.cementum.with_water(self.cementum_water),
                                x("cementum_um") * 1e-4 * x("cementum_density"),
                                1.0 if source == "cementum" else 0.0))
        layers.append(Layer(self.enamel, x("enamel_um") * 1e-4 * x("enamel_density"),
                            1.0 if source == "enamel" else 0.0))
        target = len(layers) - 1
        layers.append(Layer(self.dentine.with_water(self.dentine_water),
                            x("dentine_um") * 1e-4 * x("dentine_density"),
                            1.0 if source == "dentine" else 0.0))
        layers.append(Layer(sed, math.inf, 1.0 if source in ("sediment", "sediment_inner") else 0.0))
        return layers, target

    def fraction(self, source: str, E: float) -> float:
        """Mean dose in the dated enamel per unit infinite-matrix dose of
        ``source`` ('enamel', 'dentine', 'cementum', 'sediment' — both sides —,
        'sediment_outer', 'sediment_inner') for energy E."""
        key = (source, round(E, 6))
        if key not in self._cache:
            layers, t = self._stack(source)
            sol = solve_fluence(layers, E, self.scatter_factor)
            x = self._x
            z0 = x("strip_outer_um") * 1e-4 * x("enamel_density")
            z1 = (x("enamel_um") - x("strip_inner_um")) * 1e-4 * x("enamel_density")
            self._cache[key] = mean_dose(sol, t, z0, z1, E) / E
        return self._cache[key]

    def chain_fraction(self, source: str, chain: str) -> float:
        """Energy-weighted fraction for a whole decay chain or segment.

        ``chain``: 'U' (U-238 + U-235 in equilibrium), 'Th', 'K', or a
        U-series segment name ('U238', 'U234', 'Th230', 'Rn222', 'U235', 'Pa231').
        """
        def f(E):
            return self.fraction(source, E)
        if chain == "Th":
            return weighted_fraction(TH232_EMITTERS, f)
        if chain == "K":
            return weighted_fraction(K40_EMITTERS, f)
        if chain in U238_EMITTERS:
            return weighted_fraction(U238_EMITTERS[chain], f)
        if chain in U235_EMITTERS:
            return weighted_fraction(U235_EMITTERS[chain], f)
        if chain == "U":
            # natural U: weight segments by their beta dose (abundance 99.29 % / 0.71 %)
            parts = []
            for segs, w in ((U238_EMITTERS, 0.9929 * 0.143 / 2.28), (U235_EMITTERS, 0.0071 * 0.515 / 1.27)):
                for emitters in segs.values():
                    for e in emitters:
                        parts.append(Emitter(e.name, e.energy, w * e.weight))
            return weighted_fraction(parts, f)
        raise ValueError(f"unknown chain {chain!r}")
