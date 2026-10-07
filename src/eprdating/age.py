"""Age calculation: solve ∫_0^T Ḋ(t) dt = De for the age T.

Generic layer
-------------
:class:`DoseRateComponent` is one contribution to the dose rate with its
present-day (equilibrium) value, an optional uptake model and an optional
U-series ingrowth function. :func:`solve_age` combines any list of them.

Tooth-enamel layer
------------------
:class:`ToothSample` assembles the usual components of ESR dating of
enamel (internal alpha and beta, beta from dentine and sediment, gamma,
cosmic) and provides a Monte Carlo age with :meth:`ToothSample.age_mc`.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from ._types import ValueLike, as_value
from .alpha import escape_fractions, natural_u_k_ratio, segment_k_ratios, th232_k_ratio
from .beta import BetaGeometry
from .dose_rate import (
    DEFAULT_FACTORS,
    K_ENAMEL,
    RA226_SHARE_OF_TH230,
    ConversionFactors,
    Sediment,
    conversion_factors,
    matrix_dose_rates,
    u_series_split,
    water_correction,
)
from .history import History, integrate
from .onegroup import ToothLayers
from .series import USeries, load_partition
from .uptake import USModel


@dataclass
class DoseRateComponent:
    """One contribution to the total dose rate (Gy/ka).

    rate   : present-day dose rate the source would deliver in secular
             equilibrium with its present U content.
    uptake : None for a constant source, else a :class:`USModel`.
    G      : time-integrated activity ratio of incorporated U
             (``tau -> ∫ D/D_eq``); None means secular equilibrium.
    profile: for a constant source whose rate changed in the past,
             ``(breaks, rates)``: piecewise-constant rates in ka before
             present (see :class:`eprdating.history.History`); ``rate`` is
             then ``rates[0]``.
    """

    name: str
    rate: float
    uptake: USModel | None = None
    G: Callable[[float], float] | None = None
    profile: tuple[Sequence[float], Sequence[float]] | None = None

    def __post_init__(self) -> None:
        if self.profile is not None:
            if self.uptake is not None or self.G is not None:
                raise ValueError("a rate profile is only supported for a constant source")
            breaks, rates = self.profile
            if len(rates) != len(breaks) + 1:
                raise ValueError("a profile needs one rate more than breaks")

    def accumulated(self, T: float) -> float:
        """Dose (Gy) delivered by this component over ``T`` ka."""
        if self.profile is not None:
            return integrate(*self.profile, T)
        if self.rate == 0:
            return 0.0
        if self.uptake is None:
            if self.G is None:
                return self.rate * T
            return self.rate * self.G(T)
        return self.rate * self.uptake.accumulated(T, self.G)


@dataclass
class AgeResult:
    age: float  # ka
    De: float  # Gy
    components: dict[str, float]  # present-day dose rates, Gy/ka
    accumulated: dict[str, float]  # dose delivered by each component, Gy

    @property
    def mean_dose_rate(self) -> float:
        """Time-averaged dose rate De / T (Gy/ka)."""
        return self.De / self.age

    def summary(self) -> str:
        lines = [f"Age = {self.age:.4g} ka   (De = {self.De:.4g} Gy, <Ḋ> = {self.mean_dose_rate:.4g} Gy/ka)"]
        for k, r in self.components.items():
            if r == 0 and self.accumulated[k] == 0:
                continue
            share = 100 * self.accumulated[k] / self.De
            lines.append(f"  {k:16s} Ḋ_now = {r:8.4f} Gy/ka   contributes {share:5.1f} % of De")
        return "\n".join(lines)


def solve_age(De: float, components: Sequence[DoseRateComponent], t_max: float = 1e5) -> AgeResult:
    """Find the age T (ka) such that the accumulated dose equals ``De`` (Gy)."""
    if De <= 0:
        raise ValueError("De must be positive")
    comps = list(components)

    def f(T: float) -> float:
        return sum(c.accumulated(T) for c in comps) - De

    lo, hi = 1e-9, 1.0
    while f(hi) < 0:
        hi *= 4.0
        if hi > t_max:
            raise RuntimeError(f"no solution below {t_max} ka: dose rate too low for this De")
    T = brentq(f, lo, hi, xtol=1e-10, rtol=1e-10)
    return AgeResult(
        age=T,
        De=De,
        components={c.name: c.rate for c in comps},
        accumulated={c.name: c.accumulated(T) for c in comps},
    )


# --------------------------------------------------------------------------
# Tooth enamel
# --------------------------------------------------------------------------


@dataclass
class AgeMC:
    """Monte Carlo age distribution."""

    samples: np.ndarray
    nominal: AgeResult
    n_failed: int = 0

    @property
    def mean(self) -> float:
        return float(np.mean(self.samples))

    @property
    def std(self) -> float:
        return float(np.std(self.samples, ddof=1))

    def interval(self, level: float = 0.68) -> tuple:
        a = (1 - level) / 2
        return tuple(np.quantile(self.samples, [a, 1 - a]))

    def summary(self) -> str:
        lo, hi = self.interval(0.68)
        lo95, hi95 = self.interval(0.95)
        return (
            f"Age (nominal) = {self.nominal.age:.4g} ka\n"
            f"MC: mean = {self.mean:.4g} ± {self.std:.2g} ka "
            f"(68 %: {lo:.4g}–{hi:.4g}; 95 %: {lo95:.4g}–{hi95:.4g}; n = {self.samples.size}"
            + (f", {self.n_failed} failed" if self.n_failed else "")
            + ")"
        )


@dataclass
class ToothSample:
    """ESR dating of tooth enamel with the classical component model.

    Concentrations: U in ppm. Dose in Gy, dose rates in Gy/ka, ages in ka.

    Parameters
    ----------
    De : equivalent dose of the enamel.
    enamel_U, dentine_U : present-day U in each tissue.
    sediment : U, Th, K and water content of the surrounding sediment. The
        water may be a :class:`~eprdating.history.History` (wetter or drier
        periods); the sediment beta and gamma dose rates then follow it,
        and the present-day value (its first segment) is used elsewhere.
    beta : either fixed geometry factors (:class:`eprdating.beta.BetaGeometry`)
        or a layered geometry (:class:`eprdating.onegroup.ToothLayers`), in
        which case beta attenuation is computed with one-group theory for each
        emitter and U-series segment (as ROSY does). Uncertainties on the
        layer thicknesses, stripping and densities are sampled in
        :meth:`age_mc` together with the water contents.
    gamma : external gamma dose rate. If None, computed from ``sediment``
        as an infinite matrix (use in-situ measurements when available).
        A measured value is taken as today's; with a water history it is
        rescaled to the water of each period. A
        :class:`~eprdating.history.History` is used as given.
    cosmic : cosmic dose rate (see :func:`eprdating.dose_rate.cosmic_dose_rate`),
        or a :class:`~eprdating.history.History` of it, e.g. from a burial
        depth history with :func:`eprdating.dose_rate.cosmic_history`.
    k_alpha : alpha efficiency of enamel. With ``alpha_efficiency="energy"``
        it is the value at ``alpha_eref`` MeV.
    alpha_efficiency : ``"constant"`` (default, as in DATA) or ``"energy"``
        (as in ROSY): k varies with alpha energy as R(E)/E, so each U-series
        segment gets its own efficiency (see :mod:`eprdating.alpha`).
    alpha_eref : reference alpha energy for ``k_alpha`` in MeV (ROSY: 5.3).
    alpha_escape : account for alpha particles crossing the enamel surfaces
        (needs a :class:`~eprdating.onegroup.ToothLayers` geometry; default
        False, as DATA). Alphas born within their range (12-40 µm) of a
        surface partly leave the enamel, and those of the dentine, the
        sediment (or cementum) partly enter it; both are averaged over the
        enamel left after stripping, emitter by emitter (see
        :func:`eprdating.alpha.escape_fractions`). The incoming alphas appear
        as the components "dentine alpha" and "sediment alpha" (or
        "cementum alpha"). With 20-40 µm stripped on each side the effect
        vanishes; without stripping it is about R/(8T) per surface (~1 % for
        300 µm of enamel).
    dentine_water : water content used for the dentine beta contribution.
    enamel_water : water content of the enamel (corrects the internal alpha
        and beta dose rates, Zimmerman coefficients).
    cementum_U, cementum_water, uptake_cementum, u234_u238_cementum,
    radon_loss_cementum : the same for a cementum layer on the outer side
        of the enamel (needs a :class:`~eprdating.onegroup.ToothLayers`
        geometry with ``cementum_um > 0``).
    uptake_enamel, uptake_dentine : uptake models (EU, LU, US with p).
    u234_u238_enamel, u234_u238_dentine : 234U/238U activity ratio of each
        tissue, measured today (``u234_u238_is="present"``, default) or of the
        incoming uranium (``"initial"``, used by :mod:`eprdating.usesr`).
    radon_loss_enamel, radon_loss_dentine : fraction of 222Rn escaping each
        tissue (0 to 1).
    ingrowth : model U-series daughter ingrowth after uptake (default). With
        ``False`` the tissues are taken in secular equilibrium (warns).
    partition : optional custom U-series segment table
        (see :mod:`eprdating.series`).
    factors : name of the conversion-factor set.
    sample_geometry : with a :class:`~eprdating.onegroup.ToothLayers`
        geometry, recompute the one-group factors for every Monte Carlo draw
        (default). With ``False`` they are kept at their nominal values.
    beta_by_segment : with a :class:`~eprdating.onegroup.ToothLayers`
        geometry, attenuate the beta dose of each U-series segment with its
        own factor (default, as ROSY). ``False`` applies one factor for the
        whole chain to the dose with ingrowth, as the DATA program (Grün
        2009) does; for young teeth this lowers the dentine beta dose by up
        to ~40 % because the hard 234mPa betas dominate before 226Ra grows in.
    """

    De: ValueLike
    enamel_U: ValueLike
    dentine_U: ValueLike
    sediment: Sediment
    beta: BetaGeometry | ToothLayers
    cosmic: ValueLike | History
    gamma: ValueLike | History | None = None
    k_alpha: ValueLike = K_ENAMEL
    dentine_water: ValueLike = 0.0
    uptake_enamel: USModel = field(default_factory=lambda: USModel(-1.0))
    uptake_dentine: USModel = field(default_factory=lambda: USModel(-1.0))
    u234_u238_enamel: ValueLike = 1.0
    u234_u238_dentine: ValueLike = 1.0
    radon_loss_enamel: ValueLike = 0.0
    radon_loss_dentine: ValueLike = 0.0
    enamel_water: ValueLike = 0.0
    cementum_U: ValueLike = 0.0
    cementum_water: ValueLike = 0.0
    uptake_cementum: USModel = field(default_factory=lambda: USModel(-1.0))
    u234_u238_cementum: ValueLike = 1.0
    radon_loss_cementum: ValueLike = 0.0
    ingrowth: bool = True
    partition: dict | None = None
    factors: str = DEFAULT_FACTORS
    sample_geometry: bool = True
    alpha_efficiency: str = "constant"
    u234_u238_is: str = "present"
    alpha_eref: float = 5.3
    beta_by_segment: bool = True
    alpha_escape: bool = False

    _SCALARS = (
        "De", "enamel_U", "dentine_U", "k_alpha", "dentine_water",
        "u234_u238_enamel", "u234_u238_dentine", "radon_loss_enamel", "radon_loss_dentine",
        "enamel_water", "cementum_U", "cementum_water", "u234_u238_cementum", "radon_loss_cementum",
    )

    # ---- parameter handling -------------------------------------------
    def _environment(self) -> dict:
        """External inputs that may be histories. In the value dicts a history
        ``h`` of key ``k`` gives ``k#i`` for each segment and ``k`` = today."""
        return {"cosmic": self.cosmic, "gamma": self.gamma, "sed_water": self.sediment.water}

    def _histories(self) -> dict[str, History]:
        return {k: x for k, x in self._environment().items() if isinstance(x, History)}

    def _nominal(self) -> dict[str, float]:
        v = {k: as_value(getattr(self, k)).value for k in self._SCALARS}
        for key, x in self._environment().items():
            if isinstance(x, History):
                v.update({f"{key}#{i}": val for i, val in enumerate(x.nominal())})
                v[key] = v[f"{key}#0"]
            elif x is not None:
                v[key] = as_value(x).value
        for k in ("U", "Th", "K"):
            v["sed_" + k] = as_value(getattr(self.sediment, k)).value
        if self.sediment.U_ra226 is not None:
            v["sed_U_ra226"] = as_value(self.sediment.U_ra226).value
        if isinstance(self.beta, BetaGeometry):
            for k, val in self.beta.values().items():
                v["beta_" + k] = val.value
        else:
            for k, val in self.beta.nominal_values().items():
                v["geo_" + k] = val
        return v

    def _sample(self, rng: np.random.Generator, n: int) -> dict[str, np.ndarray]:
        s = {k: as_value(getattr(self, k)).sample(rng, n) for k in self._SCALARS}
        for key, x in self._environment().items():
            if isinstance(x, History):
                s.update({f"{key}#{i}": arr for i, arr in enumerate(x.sample(rng, n))})
                s[key] = s[f"{key}#0"]
            elif x is not None:
                s[key] = as_value(x).sample(rng, n)
        for k in ("U", "Th", "K"):
            s["sed_" + k] = as_value(getattr(self.sediment, k)).sample(rng, n)
        if self.sediment.U_ra226 is not None:
            s["sed_U_ra226"] = as_value(self.sediment.U_ra226).sample(rng, n)
        if isinstance(self.beta, BetaGeometry):
            for k, val in self.beta.values().items():
                s["beta_" + k] = val.sample(rng, n)
        elif self.sample_geometry:
            for k in self.beta.GEOMETRY:
                s["geo_" + k] = as_value(getattr(self.beta, k)).sample(rng, n)
        # physical constraints: no negative concentrations / fractions
        for k, arr in s.items():
            np.clip(arr, 0.0, None, out=arr)
        for k in ("radon_loss_enamel", "radon_loss_dentine", "radon_loss_cementum"):
            np.clip(s[k], 0.0, 1.0, out=s[k])
        return s

    def _useries(self, tissue: str, v: dict[str, float]) -> USeries | None:
        if not self.ingrowth:
            return None
        return USeries(
            ratio=v[f"u234_u238_{tissue}"],
            ratio_is=self.u234_u238_is,
            radon_loss=v[f"radon_loss_{tissue}"],
            partition=self._partition(),
        )

    def _partition(self) -> dict:
        if self.partition is None:
            self.partition = load_partition()
        return self.partition

    # ---- one-group beta factors ------------------------------------------
    def _onegroup(self, v: dict[str, float], seg: int | None = None) -> dict:
        """Per-source, per-chain/segment one-group factors for the geometry and
        water contents in ``v`` (nominal geometry if ``v`` has none);
        ``seg`` selects a segment of a sediment-water history."""
        wkey = "sed_water" if seg is None else f"sed_water#{seg}"
        geo_vals = {k: v.get("geo_" + k, as_value(getattr(self.beta, k)).value) for k in self.beta.GEOMETRY}
        if not self.sample_geometry:
            geo_vals = self.beta.nominal_values()
            nom = self._nominal()
            water = (nom[wkey], nom["dentine_water"], nom["cementum_water"])
        else:
            water = (v[wkey], v["dentine_water"], v["cementum_water"])
        key = (tuple(geo_vals.values()), water)
        cache = self.__dict__.setdefault("_og_cache", {})
        if key not in cache:
            geo = self.beta.at(**geo_vals, sediment_water=water[0], dentine_water=water[1],
                               cementum_water=water[2])
            segs = ("U238", "U234", "Th230", "Rn222", "U235", "Pa231")
            if len(cache) > 8:  # keep the nominal entries, not every MC draw
                cache.clear()
            cache[key] = {
                "enamel": {s: geo.chain_fraction("enamel", s) for s in segs},
                "dentine": {s: geo.chain_fraction("dentine", s) for s in segs},
                "enamel_U": geo.chain_fraction("enamel", "U"),
                "dentine_U": geo.chain_fraction("dentine", "U"),
                "cementum": ({s: geo.chain_fraction("cementum", s) for s in segs}
                             if geo_vals["cementum_um"] > 0 else None),
                "cementum_U": geo.chain_fraction("cementum", "U") if geo_vals["cementum_um"] > 0 else 0.0,
                "sediment": {c: geo.chain_fraction("sediment", c) for c in ("U", "Th", "K")},
                "sediment_seg": {s: geo.chain_fraction("sediment", s) for s in segs},
            }
        return cache[key]

    # ---- model -----------------------------------------------------------
    def components(
        self,
        v: dict[str, float] | None = None,
        uptake_enamel: USModel | None = None,
        uptake_dentine: USModel | None = None,
        uptake_cementum: USModel | None = None,
    ) -> list[DoseRateComponent]:
        """Dose-rate components for the input values ``v`` (nominal if None).

        ``uptake_enamel`` / ``uptake_dentine`` override the sample's uptake
        models (used by the US-ESR solver).
        """
        v = v or self._nominal()
        up_e = uptake_enamel or self.uptake_enamel
        up_d = uptake_dentine or self.uptake_dentine
        up_c = uptake_cementum or self.uptake_cementum
        cf: ConversionFactors = conversion_factors(self.factors)
        cU = {r: cf.get("U", r).value for r in ("alpha", "beta", "gamma")}
        use, usd, usc = (self._useries(t, v) for t in ("enamel", "dentine", "cementum"))
        if self.alpha_efficiency == "constant":
            k_scale, Ga_e = 1.0, (use.G("alpha") if use else None)
        elif self.alpha_efficiency == "energy":
            seg_k = segment_k_ratios(self.alpha_eref)
            k_scale = natural_u_k_ratio(self.alpha_eref)
            Ga_e = use.G("alpha", {s: r / k_scale for s, r in seg_k.items()}) if use else None
        else:
            raise ValueError("alpha_efficiency must be 'constant' or 'energy'")

        if isinstance(self.beta, BetaGeometry):
            Gb_e = use.G("beta") if use else None
            Gb_d = usd.G("beta") if usd else None
            en_beta = v["beta_internal"] * v["enamel_U"] * cU["beta"]
            den_beta = v["beta_dentine"] * water_correction(v["dentine_U"] * cU["beta"], v["dentine_water"], "beta")
            if v["cementum_U"] > 0:
                raise ValueError("cementum U needs a ToothLayers geometry with a cementum layer")
            cem_beta, Gb_c = 0.0, None
        else:
            og = self._onegroup(v)
            # one-group factors are relative to the wet medium's own infinite-matrix
            # dose, i.e. the dry dose rate diluted by (1 + water)
            # rates are the attenuated equilibrium values; the per-segment factors
            # enter G relative to the whole-chain factor
            w_e = {s: f / og["enamel_U"] for s, f in og["enamel"].items()} if self.beta_by_segment else None
            w_d = {s: f / og["dentine_U"] for s, f in og["dentine"].items()} if self.beta_by_segment else None
            Gb_e = use.G("beta", w_e) if use else None
            Gb_d = usd.G("beta", w_d) if usd else None
            en_beta = v["enamel_U"] * cU["beta"] * og["enamel_U"]
            den_beta = v["dentine_U"] * cU["beta"] / (1.0 + v["dentine_water"]) * og["dentine_U"]
            if v["cementum_U"] > 0:
                if og["cementum"] is None:
                    raise ValueError("cementum U needs a ToothLayers geometry with cementum_um > 0")
                w_c = ({s: f / og["cementum_U"] for s, f in og["cementum"].items()}
                       if self.beta_by_segment else None)
                Gb_c = usc.G("beta", w_c) if usc else None
                cem_beta = v["cementum_U"] * cU["beta"] / (1.0 + v["cementum_water"]) * og["cementum_U"]
            else:
                cem_beta, Gb_c = 0.0, None
        alpha_in: list[DoseRateComponent] = []
        if self.alpha_escape:
            k_scale, Ga_e, alpha_in = self._alpha_escape(v, cf, use, usd, usc, up_d, up_c)
        en_alpha = water_correction(k_scale * v["k_alpha"] * v["enamel_U"] * cU["alpha"], v["enamel_water"], "alpha")
        en_beta = water_correction(en_beta, v["enamel_water"], "beta")
        return [
            DoseRateComponent("enamel alpha", en_alpha, up_e, Ga_e),
            *alpha_in,
            DoseRateComponent("enamel beta", en_beta, up_e, Gb_e),
            DoseRateComponent("dentine beta", den_beta, up_d, Gb_d),
            DoseRateComponent("cementum beta", cem_beta, up_c, Gb_c),
            *self._external(v, cf),
        ]

    def _alpha_escape(self, v, cf, use, usd, usc, up_d, up_c):
        """Enamel alpha with escape through its surfaces, and the alpha dose
        entering it from the neighbouring layers: ``(k_scale, G, components)``."""
        if not isinstance(self.beta, ToothLayers):
            raise TypeError("alpha_escape needs a ToothLayers geometry (thicknesses and stripping)")
        geo = (self.beta.nominal_values() if not self.sample_geometry else
               {k: v.get("geo_" + k, as_value(getattr(self.beta, k)).value) for k in self.beta.GEOMETRY})
        energy = self.alpha_efficiency == "energy"
        outer = self.beta.cementum if geo["cementum_um"] > 0 else self.beta.sediment
        key = (geo["enamel_um"], geo["strip_outer_um"], geo["strip_inner_um"], geo["enamel_density"],
               energy, self.alpha_eref, outer.key)
        cache = self.__dict__.setdefault("_alpha_cache", {})
        if key not in cache:
            if len(cache) > 8:
                cache.clear()
            cache[key] = escape_fractions(*key[:4], energy=energy, e_ref=self.alpha_eref,
                                          material=self.beta.enamel, outer_material=outer,
                                          inner_material=self.beta.dentine)
        esc = cache[key]
        seg_k = segment_k_ratios(self.alpha_eref) if energy else dict.fromkeys(esc["own"], 1.0)
        part = self._partition()["alpha"]

        def weights(side):
            w = {s: seg_k[s] * esc[side][s] for s in part}
            norm = sum(part[s] * w[s] for s in part)
            return norm, ({s: x / norm for s, x in w.items()} if norm > 0 else None)

        k_alpha, cU = v["k_alpha"], cf.get("U", "alpha").value
        k_own, w_own = weights("own")
        G_own = use.G("alpha", w_own) if use else None
        comps = []
        k_in, w_in = weights("inner")
        rate = water_correction(k_alpha * k_in * v["dentine_U"] * cU, v["dentine_water"], "alpha")
        comps.append(DoseRateComponent("dentine alpha", rate, up_d, usd.G("alpha", w_in) if usd and w_in else None))
        k_out, w_out = weights("outer")
        if geo["cementum_um"] > 0:
            rate = water_correction(k_alpha * k_out * v["cementum_U"] * cU, v["cementum_water"], "alpha")
            comps.append(DoseRateComponent("cementum alpha", rate, up_c,
                                           usc.G("alpha", w_out) if usc and w_out else None))
        else:
            U, ura = v["sed_U"], v.get("sed_U_ra226")
            mult = dict.fromkeys(part, 1.0)
            if ura is not None and U > 0:
                mult["Rn222"] = ura / U
                mult["Th230"] = 1.0 + (ura / U - 1.0) * RA226_SHARE_OF_TH230["alpha"]
            u_frac = sum(part[s] * seg_k[s] * esc["outer"][s] * mult[s] for s in part)
            th_k = th232_k_ratio(self.alpha_eref) if energy else 1.0
            dry = (U * cU * u_frac
                   + v["sed_Th"] * cf.get("Th", "alpha").value * th_k * esc["outer"]["Th232"])
            comps.append(DoseRateComponent("sediment alpha",
                                           water_correction(k_alpha * dry, v["sed_water"], "alpha")))
        return k_own, G_own, comps

    def _external(self, v: dict[str, float], cf: ConversionFactors) -> list[DoseRateComponent]:
        """Sediment beta, gamma and cosmic components, following any history."""
        hist = self._histories()

        def segments(key):
            return [v[f"{key}#{i}"] for i in range(len(hist[key]))]

        ura = v.get("sed_U_ra226")
        dry = matrix_dose_rates(v["sed_U"], v["sed_Th"], v["sed_K"], cf, U_ra226=ura)
        w_now = v["sed_water"]
        if "sed_water" in hist:
            waters, wbreaks = segments("sed_water"), hist["sed_water"].breaks
        else:
            waters, wbreaks = [w_now], ()
        sed_beta = [self._sediment_beta(v, cf, dry, w, i if "sed_water" in hist else None)
                    for i, w in enumerate(waters)]
        if "gamma" in hist:
            gamma, gbreaks = segments("gamma"), hist["gamma"].breaks
        elif "gamma" in v:  # measured today: rescale to the water of each period
            gamma = [v["gamma"] * water_correction(1.0, w, "gamma") / water_correction(1.0, w_now, "gamma")
                     for w in waters]
            gbreaks = wbreaks
        else:
            gamma, gbreaks = [water_correction(dry["gamma"], w, "gamma") for w in waters], wbreaks
        if "cosmic" in hist:
            cosmic, cbreaks = segments("cosmic"), hist["cosmic"].breaks
        else:
            cosmic, cbreaks = [v["cosmic"]], ()

        def component(name, rates, breaks):
            return DoseRateComponent(name, rates[0], profile=(breaks, rates) if breaks else None)

        return [component("sediment beta", sed_beta, wbreaks), component("gamma", gamma, gbreaks),
                component("cosmic", cosmic, cbreaks)]

    def _sediment_beta(self, v: dict[str, float], cf: ConversionFactors, dry: dict, water: float,
                       seg: int | None) -> float:
        """Beta dose rate from the sediment into the enamel for one water content."""
        if isinstance(self.beta, BetaGeometry):
            return v["beta_external"] * water_correction(dry["beta"], water, "beta")  # diseq. via dry
        og = self._onegroup(v, seg)
        U, ura = v["sed_U"], v.get("sed_U_ra226")
        u_beta = U * og["sediment"]["U"]
        if ura is not None and U > 0:
            # U chain out of equilibrium: rescale with the per-segment attenuation
            pb = self._partition()["beta"]
            sh = RA226_SHARE_OF_TH230["beta"]
            segf = og["sediment_seg"]
            mult = {s: 1.0 for s in segf}
            mult["Rn222"] = ura / U
            mult["Th230"] = 1.0 + (ura / U - 1.0) * sh
            eq = sum(pb[s] * segf[s] for s in segf)
            u_beta *= sum(pb[s] * segf[s] * mult[s] for s in segf) / eq
        elif ura is not None:
            u_beta = ura * og["sediment"]["U"] * u_series_split("beta")[1]
        return (
            u_beta * cf.get("U", "beta").value
            + sum(v["sed_" + nuc] * cf.get(nuc, "beta").value * og["sediment"][nuc] for nuc in ("Th", "K"))
        ) / (1.0 + water)

    def age(self) -> AgeResult:
        """Nominal age with the central value of every input."""
        if not self.ingrowth and (as_value(self.enamel_U).value > 0 or as_value(self.dentine_U).value > 0):
            warnings.warn(
                "Secular equilibrium assumed for U in the dental tissues (ingrowth=False). "
                "This ignores 230Th ingrowth and overestimates the internal dose rate, "
                "especially for samples younger than a few hundred ka.",
                stacklevel=2,
            )
        v = self._nominal()
        return solve_age(v["De"], self.components(v))

    def age_mc(self, n: int = 2000, seed: int | None = None) -> AgeMC:
        """Monte Carlo age: every input is sampled from its Gaussian uncertainty."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nominal = self.age()
        rng = np.random.default_rng(seed)
        s = self._sample(rng, n)
        ages = np.full(n, np.nan)
        for i in range(n):
            v = {k: float(arr[i]) for k, arr in s.items()}
            if v["De"] <= 0:
                continue
            try:
                ages[i] = solve_age(v["De"], self.components(v)).age
            except (RuntimeError, ValueError):
                pass
        ok = np.isfinite(ages)
        return AgeMC(samples=ages[ok], nominal=nominal, n_failed=int((~ok).sum()))
