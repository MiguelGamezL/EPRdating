"""Environmental dose-rate building blocks.

* Infinite-matrix dose rates from U, Th, K concentrations (conversion factors).
* Water correction (Zimmerman 1971; Aitken 1985).
* Alpha efficiency (k-value).
* Cosmic dose rate (Prescott & Hutton 1994).

All dose rates are in Gy/ka.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from functools import cache
from importlib import resources

import numpy as np

from ._types import Value, ValueLike, as_value
from .history import History

#: Water-correction coefficients for alpha, beta and gamma radiation
#: (Zimmerman 1971; Aitken 1985): D_wet = D_dry / (1 + c * W).
WATER_COEFF = {"alpha": 1.49, "beta": 1.25, "gamma": 1.14}

#: Typical alpha efficiency for tooth enamel (Grün & Katzenberger-Apel 1994).
K_ENAMEL = Value(0.13, 0.02)

DEFAULT_FACTORS = "guerin_2011"


@dataclass(frozen=True)
class ConversionFactors:
    """Dose-rate conversion factors in Gy/ka per ppm (U, Th) or per % (K)."""

    key: str
    reference: str
    table: dict[str, dict[str, Value]]

    def get(self, nuclide: str, radiation: str) -> Value:
        try:
            return self.table[nuclide][radiation]
        except KeyError:
            return Value(0.0, 0.0)


@cache
def _raw_factors() -> dict:
    with resources.files("eprdating.data").joinpath("conversion_factors.json").open(encoding="utf-8") as fh:
        return json.load(fh)


def available_factor_sets() -> list:
    return [k for k in _raw_factors() if not k.startswith("_")]


def conversion_factors(key: str = DEFAULT_FACTORS) -> ConversionFactors:
    """Load a published set of conversion factors (see :func:`available_factor_sets`)."""
    raw = _raw_factors()
    if key not in raw or key.startswith("_"):
        raise KeyError(f"unknown factor set {key!r}; available: {available_factor_sets()}")
    table = {
        nuc: {rad: Value(*vals) for rad, vals in rads.items()}
        for nuc, rads in raw[key]["factors"].items()
    }
    return ConversionFactors(key, raw[key]["reference"], table)


#: share of the 230Th+226Ra segment carried by 226Ra (A&A 1998, Table 2 energies:
#: 230Th 4.58 / 0.013 / 0.0014 MeV, 226Ra 4.77 / 0.0038 / 0.0074 MeV for alpha / beta / gamma)
RA226_SHARE_OF_TH230 = {"alpha": 4.77 / (4.58 + 4.77), "beta": 0.0038 / (0.013 + 0.0038),
                        "gamma": 0.0074 / (0.0014 + 0.0074)}


def u_series_split(radiation: str) -> tuple[float, float]:
    """Fractions of the natural-U dose rate emitted before and from 226Ra on.

    "Before" is 238U, 234Th, 234Pa, 234U, 230Th and the 235U chain; "from
    226Ra" is 226Ra, 222Rn and its daughters. Returns ``(pre, post)``.
    """
    from .series import load_partition

    p = load_partition()[radiation]
    share = RA226_SHARE_OF_TH230[radiation]
    post = p["Th230"] * share + p["Rn222"]
    return 1.0 - post, post


def u_disequilibrium_factor(radiation: str, U: float, U_ra226: float | None) -> float:
    """U-equivalent content giving the dose rate of a chain whose 226Ra and
    daughters correspond to ``U_ra226`` ppm while the rest follows ``U``."""
    if U_ra226 is None:
        return U
    pre, post = u_series_split(radiation)
    return U * pre + U_ra226 * post


def matrix_dose_rates(
    U: float = 0.0,
    Th: float = 0.0,
    K: float = 0.0,
    factors: ConversionFactors | None = None,
    U_ra226: float | None = None,
) -> dict[str, float]:
    """Dry infinite-matrix alpha/beta/gamma dose rates (Gy/ka).

    ``U`` and ``Th`` in ppm, ``K`` in %. The Th chain is taken in secular
    equilibrium; so is the U chain unless ``U_ra226`` (226Ra and daughters
    expressed as ppm of U in equilibrium, as measured by gamma spectrometry
    through 214Pb/214Bi) differs from ``U`` (238U).
    """
    cf = factors or conversion_factors()
    out = {}
    for rad in ("alpha", "beta", "gamma"):
        out[rad] = (
            u_disequilibrium_factor(rad, U, U_ra226) * cf.get("U", rad).value
            + Th * cf.get("Th", rad).value
            + K * cf.get("K", rad).value
        )
    return out


def water_correction(rate: float, water: float, radiation: str) -> float:
    """Correct a dry dose rate for water content.

    ``water`` is the ratio mass of water / mass of dry material (e.g. 0.15).
    """
    if water < 0:
        raise ValueError("water content must be >= 0")
    return rate / (1.0 + WATER_COEFF[radiation] * water)


# --------------------------------------------------------------------------
# Cosmic dose rate
# --------------------------------------------------------------------------


@cache
def _fhj_table() -> dict[int, tuple]:
    path = resources.files("eprdating.data").joinpath("cosmic_prescott_hutton_FHJ.csv")
    with path.open(encoding="utf-8") as fh:
        return {int(r["geomag_lat_deg"]): (float(r["F"]), float(r["H"]), float(r["J"])) for r in csv.DictReader(fh)}


def geomagnetic_latitude(lat_deg: float, lon_deg: float) -> float:
    """Geomagnetic latitude (degrees) from geographic coordinates.

    Dipole approximation with the pole at 78.3°N, 291°E, as used by
    Prescott & Hutton (1994) and DRAC.
    """
    p = math.pi / 180
    s = 0.203 * math.cos(p * lat_deg) * math.cos(p * lon_deg - 291 * p) + 0.979 * math.sin(p * lat_deg)
    return math.asin(s) / p


def cosmic_dose_rate_sea_level(depth_m: float, density: float) -> float:
    """Cosmic dose rate at 55°N geomagnetic latitude and sea level (Gy/ka).

    ``depth_m`` is the overburden thickness (m) and ``density`` its mean
    density (g/cm³); their product is the shielding in hg/cm².
    """
    a = depth_m * density  # hg/cm^2
    if a < 1.67:
        return 0.0321 * a**4 - 0.135 * a**3 + 0.221 * a**2 - 0.207 * a + 0.295
    return (6072.0 / (((a + 11.6) ** 1.68 + 75.0) * (a + 212.0))) * math.exp(-0.00055 * a)


def cosmic_dose_rate(
    depth_m: float,
    density: float,
    lat_deg: float,
    lon_deg: float,
    altitude_m: float,
    rel_sigma: float = 0.10,
) -> Value:
    """Cosmic dose rate after Prescott & Hutton (1994), Gy/ka.

    Corrected for altitude and geomagnetic latitude with the F, J, H
    factors of Prescott & Stefan (1982). A 10 % relative uncertainty is
    assigned by default, as is common practice.
    """
    d0 = cosmic_dose_rate_sea_level(depth_m, density)
    lam = round(geomagnetic_latitude(lat_deg, lon_deg))
    F, H, J = _fhj_table()[lam]
    dc = d0 * (F + J * math.exp((altitude_m / 1000.0) / H))
    return Value(dc, dc * rel_sigma)


def cosmic_history(
    depth_m: History,
    density: float,
    lat_deg: float,
    lon_deg: float,
    altitude_m: float,
    rel_sigma: float = 0.10,
) -> History:
    """Cosmic dose-rate history from a burial-depth history (m, ka before present).

    Each segment gets :func:`cosmic_dose_rate` at its depth, with the
    ``rel_sigma`` uncertainty combined in quadrature with that of the depth
    (propagated through the local slope of the depth curve). A gradual
    burial is approximated by several short segments.
    """
    out = []
    for d in depth_m.values:
        dc = cosmic_dose_rate(d.value, density, lat_deg, lon_deg, altitude_m, rel_sigma)
        slope = 0.0
        if d.sigma > 0:
            h = max(1e-3, 0.01 * d.sigma)
            lo, hi = max(d.value - h, 0.0), d.value + h
            slope = (cosmic_dose_rate(hi, density, lat_deg, lon_deg, altitude_m).value
                     - cosmic_dose_rate(lo, density, lat_deg, lon_deg, altitude_m).value) / (hi - lo)
        out.append(Value(dc.value, math.hypot(dc.sigma, slope * d.sigma)))
    return History(out, depth_m.breaks)


# --------------------------------------------------------------------------
# Convenience for a sediment/soil source
# --------------------------------------------------------------------------


@dataclass
class Sediment:
    """Radionuclide content of a sediment or soil (U, Th in ppm, K in %).

    ``U_ra226``: 226Ra and its daughters as ppm of U in equilibrium, when it
    differs from the 238U content ``U`` (e.g. from gamma spectrometry, 214Pb
    and 214Bi lines vs 234Th and 234mPa). ``None`` means equilibrium.
    ``water`` (mass of water / dry mass) may be a
    :class:`~eprdating.history.History` when it changed during burial; its
    first value is the present-day one.
    """

    U: ValueLike = 0.0
    Th: ValueLike = 0.0
    K: ValueLike = 0.0
    water: ValueLike | History = 0.0
    U_ra226: ValueLike | None = None

    def dose_rate(self, radiation: str, factors: ConversionFactors | None = None, values=None) -> float:
        """Wet dose rate for one radiation type, with the present-day water.
        ``values`` overrides the inputs (used by the Monte Carlo engine)."""
        if values is None and isinstance(self.water, History):
            values = {k: as_value(getattr(self, k)).value for k in ("U", "Th", "K")}
            values["water"] = self.water.nominal()[0]
        v = values or {k: as_value(getattr(self, k)).value for k in ("U", "Th", "K", "water")}
        ura = v.get("U_ra226", None if self.U_ra226 is None else as_value(self.U_ra226).value)
        dry = matrix_dose_rates(v["U"], v["Th"], v["K"], factors, U_ra226=ura)[radiation]
        return water_correction(dry, v["water"], radiation)


def samples_from(obj, names, rng: np.random.Generator, n: int) -> dict[str, np.ndarray]:
    """Draw Gaussian samples for the listed attributes of ``obj``."""
    return {k: as_value(getattr(obj, k)).sample(rng, n) for k in names}
