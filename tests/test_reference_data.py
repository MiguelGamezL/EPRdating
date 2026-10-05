"""Reference data checked against their sources.

Conversion factors: the central values of the three sets are compared with
the papers and with an independent transcription, the
``BaseDataSet.ConversionFactors`` table of the R package Luminescence
(github.com/R-Lum/Luminescence, data/). Guérin et al. (2011, Ancient TL 29,
5-8): Th and K checked against Tables 1 and 3 of the paper; natural U is the
abundance-weighted sum of its Table 2 and agrees with Luminescence. Guérin et
al. give no uncertainties; EPRdating carries the relative uncertainties of
Adamiec & Aitken (1998), as DRAC does.

Cosmic dose rate: compared with ``calc_CosmicDoseRate`` of Luminescence, an
implementation of Prescott & Hutton (1994) independent of DRAC.
"""

import math

import pytest

from eprdating import conversion_factors
from eprdating.dose_rate import _fhj_table, cosmic_dose_rate, cosmic_dose_rate_sea_level, geomagnetic_latitude

PUBLISHED = {
    "adamiec_aitken_1998": {"U": (2.78, 0.146, 0.113), "Th": (0.732, 0.0273, 0.0476), "K": (None, 0.782, 0.243)},
    "guerin_2011": {"U": (2.795, 0.1457, 0.1116), "Th": (0.7375, 0.0277, 0.0479), "K": (None, 0.7982, 0.2491)},
    "liritzis_2013": {"U": (2.793, 0.1459, 0.1118), "Th": (0.7375, 0.0275, 0.0481), "K": (None, 0.8011, 0.2498)},
}
LIRITZIS_SIGMA = {"U": (0.011, 0.0004, 0.0002), "Th": (0.0026, 0.0009, 0.0002), "K": (None, 0.0073, 0.0048)}


@pytest.mark.parametrize("key", PUBLISHED)
def test_conversion_factors_match_the_sources(key):
    cf = conversion_factors(key)
    for nuc, vals in PUBLISHED[key].items():
        for rad, v in zip(("alpha", "beta", "gamma"), vals, strict=True):
            if v is not None:
                assert cf.get(nuc, rad).value == pytest.approx(v, abs=1e-6), (nuc, rad)


def test_liritzis_uncertainties():
    cf = conversion_factors("liritzis_2013")
    for nuc, vals in LIRITZIS_SIGMA.items():
        for rad, s in zip(("alpha", "beta", "gamma"), vals, strict=True):
            if s is not None:
                assert cf.get(nuc, rad).sigma == pytest.approx(s, abs=1e-6), (nuc, rad)


# --- cosmic dose rate against Luminescence::calc_CosmicDoseRate ------------------
def _lum_d0(hg):
    d0 = 6072 / (((hg + 11.6) ** 1.68 + 75) * (hg + 212)) * math.exp(-0.00055 * hg)
    t = hg * 100
    if t < 167:
        soft = -6e-8 * t**3 + 2e-5 * t**2 - 0.0025 * t + 0.2969 if t < 40 else 2e-6 * t**2 - 0.0008 * t + 0.2535
        d0 = max(soft, d0)
    return d0


def _lum_gml(lat, lon):
    r = math.pi / 180
    return math.asin(0.203 * math.cos(r * lat) * math.cos(r * (lon - 291)) + 0.979 * math.sin(r * lat)) / r


def _lum_fjh(g):
    F = -7e-7 * g**3 - 8e-5 * g**2 - 0.0009 * g + 0.3988 if g < 36.5 else -0.0001 * g + 0.2347
    J = 5e-6 * g**3 - 5e-5 * g**2 + 0.0026 * g + 0.5177 if g < 34 else 0.0005 * g + 0.7388
    H = -3e-6 * g**3 - 5e-5 * g**2 - 0.0031 * g + 4.398 if g < 36 else 0.0002 * g + 4.0914
    return F, J, H


@pytest.mark.parametrize("hg", [1.7, 2.0, 4.0, 10.0, 30.0, 100.0])
def test_cosmic_deep_burial_identical(hg):
    assert cosmic_dose_rate_sea_level(hg, 1.0) == pytest.approx(_lum_d0(hg), rel=1e-4)


@pytest.mark.parametrize("hg", [0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 1.0, 1.5])
def test_cosmic_shallow_burial_within_seven_percent(hg):
    # two fits to the soft component of Prescott & Hutton (1994): DRAC's
    # polynomial (used here) and Luminescence's; they differ by up to 6.7 %
    # between 0.2 and 1 hg/cm2, inside the usual 10 % uncertainty
    assert cosmic_dose_rate_sea_level(hg, 1.0) == pytest.approx(_lum_d0(hg), rel=0.07)


@pytest.mark.parametrize("lat,lon", [(4.6, -74.1), (-26.0, 27.7), (45.0, 10.0), (-7.0, 111.0), (60.0, 30.0)])
def test_geomagnetic_latitude(lat, lon):
    assert geomagnetic_latitude(lat, lon) == pytest.approx(_lum_gml(lat, lon), abs=1e-9)


@pytest.mark.parametrize("g", range(0, 91, 5))
def test_fjh_table_matches_an_independent_reading(g):
    # two readings of Fig. 3 of Prescott & Stephan (1982): DRAC's table (used
    # here) and Luminescence's fits; F differs by up to 4 % near 15 degrees and
    # J by 3 % at high latitude, where Luminescence extrapolates a line and the
    # table stays flat; the dose rate itself agrees within 2.5 % (test below)
    F, H, J = _fhj_table()[g]
    lF, lJ, lH = _lum_fjh(g)
    assert F == pytest.approx(lF, rel=0.05)
    assert J == pytest.approx(lJ, rel=0.05)
    assert H == pytest.approx(lH, rel=0.01)


@pytest.mark.parametrize("lat,lon", [(4.6, -74.1), (-26.0, 27.7), (45.0, 10.0), (-7.0, 111.0)])
@pytest.mark.parametrize("alt", [0.0, 1000.0, 2600.0])
def test_cosmic_dose_rate_end_to_end(lat, lon, alt):
    F, J, H = _lum_fjh(abs(_lum_gml(lat, lon)))
    lum = _lum_d0(2.0) * (F + J * math.exp(alt / 1000 / H))
    assert cosmic_dose_rate(1.0, 2.0, lat, lon, alt).value == pytest.approx(lum, rel=0.025)
