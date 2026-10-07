"""Alpha escape at layer surfaces (eprdating.alpha, ToothSample(alpha_escape=True))."""

import math

import numpy as np
import pytest

from eprdating import Sediment, ToothLayers, ToothSample
from eprdating.alpha import escape_fractions, slab_fractions, surface_fraction
from eprdating.beta import BetaGeometry


def test_escape_from_a_surface_matches_straight_isotropic_tracks():
    """e(u) = [(1-u) + u ln u]/2 is the track-length fraction leaving through a
    plane at distance u·R, for isotropic straight tracks of length R."""
    rng = np.random.default_rng(0)
    mu = rng.uniform(-1, 1, 400_000)  # direction cosine towards the surface
    for u in (0.05, 0.3, 0.7):
        out = np.where(mu > u, 1 - u / np.where(mu > 0, mu, 1), 0.0)
        e = 0.5 * ((1 - u) + u * math.log(u))
        assert out.mean() == pytest.approx(e, abs=2e-3)


def test_layer_averages():
    R = 20.0
    # whole layer: R/(8T) lost per surface, and as much received from each side
    own, outer, inner = slab_fractions(1000, 0, 0, R)
    assert own == pytest.approx(1 - R / 4000) and outer == inner == pytest.approx(R / 8000)
    # stripping at least one range removes both effects
    assert slab_fractions(1000, 25, 25, R) == pytest.approx((1.0, 0.0, 0.0))
    # stripping one side only
    own, outer, inner = slab_fractions(1000, 25, 0, R)
    assert outer == 0.0 and inner == pytest.approx(R / 8 / 975)
    assert surface_fraction(0, 0, R) == pytest.approx(0.5)  # at the surface half the directions leave
    with pytest.raises(ValueError):
        slab_fractions(100, 60, 60, R)


def test_segment_fractions():
    f = escape_fractions(300)
    # the post-radon alphas are the most energetic: they escape most
    assert f["own"]["Rn222"] < f["own"]["U238"] < 1
    assert f["own"]["Rn222"] == pytest.approx(1 - 2 * f["outer"]["Rn222"])
    g = escape_fractions(300, outer_material=ToothLayers(300).sediment,
                         inner_material=ToothLayers(300).dentine)
    # alphas from media with a shorter mass range bring less into the enamel
    assert g["inner"]["Rn222"] < g["outer"]["Rn222"] < f["outer"]["Rn222"]


def _tooth(strip, **kw):
    return ToothSample(
        De=200, enamel_U=1.0, dentine_U=20.0, sediment=Sediment(U=3.0, Th=10.0, K=1.0, water=0.1),
        beta=ToothLayers(enamel_um=400, strip_outer_um=strip, strip_inner_um=strip), cosmic=0.15,
        alpha_efficiency="energy", **kw)


def test_tooth_sample():
    plain, esc = _tooth(0).age(), _tooth(0, alpha_escape=True).age()
    assert {"dentine alpha", "sediment alpha"} <= set(esc.components)
    assert esc.components["enamel alpha"] < plain.components["enamel alpha"]
    assert esc.age < plain.age  # here the alphas entering outweigh those leaving
    # with the surfaces stripped beyond the alpha range nothing changes
    assert _tooth(45, alpha_escape=True).age().age == pytest.approx(_tooth(45).age().age, rel=1e-9)
    with pytest.raises(TypeError):
        ToothSample(De=100, enamel_U=1, dentine_U=10, sediment=Sediment(U=1, Th=1, K=1),
                    beta=BetaGeometry(internal=0.7, dentine=0.1, external=0.1), cosmic=0.1,
                    alpha_escape=True).age()
    mc = _tooth(0, alpha_escape=True).age_mc(n=50, seed=1)
    assert mc.samples.size == 50
