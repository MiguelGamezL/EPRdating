from itertools import pairwise

import pytest

from eprdating.alpha import alpha_range, k_ratio, natural_u_k_ratio, segment_k_ratios


def test_reference_energy_has_ratio_one():
    assert k_ratio(5.3) == pytest.approx(1.0)
    assert k_ratio(6.0, e_ref=6.0) == pytest.approx(1.0)


def test_efficiency_grows_with_energy():
    es = [4.0, 4.8, 5.5, 6.5, 7.7]
    ks = [k_ratio(e) for e in es]
    assert all(a < b for a, b in pairwise(ks))


def test_range_plausible():
    # a 5.3 MeV alpha travels a few tens of µm in enamel (density ~3 g/cm3)
    um = alpha_range(5.3) / 3.0 * 1e4
    assert 10 < um < 30


def test_segments():
    seg = segment_k_ratios()
    assert seg["U238"] < seg["U234"] < seg["Rn222"]
    assert seg["Rn222"] == pytest.approx(seg["Pa231"], rel=0.02)
    assert 1.0 < natural_u_k_ratio() < 1.05
