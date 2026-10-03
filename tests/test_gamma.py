import numpy as np
import pytest

from eprdating.gamma import (
    BQ_PER_KG,
    IAEA_RGK_1,
    IAEA_RGTH_1,
    IAEA_RGU_1,
    Calibration,
    GammaSpectrum,
    Reference,
    analyse,
    calibrate_natural,
    read_spectrum_txt,
    water_content,
)
from eprdating.gamma.comparative import LINES

TRUE = Calibration((2.2, 0.35175), (0.33, 0.0004))
N = 8192
# line energy -> counts per unit content per second (arbitrary efficiencies)
EFF = {ln.energy: np.exp(-ln.energy / 900) for ln in LINES}
for ln in LINES:
    for nb in ln.neighbours:
        EFF[nb] = 0.5 * np.exp(-nb / 900)
EXTRA = {1459.14: 0.0002}  # 228Ac 1459.1 keV under 40K: ~0.2 % of the K line in a soil
GROUP_OF = {ln.energy: ln.group for ln in LINES}
NEIGHBOUR_GROUP = {300.087: "Th232", 241.997: "Ra226", 964.766: "Th232", 1459.14: "Th232"}


ROOM = {"K": 0.01, "Ra226": 0.01, "U238": 0.01, "Th232": 0.01}  # laboratory background lines


def synth(activity: dict, live=86400.0, cal=TRUE, seed=0, continuum=2.0, name="s"):
    """Spectrum with line intensities ∝ group activity (counts/s per unit),
    on top of the laboratory background (present in every spectrum)."""
    activity = {g: activity.get(g, 0.0) + ROOM[g] for g in ROOM}
    ch = np.arange(1, N + 1, dtype=float)
    E = cal.energy(ch)
    lam = continuum * np.exp(-E / 700) * live / 86400 + 0.2
    for e, eff in {**EFF, **EXTRA}.items():
        g = GROUP_OF.get(e) or NEIGHBOUR_GROUP[e]
        a = activity.get(g, 0.0)
        if a <= 0:
            continue
        s = cal.sigma_keV(e) / cal.gain(cal.channel(e))
        area = a * eff * live
        lam = lam + area / (s * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((ch - cal.channel(e)) / s) ** 2)
    rng = np.random.default_rng(seed)
    return GammaSpectrum(rng.poisson(lam).astype(float), live, name=name, channels=ch)


def _write(tmp_path, spec, cal):
    lines = ["#", "# Sample name: test", "# Start time:    2024-08-01, 12:28:55", "# Real time (s): 86565.170",
             f"# Live time (s): {spec.live_time:.3f}", "#",
             "# Energy calibration coefficients ( E = sum(Ai * n**i) )",
             f"#     A0: {cal.coef[0]:.6f}", f"#     A1: {cal.coef[1]:.6f}", "#     A2: 0.000000",
             "# Energy unit: keV", "#", "# Channel data", "# n\tenergy(keV)\tcounts\trate(1/s)",
             "#" + "-" * 20]
    for n, c in zip(spec.channels, spec.counts):
        lines.append(f"{int(n)}\t{cal.energy(n):.3f}\t{int(c)}\t{c / spec.live_time:.6g}")
    p = tmp_path / "x_24h.txt"
    p.write_text("\n".join(lines), encoding="latin-1")
    return p


def test_reader(tmp_path):
    s = synth({"K": 1.0})
    p = _write(tmp_path, s, Calibration((0.0, 0.25)))
    r = read_spectrum_txt(p)
    assert r.live_time == 86400.0 and r.real_time == pytest.approx(86565.17)
    assert r.calibration.coef[:2] == (0.0, 0.25)
    np.testing.assert_array_equal(r.counts, s.counts)
    assert r.dead_time_fraction == pytest.approx(1 - 86400 / 86565.17)


def test_calibration_recovers_drift_and_drops_blended_line():
    s = synth({"K": 0.5, "Ra226": 1.0, "Th232": 3.0}, cal=Calibration((2.9, 0.35160), (0.33, 0.0004)))
    cal = calibrate_natural(s, TRUE)
    E = np.array([300.0, 1500.0, 2600.0])
    np.testing.assert_allclose(cal.energy(cal.channel(E)), E)
    true_ch = (E - 2.9) / 0.35160
    np.testing.assert_allclose(cal.energy(true_ch), E, atol=0.1)
    assert np.all(np.abs(cal.residuals_keV) < 0.4)
    assert cal.sigma_keV(1332) == pytest.approx(TRUE.sigma_keV(1332), rel=0.1)


def test_comparative_method_recovers_contents_and_equilibrium():
    ref_U, ref_Th, ref_K = 400.0, 800.0, 44.8
    # activity units per unit content (same for sample and reference)
    U = synth({"Ra226": ref_U * 0.01, "U238": ref_U * 0.01}, seed=1,
              cal=Calibration((2.34, 0.35173), (0.35, 0.0004)), name="U")
    Th = synth({"Th232": ref_Th * 0.01}, seed=2, cal=Calibration((2.35, 0.35177), (0.36, 0.0004)), name="Th")
    K = synth({"K": ref_K * 0.01}, seed=3, cal=Calibration((3.05, 0.35156), (0.37, 0.0004)), name="K")
    bg = synth({}, seed=4, cal=Calibration((2.43, 0.35154), (0.31, 0.0004)), continuum=0.5, name="bg")
    m = 600.0  # sample mass; references 500 g
    f = m / 500.0
    s = synth({"Ra226": 2.0 * 0.01 * f, "U238": 2.0 * 0.01 * f, "Th232": 8.0 * 0.01 * f, "K": 3.0 * 0.01 * f},
              seed=5, cal=Calibration((2.15, 0.35177), (0.30, 0.0004)), name="soil", live=86400 * 4)
    for sp in (U, Th, K, bg, s):
        calibrate_natural(sp, TRUE, min_significance=5)
    refs = {"U": (U, IAEA_RGU_1), "Th": (Th, IAEA_RGTH_1), "K": (K, IAEA_RGK_1)}
    r = analyse(s, m, refs, background=bg, groups=["K", "Ra226", "Th232"])
    for got, true in ((r.K, 3.0), (r.U, 2.0), (r.Th, 8.0)):
        assert abs(got.value - true) < 3 * got.sigma
        assert got.sigma / got.value < 0.05
    # reference uncertainty (2 % for RGU-1) is part of the total
    assert r.U.sigma / r.U.value >= 0.02
    sed = r.sediment(water=0.1)
    assert sed.K is r.K and sed.water == 0.1
    assert "Ra226" in r.summary()


def test_disequilibrium_shows_in_activity_ratio():
    U = synth({"Ra226": 4.0, "U238": 4.0}, seed=11, name="U")
    s = synth({"Ra226": 2.0, "U238": 4.0}, seed=12, name="soil")
    for sp in (U, s):
        calibrate_natural(sp, TRUE, min_significance=5)
    r = analyse(s, 500.0, {"U": (U, Reference("ref", {"U": 400.0}))}, groups=["Ra226", "U238"])
    assert r.activity_ratio_ra226_u238.value == pytest.approx(0.5, abs=0.1)


def test_water_content_and_specific_activities():
    w = water_content((700.0, 0.5), (600.0, 0.5))
    assert w.value == pytest.approx(100 / 600)
    assert w.sigma > 0
    # consistent with the IAEA certificates
    assert BQ_PER_KG["K"] * 44.8 == pytest.approx(14180, rel=0.01)
    assert BQ_PER_KG["U"] * 400 == pytest.approx(4941, rel=0.01)
    assert BQ_PER_KG["Th"] * 800 == pytest.approx(3250, rel=0.01)
