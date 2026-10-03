import numpy as np
import pytest

mpl = pytest.importorskip("matplotlib")
mpl.use("Agg")
import matplotlib.pyplot as plt

from eprdating import Sediment, ToothSample, fit_dose_response, plot
from eprdating.beta import BetaGeometry
from eprdating.spectra import ComponentBasis, Spectrum


def _ld(B, c=336.8, w=0.6):
    x = (B - c) / (w * np.sqrt(3) / 2)
    return -x / (1 + x**2) ** 2


@pytest.fixture(autouse=True)
def _close():
    yield
    plt.close("all")


def test_spectra_and_spectrum():
    B = np.linspace(332, 342, 200)
    rng = np.random.default_rng(0)
    sp = [Spectrum(B, a * _ld(B)[None, :] + 0.05 * rng.standard_normal((3, B.size)), name=f"s{a}")
          for a in (1, 2, 3)]
    basis = ComponentBasis(B, {"c": _ld(B)})
    fits = [basis.fit(s.y) for s in sp]
    ax = plot.plot_spectra(sp, ["0 Gy", "10 Gy", "20 Gy"], fits=fits)
    assert len(ax.lines) == 6 and len(ax.texts) == 3
    ax = plot.plot_spectrum(sp[2], fit=fits[2])
    assert ax.get_legend() is not None


def test_dose_response_age_and_dose_rate():
    D = np.array([0, 20, 40, 60, 80.0])
    drc = fit_dose_response(D, 2 * (D + 10) + np.array([1, -1, 1, -1, 1.0]), "LIN", sigma=np.ones(5))
    ax = plot.plot_dose_response(drc, excluded=[(100, 50, 3)])
    assert any("De =" in t.get_text() for t in ax.texts)
    t = ToothSample(De=(5, 1), enamel_U=0, dentine_U=0, sediment=Sediment(U=1.5, Th=6, K=2.5, water=0.1),
                    beta=BetaGeometry(0, 0, 0), cosmic=0.18)
    plot.plot_dose_rate(t.age())
    ax = plot.plot_age_distribution(t.age_mc(n=200, seed=1), reference=(1, 3))
    assert ax.get_xlabel() == "Age (ka)"


def test_gamma_plots():
    from test_gamma import TRUE, synth

    from eprdating.gamma import IAEA_RGK_1, IAEA_RGTH_1, IAEA_RGU_1, analyse, calibrate_natural

    s = synth({"K": 0.03, "Ra226": 0.02, "U238": 0.02, "Th232": 0.08}, seed=3, live=86400 * 4)
    refs = {"U": synth({"Ra226": 4.0, "U238": 4.0}, seed=4), "Th": synth({"Th232": 8.0}, seed=5),
            "K": synth({"K": 0.448}, seed=6)}
    for x in (s, *refs.values()):
        calibrate_natural(x, TRUE, min_significance=5)
    r = analyse(s, 500, {"U": (refs["U"], IAEA_RGU_1), "Th": (refs["Th"], IAEA_RGTH_1),
                         "K": (refs["K"], IAEA_RGK_1)})
    plot.plot_gamma_spectrum(s)
    axs = plot.plot_gamma_lines(r)
    assert len(axs) == 4
