"""Generate the figures of the documentation (docs/img) from synthetic data.

    python docs/make_figures.py

Everything is simulated so the figures can be published: a dose series of
an enamel sample (De = 40 Gy), its dose-response curve, a tooth age with
Monte Carlo, and HPGe spectra of a soil and the IAEA reference materials.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from eprdating import (
    LinearUptake,
    Sediment,
    ToothLayers,
    ToothSample,
    fit_dose_response,
    plot,
)
from eprdating.gamma import (
    IAEA_RGK_1,
    IAEA_RGTH_1,
    IAEA_RGU_1,
    Calibration,
    GammaSpectrum,
    analyse,
    auto_calibrate,
)
from eprdating.gamma.comparative import LINES
from eprdating.spectra import ComponentBasis, Spectrum, pseudo_modulation

OUT = Path(__file__).parent / "img"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(7)


def save(ax_or_fig, name):
    fig = ax_or_fig if hasattr(ax_or_fig, "savefig") else np.atleast_1d(ax_or_fig)[0].figure
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=110, facecolor=fig.get_facecolor())
    plt.close(fig)


# ------------------------------------------------------------------ EPR
B = np.linspace(332.0, 342.0, 400)


def co2_like(B):
    """Asymmetric derivative line resembling the enamel CO2- signal."""
    def ld(c, w):
        x = (B - c) / (w * np.sqrt(3) / 2)
        return -x / (1 + x**2) ** 2
    y = ld(336.9, 0.25) + 0.35 * ld(337.35, 0.35)
    return pseudo_modulation(B, y / np.ptp(y), 0.1)


shape = co2_like(B)
doses = np.array([0, 25, 50, 100, 200, 400, 800, 1600])
De_true, D0, Imax = 40.0, 1200.0, 100.0
spectra = []
for d in doses:
    amp = Imax * (1 - np.exp(-(d + De_true) / D0)) * rng.normal(1, 0.02)  # 2 % aliquot scatter
    scans = amp * shape[None, :] + 0.6 * rng.standard_normal((4, B.size)) + rng.normal(0, 0.3)
    spectra.append(Spectrum(B, scans, name=f"{d} Gy", freq_GHz=9.43, power_mW=2.0))

basis = ComponentBasis(B, {"CO2-": shape}, baseline_order=1)
fits = [basis.fit(s.y, nonnegative=False) for s in spectra]
amps = np.array([f.amplitudes["CO2-"] for f in fits])
errs = np.array([f.errors["CO2-"] for f in fits])

save(plot.plot_spectra(spectra, [f"{d} Gy" for d in doses], fits=fits, title="Additive-dose series"), "spectra")
save(plot.plot_spectrum(spectra[3], title="One aliquot, four scans"), "spectrum")
drc = fit_dose_response(doses, amps, "SSE", sigma=errs)
save(plot.plot_dose_response(drc, title="Dose-response (SSE)"), "dose_response")

# ------------------------------------------------------------------ age
tooth = ToothSample(
    De=(drc.De, drc.De_sigma), enamel_U=(0.6, 0.06), dentine_U=(12.0, 1.2),
    sediment=Sediment(U=(2.0, 0.1), Th=(7.0, 0.3), K=(1.2, 0.05), water=(0.15, 0.05)),
    beta=ToothLayers(enamel_um=(1100, 100), dentine_um=2000, strip_outer_um=(50, 20), strip_inner_um=(50, 20)),
    cosmic=(0.15, 0.02), uptake_enamel=LinearUptake(), uptake_dentine=LinearUptake(),
)
save(plot.plot_dose_rate(tooth.age()), "dose_rate")
save(plot.plot_age_distribution(tooth.age_mc(n=2000, seed=3), title="Monte Carlo age"), "age")

# ------------------------------------------------------------------ gamma
TRUE = Calibration((2.2, 0.35175), (0.33, 0.0004))
GROUP_OF = {ln.energy: ln.group for ln in LINES}
NEIGH = {300.087: "Th232", 241.997: "Ra226", 964.766: "Th232"}
EFF = {e: np.exp(-e / 900) for e in [*GROUP_OF, *NEIGH]}
ROOM = {"K": 0.01, "Ra226": 0.01, "U238": 0.01, "Th232": 0.01}


def gamma(activity, live=86400.0, cal=TRUE, name="s"):
    ch = np.arange(1, 8193, dtype=float)
    E = cal.energy(ch)
    lam = (300.0 * np.exp(-E / 450) + 25.0 * np.exp(-E / 1500)) * live / 86400 + 0.5  # Compton continuum
    for e, eff in EFF.items():
        g = GROUP_OF.get(e) or NEIGH[e]
        a = activity.get(g, 0.0) + ROOM[g]
        s = cal.sigma_keV(e) / cal.gain(cal.channel(e))
        lam = lam + 3 * a * eff * live / (s * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((ch - cal.channel(e)) / s) ** 2)
    return GammaSpectrum(rng.poisson(lam).astype(float), live, name=name, channels=ch)


soil = gamma({"K": 0.03, "Ra226": 0.02, "U238": 0.02, "Th232": 0.08}, live=86400 * 2, name="soil")
refs = {"U": gamma({"Ra226": 4.0, "U238": 4.0}, name="RGU-1"), "Th": gamma({"Th232": 8.0}, name="RGTh-1"),
        "K": gamma({"K": 0.448}, name="RGK-1")}
bg = gamma({}, name="background")
for s in (soil, bg, *refs.values()):
    auto_calibrate(s)
res = analyse(soil, 500.0, {"U": (refs["U"], IAEA_RGU_1), "Th": (refs["Th"], IAEA_RGTH_1),
                            "K": (refs["K"], IAEA_RGK_1)}, background=bg)
save(plot.plot_gamma_spectrum(soil, title="Soil sample, HPGe"), "gamma_spectrum")
save(plot.plot_gamma_lines(res), "gamma_lines")
print(drc.summary())
print(res.summary())
