"""Age panel: De, tooth, sediment and site → dose rate and age (Monte Carlo)."""

from __future__ import annotations

import html
import json
import warnings

import ipywidgets as w

from .. import __version__
from .._types import as_value
from ..age import ToothSample
from ..dose_rate import Sediment, available_factor_sets, cosmic_dose_rate
from ..history import History
from ..onegroup import ToothLayers
from ..uptake import USModel
from ._common import (
    STYLE,
    button,
    checkbox,
    download_link,
    message,
    read_value,
    set_value,
    show_figure,
    value_field,
)

UPTAKE = {"EU (p = -1)": -1.0, "LU (p = 0)": 0.0, "US (give p)": None}


class AgePanel:
    """Dose rate and age of a tooth (one tab of :func:`eprdating.gui.app`).
    :meth:`tooth` builds the :class:`~eprdating.ToothSample` from the form,
    :meth:`compute` gives the nominal and Monte Carlo ages."""

    def __init__(self) -> None:
        self.result = None
        self.mc = None
        self._build()

    def _build(self) -> None:
        v = value_field
        self.De = v("De", 0.0, 0.0, "Gy")
        self.de_note = w.HTML("<small>Filled from the De tab after a fit.</small>")

        self.enamel_U = v("Enamel U", 0.5, 0.05, "µg/g")
        self.dentine_U = v("Dentine U", 10.0, 1.0, "µg/g")
        self.enamel_um = v("Enamel thickness", 1000.0, 100.0, "µm")
        self.dentine_um = v("Dentine thickness", 2000.0, 0.0, "µm")
        self.strip_out = v("Removed, outer side", 50.0, 20.0, "µm")
        self.strip_in = v("Removed, inner side", 50.0, 20.0, "µm")
        self.dentine_water = v("Dentine water", 0.05, 0.02, "")
        self.enamel_water = v("Enamel water", 0.0, 0.0, "")
        self.k_alpha = v("Alpha efficiency k", 0.13, 0.02, "")
        self.alpha_mode = w.Dropdown(options=["constant", "energy"], value="constant",
                                     description="k with alpha energy", style=STYLE)
        self.up_e = w.Dropdown(options=list(UPTAKE), value="LU (p = 0)", description="Uptake enamel", style=STYLE)
        self.p_e = w.FloatText(value=0.0, description="p", style={"description_width": "14px"},
                               layout=w.Layout(width="100px"))
        self.up_d = w.Dropdown(options=list(UPTAKE), value="LU (p = 0)", description="Uptake dentine", style=STYLE)
        self.p_d = w.FloatText(value=0.0, description="p", style={"description_width": "14px"},
                               layout=w.Layout(width="100px"))
        self.r_e = v("²³⁴U/²³⁸U enamel", 1.0, 0.0, "")
        self.r_d = v("²³⁴U/²³⁸U dentine", 1.0, 0.0, "")
        self.rn_d = v("Radon loss dentine", 0.0, 0.0, "(0-1)")
        self.by_segment = checkbox("beta attenuation per U-series segment (as ROSY; untick for DATA's single factor)")
        self.alpha_escape = checkbox("alpha escape at the enamel surfaces (as ROSY; matters without stripping)",
                                     value=False)
        self.factors = w.Dropdown(options=available_factor_sets(), value="guerin_2011",
                                  description="Conversion factors", style=STYLE)

        self.sed_U = v("Sediment U", 2.0, 0.1, "µg/g")
        self.sed_Th = v("Sediment Th", 6.0, 0.3, "µg/g")
        self.sed_K = v("Sediment K", 1.0, 0.05, "%")
        self.sed_Ura = v("U from ²²⁶Ra (0: equil.)", 0.0, 0.0, "µg/g")
        self.sed_water = v("Sediment water", 0.10, 0.05, "")
        self.sed_note = w.HTML("<small>Filled from the Sediment tab after an analysis.</small>")
        self.gamma_mode = w.Dropdown(options=["from the sediment", "measured in situ"], value="from the sediment",
                                     description="Gamma", style=STYLE)
        self.gamma = v("In-situ gamma", 0.0, 0.0, "Gy/ka")

        self.cosmic_mode = w.Dropdown(options=["from depth and site", "given"], value="from depth and site",
                                      description="Cosmic", style=STYLE)
        f = lambda d, x: w.FloatText(value=x, description=d, style=STYLE, layout=w.Layout(width="230px"))
        self.depth = f("Depth (m)", 1.0)
        self.density = f("Overburden density (g/cm³)", 1.9)
        self.lat = f("Latitude (°)", 4.6)
        self.lon = f("Longitude (°)", -74.1)
        self.alt = f("Altitude (m)", 2600.0)
        self.cosmic = v("Cosmic", 0.15, 0.015, "Gy/ka")

        self.n_mc = w.BoundedIntText(value=1000, min=0, max=100000, description="Monte Carlo draws",
                                     style=STYLE, layout=w.Layout(width="230px"))
        self.run = button("Compute age", "play", primary=True, width="160px")
        self.run.on_click(lambda _: self._guard(self.compute))
        self.status = w.HTML()
        self.summary = w.HTML()
        self.out = w.Output()
        self.export_btn = button("Prepare downloads", "download")
        self.export_btn.on_click(lambda _: self._guard(self._export))
        self.export_html = w.HTML()

        h = lambda t: w.HTML(f"<h4 style='margin:10px 0 4px'>{t}</h4>")
        self.widget = w.VBox([
            h("1. Equivalent dose"), self.De, self.de_note,
            h("2. Tooth"),
            w.HBox([w.VBox([self.enamel_U, self.dentine_U, self.r_e, self.r_d, self.rn_d]),
                    w.VBox([self.enamel_um, self.dentine_um, self.strip_out, self.strip_in])]),
            w.HBox([w.VBox([self.dentine_water, self.enamel_water]),
                    w.VBox([self.k_alpha, self.alpha_mode])]),
            w.HBox([self.up_e, self.p_e]), w.HBox([self.up_d, self.p_d]),
            self.by_segment, self.alpha_escape, self.factors,
            h("3. Sediment and site"),
            w.HBox([w.VBox([self.sed_U, self.sed_Th, self.sed_K]), w.VBox([self.sed_Ura, self.sed_water])]),
            self.sed_note,
            w.HBox([self.gamma_mode, self.gamma]),
            self.cosmic_mode,
            w.HBox([w.VBox([self.depth, self.density]), w.VBox([self.lat, self.lon, self.alt])]),
            self.cosmic,
            h("4. Age"),
            w.HBox([self.n_mc, self.run]), self.status, self.summary, self.out,
            h("5. Export"),
            w.HBox([self.export_btn, self.export_html]),
        ])

    def _ipython_display_(self):
        from IPython.display import display

        display(self.widget)

    def _guard(self, fn) -> None:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fn()
        except Exception as e:  # noqa: BLE001 - shown in the panel
            self.status.value = message(f"{type(e).__name__}: {e}", "error")

    # ---- links from the other tabs ------------------------------------------
    def set_De(self, value: float, sigma: float) -> None:
        set_value(self.De, value, sigma)

    def set_sediment(self, sed: Sediment) -> None:
        for box, x in ((self.sed_U, sed.U), (self.sed_Th, sed.Th), (self.sed_K, sed.K)):
            q = as_value(x)
            set_value(box, q.value, q.sigma)
        if sed.U_ra226 is not None:
            q = as_value(sed.U_ra226)
            set_value(self.sed_Ura, q.value, q.sigma)
        else:
            set_value(self.sed_Ura, 0.0, 0.0)
        if not isinstance(sed.water, History):
            q = as_value(sed.water)
            set_value(self.sed_water, q.value, q.sigma)

    # ---- model ---------------------------------------------------------------
    def _uptake(self, dd, p) -> USModel:
        val = UPTAKE[dd.value]
        return USModel(float(p.value) if val is None else val)

    def tooth(self) -> ToothSample:
        rv = read_value
        de = rv(self.De)
        if de[0] <= 0:
            raise ValueError("give a positive De")
        ura = rv(self.sed_Ura)
        sed = Sediment(U=rv(self.sed_U), Th=rv(self.sed_Th), K=rv(self.sed_K), water=rv(self.sed_water),
                       U_ra226=ura if ura[0] > 0 else None)
        if self.cosmic_mode.value == "given":
            cosmic = rv(self.cosmic)
        else:
            cosmic = cosmic_dose_rate(float(self.depth.value), float(self.density.value), float(self.lat.value),
                                      float(self.lon.value), float(self.alt.value))
        gamma = rv(self.gamma) if self.gamma_mode.value == "measured in situ" else None
        geo = ToothLayers(enamel_um=rv(self.enamel_um), dentine_um=rv(self.dentine_um),
                          strip_outer_um=rv(self.strip_out), strip_inner_um=rv(self.strip_in))
        return ToothSample(
            De=de, enamel_U=rv(self.enamel_U), dentine_U=rv(self.dentine_U), sediment=sed, beta=geo,
            cosmic=cosmic, gamma=gamma, k_alpha=rv(self.k_alpha), alpha_efficiency=self.alpha_mode.value,
            dentine_water=rv(self.dentine_water), enamel_water=rv(self.enamel_water),
            uptake_enamel=self._uptake(self.up_e, self.p_e), uptake_dentine=self._uptake(self.up_d, self.p_d),
            u234_u238_enamel=rv(self.r_e), u234_u238_dentine=rv(self.r_d), radon_loss_dentine=rv(self.rn_d),
            factors=self.factors.value, beta_by_segment=bool(self.by_segment.value),
            alpha_escape=bool(self.alpha_escape.value),
        )

    def compute(self):
        tooth = self.tooth()
        self.result = tooth.age()
        n = int(self.n_mc.value)
        self.mc = tooth.age_mc(n=n, seed=1) if n > 0 else None
        text = self.result.summary() + ("\n\n" + self.mc.summary() if self.mc is not None else "")
        self.summary.value = f"<pre style='font-size:12px'>{html.escape(text)}</pre>"
        import matplotlib.pyplot as plt

        from .. import plot

        fig, axes = plt.subplots(1, 2 if self.mc is not None else 1, figsize=(12, 3.8))
        axes = list(axes) if self.mc is not None else [axes]
        plot.plot_dose_rate(self.result, ax=axes[0])
        if self.mc is not None:
            plot.plot_age_distribution(self.mc, ax=axes[1])
        fig.tight_layout()
        show_figure(self.out, fig)
        self.status.value = message(f"Age = {self.result.age:.4g} ka", "ok")
        return self.result

    def settings(self) -> dict:
        rv = read_value
        boxes = {k: getattr(self, k) for k in (
            "De", "enamel_U", "dentine_U", "enamel_um", "dentine_um", "strip_out", "strip_in", "dentine_water",
            "enamel_water", "k_alpha", "r_e", "r_d", "rn_d", "sed_U", "sed_Th", "sed_K", "sed_Ura", "sed_water",
            "gamma", "cosmic")}
        return {
            "eprdating": __version__,
            **{k: rv(b) for k, b in boxes.items()},
            "alpha_efficiency": self.alpha_mode.value, "uptake_enamel": [self.up_e.value, self.p_e.value],
            "uptake_dentine": [self.up_d.value, self.p_d.value], "beta_by_segment": self.by_segment.value,
            "alpha_escape": self.alpha_escape.value,
            "factors": self.factors.value, "gamma_mode": self.gamma_mode.value,
            "cosmic_mode": self.cosmic_mode.value,
            "site": {"depth_m": self.depth.value, "density": self.density.value, "lat": self.lat.value,
                     "lon": self.lon.value, "altitude_m": self.alt.value},
            "monte_carlo": self.n_mc.value,
            "age_ka": None if self.result is None else self.result.age,
        }

    def _export(self) -> None:
        if self.result is None:
            raise ValueError("compute the age first")
        text = self.result.summary() + ("\n\n" + self.mc.summary() if self.mc is not None else "")
        self.export_html.value = (download_link("eprdating_age.txt", text, "result (text)")
                                  + download_link("eprdating_age_settings.json", json.dumps(self.settings(), indent=2),
                                                  "settings (JSON)"))
