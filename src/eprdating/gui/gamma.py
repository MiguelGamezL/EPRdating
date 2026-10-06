"""Sediment panel: HPGe spectra → U, Th, K → infinite-matrix dose rates."""

from __future__ import annotations

import html
import json
import re
import warnings
from collections.abc import Callable

import ipywidgets as w

from .. import __version__
from .._types import as_value
from ..dose_rate import Sediment, matrix_dose_rates, water_correction
from ..gamma import NATURAL_LINES, Reference, analyse, auto_calibrate, read_gamma
from ._common import (
    STYLE,
    FilePool,
    button,
    checkbox,
    clear_upload,
    colab_hint,
    download_link,
    message,
    read_value,
    show_figure,
    table_html,
    uploaded,
    value_field,
)

SUFFIXES = (".spe", ".chn", ".n42", ".xml", ".cnf", ".spc", ".iec", ".txt", ".csv")
NONE = "(none)"
#: words in a file name that suggest its role
_ROLE_HINTS = {
    "background": r"fondo|background|bkg|blank|empty|vac",
    "U": r"rgu|(^|[_\-\s])u([_\-\s.]|$)",
    "Th": r"rgth|(^|[_\-\s])th([_\-\s.]|$)",
    "K": r"rgk|(^|[_\-\s])k([_\-\s.]|$)",
}


def guess_role(name: str) -> str | None:
    """"background", "U", "Th", "K" or None, from a file name."""
    low = name.lower()
    for role, pat in _ROLE_HINTS.items():
        if re.search(pat, low):
            return role
    return None


class GammaPanel:
    """Sediment U, Th and K by the comparative method (one tab of
    :func:`eprdating.gui.app`). Drive it from code with :meth:`load`,
    :meth:`set_roles` and :meth:`analyse`."""

    def __init__(self, on_sediment: Callable[[Sediment], None] | None = None) -> None:
        self.on_sediment = on_sediment
        self.pool = FilePool()
        self.result = None
        self.sediment = None
        self._build()

    def _build(self) -> None:
        self.upload = w.FileUpload(accept=",".join(SUFFIXES), multiple=True, description="Upload spectra",
                                   layout=w.Layout(width="180px"))
        self.upload.observe(self._on_upload, names="value")
        self.folder = w.Text(placeholder="folder with the HPGe spectra", layout=w.Layout(width="420px"))
        self.load_folder = button("Load folder", "folder-open", width="140px")
        self.load_folder.on_click(lambda _: self._guard(self._load_folder))
        self.live_time = w.FloatText(value=0.0, description="Live time for plain column files (s)", style=STYLE,
                                     layout=w.Layout(width="320px"))

        dd = lambda d: w.Dropdown(options=[NONE], description=d, style=STYLE,
                                  layout=w.Layout(width="420px"))
        self.sample = dd("Sample")
        self.background = dd("Background")
        self.ref = {el: dd(f"{name} ({el})") for el, name in
                    (("U", "IAEA-RGU-1"), ("Th", "IAEA-RGTh-1"), ("K", "IAEA-RGK-1"))}
        self.sample_mass = w.FloatText(value=500.0, description="Sample dry mass (g)", style=STYLE,
                                       layout=w.Layout(width="240px"))
        self.ref_mass = w.FloatText(value=500.0, description="Reference masses (g)", style=STYLE,
                                    layout=w.Layout(width="240px"))
        self.ref_content = {"U": value_field("RGU-1 U", 400.0, 8.0, "µg/g"),
                            "Th": value_field("RGTh-1 Th", 800.0, 16.0, "µg/g"),
                            "K": value_field("RGK-1 K", 44.8, 0.3, "%")}
        self.water = value_field("Water (water/dry mass)", 0.0, 0.0, "")
        self.equilibrium = checkbox("U chain in equilibrium (U from the 226Ra daughters)")
        self.analyse_btn = button("Analyse", "play", primary=True, width="140px")
        self.analyse_btn.on_click(lambda _: self._guard(self.analyse))
        self.status = w.HTML()
        self.summary = w.HTML()
        self.out = w.Output()
        self.export_btn = button("Prepare downloads", "download")
        self.export_btn.on_click(lambda _: self._guard(self._export))
        self.export_html = w.HTML()

        h = lambda t: w.HTML(f"<h4 style='margin:10px 0 4px'>{t}</h4>")
        self.widget = w.VBox([
            h("1. Spectra"),
            w.HBox([self.upload, self.folder, self.load_folder]),
            w.HTML(colab_hint()),
            w.HTML("<small>ORTEC .Spe/.Chn, N42, ASCII with a header; Canberra .cnf, ORTEC .spc and "
                   "IEC need the optional package becquerel. Every spectrum is calibrated on its own "
                   "natural lines.</small>"),
            self.live_time,
            h("2. Roles and masses"),
            self.sample, self.background, *self.ref.values(),
            w.HBox([self.sample_mass, self.ref_mass]),
            *self.ref_content.values(),
            self.water, self.equilibrium,
            self.analyse_btn, self.status,
            h("3. Result"),
            self.summary, self.out,
            h("4. Export"),
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

    # ---- files -------------------------------------------------------------
    def _on_upload(self, change) -> None:
        files = uploaded(self.upload)
        if not files:
            return
        paths = [self.pool.add_bytes(n, c) for n, c in files]
        clear_upload(self.upload)
        self._guard(lambda: self.load(paths))

    def _load_folder(self) -> None:
        self.load(self.pool.add_folder(self.folder.value.strip(), SUFFIXES))

    def load(self, paths) -> list[str]:
        """Add files and guess their roles from their names."""
        for p in paths:
            self.pool.add_path(p)
        names = self.pool.names()
        for d in (self.sample, self.background, *self.ref.values()):
            keep = d.value
            d.options = [NONE, *names]
            d.value = keep if keep in d.options else NONE
        for n in names:
            role = guess_role(n)
            if role == "background" and self.background.value == NONE:
                self.background.value = n
            elif role in self.ref and self.ref[role].value == NONE:
                self.ref[role].value = n
        taken = {self.background.value, *(d.value for d in self.ref.values())}
        if self.sample.value == NONE:
            rest = [n for n in names if n not in taken]
            if rest:
                self.sample.value = rest[0]
        self.status.value = message(f"{len(names)} files; check the roles below.", "ok")
        return names

    def set_roles(self, sample: str, U: str, Th: str, K: str, background: str | None = None) -> None:
        self.sample.value, self.ref["U"].value, self.ref["Th"].value, self.ref["K"].value = sample, U, Th, K
        self.background.value = background or NONE

    # ---- analysis ----------------------------------------------------------
    def _read(self, name: str):
        kw = {"live_time": float(self.live_time.value)} if self.live_time.value > 0 else {}
        sp = read_gamma(self.pool.paths[name], **kw)
        auto_calibrate(sp, NATURAL_LINES)
        return sp

    def analyse(self):
        if self.sample.value == NONE:
            raise ValueError("choose the sample spectrum")
        refs = {}
        for el, d in self.ref.items():
            if d.value == NONE:
                continue
            v, s = read_value(self.ref_content[el])
            refs[el] = (self._read(d.value), Reference(d.description, {el: (v, s)}, float(self.ref_mass.value)))
        if not refs:
            raise ValueError("choose at least one reference spectrum")
        bg = None if self.background.value == NONE else self._read(self.background.value)
        sample = self._read(self.sample.value)
        self.result = analyse(sample, float(self.sample_mass.value), refs, background=bg)
        water = read_value(self.water)
        self.sediment = self.result.sediment(water=water, equilibrium=bool(self.equilibrium.value))
        self._show(sample)
        if self.on_sediment is not None:
            self.on_sediment(self.sediment)
        return self.result

    def _show(self, sample) -> None:
        from .. import plot

        r, sed = self.result, self.sediment
        rows = []
        ura = None if sed.U_ra226 is None else as_value(sed.U_ra226).value
        dry = matrix_dose_rates(as_value(sed.U).value, as_value(sed.Th).value, as_value(sed.K).value, U_ra226=ura)
        w_now = read_value(self.water)[0]
        for rad in ("alpha", "beta", "gamma"):
            rows.append({"radiation": rad, "dry": f"{dry[rad]:.3f}",
                         "wet": f"{water_correction(dry[rad], w_now, rad):.3f}"})
        self.summary.value = (f"<pre style='font-size:12px'>{html.escape(r.summary())}</pre>"
                              "<b>Infinite-matrix dose rates (Gy/ka)</b>"
                              + table_html(rows, [("radiation", "radiation"), ("dry", "dry"),
                                                  ("wet", f"water {w_now:.0%}")]))
        import matplotlib.pyplot as plt

        fig, a1 = plt.subplots(figsize=(11, 3.8))
        plot.plot_gamma_spectrum(sample, ax=a1)
        fig.tight_layout()
        lines = plot.plot_gamma_lines(r)[0].figure
        show_figure(self.out, fig, lines)
        self.status.value = message("Done; the sediment is now in the Age tab.", "ok")

    def settings(self) -> dict:
        return {
            "eprdating": __version__,
            "sample": self.sample.value, "background": self.background.value,
            "references": {el: d.value for el, d in self.ref.items()},
            "reference_contents": {el: read_value(b) for el, b in self.ref_content.items()},
            "sample_mass_g": self.sample_mass.value, "reference_mass_g": self.ref_mass.value,
            "water": read_value(self.water), "equilibrium": self.equilibrium.value,
            "live_time_s": self.live_time.value,
        }

    def _export(self) -> None:
        if self.result is None:
            raise ValueError("analyse first")
        self.export_html.value = (download_link("eprdating_sediment.txt", self.result.summary(), "result (text)")
                                  + download_link("eprdating_sediment_settings.json",
                                                  json.dumps(self.settings(), indent=2), "settings (JSON)"))
