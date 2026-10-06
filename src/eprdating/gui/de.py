"""Equivalent-dose panel: spectra → intensities → dose response → De."""

from __future__ import annotations

import csv
import dataclasses
import io
import json
import math
import re
import warnings
from collections.abc import Callable

import ipywidgets as w
import numpy as np

from .. import __version__
from ..doseresponse import fit_dose_response
from ..spectra import (
    IntensityWindow,
    combined_intensity,
    empirical_template,
    field_for_g,
    g_for_field,
    read_epr,
)
from ..spectra.measure import METHODS
from ._common import (
    STYLE,
    FilePool,
    button,
    checkbox,
    clear_upload,
    colab_hint,
    download_link,
    message,
    scale_of,
    show_figure,
    uploaded,
)

#: files read as spectra; companions (.par of a .dat, .DTA of a .DSC) are found by the readers
MAIN_SUFFIXES = (".dat", ".dsc", ".spc", ".txt", ".csv")
ALL_SUFFIXES = (*MAIN_SUFFIXES, ".par", ".dta")
MODELS = ("LIN", "SSE", "EXPLIN", "DSE")
_DOSE_RE = re.compile(r"(\d+(?:[.,]\d+)?)\s*gy", re.IGNORECASE)


def dose_from_name(name: str) -> str:
    """Added dose written in a file name ("…_200Gy…"), or "" if there is none."""
    m = _DOSE_RE.search(name)
    return m.group(1).replace(",", ".") if m else ""


class DePanel:
    """Interactive equivalent-dose analysis (one tab of :func:`eprdating.gui.app`).

    Every step is also a method, so the panel can be driven from code:
    :meth:`load`, :meth:`compute`, :meth:`fit`, :meth:`settings`.
    """

    def __init__(self, on_de: Callable[[float, float], None] | None = None) -> None:
        self.on_de = on_de
        self.pool = FilePool()
        self.spectra: dict = {}
        self.rows: dict[str, dict] = {}
        self.points: list[dict] = []
        self.template = None
        self.drc = None
        self._build()

    # ---- widgets ---------------------------------------------------------
    def _build(self) -> None:
        self.upload = w.FileUpload(accept=",".join(ALL_SUFFIXES + (".DSC", ".DTA")), multiple=True,
                                   description="Upload spectra", layout=w.Layout(width="180px"))
        self.upload.observe(self._on_upload, names="value")
        self.folder = w.Text(placeholder="folder, e.g. /content/drive/MyDrive/M18", layout=w.Layout(width="420px"))
        self.load_folder = button("Load folder", "folder-open", width="140px")
        self.load_folder.on_click(lambda _: self._guard(self._load_folder))
        self.files_box = w.VBox()
        self.group_btn = w.Button(description="Group repeats by dose", layout=w.Layout(width="200px"),
                                  tooltip="Give the same aliquot label to files with the same dose")
        self.group_btn.on_click(lambda _: self.group_by_dose())

        self.width = w.FloatText(value=100.0, description="Window width", style=STYLE, layout=w.Layout(width="190px"))
        self.unit = w.Dropdown(options=["G", "mT"], value="G", layout=w.Layout(width="70px"))
        self.center_mode = w.Dropdown(options=["g", "G", "mT"], value="g", description="centred on",
                                      style=STYLE, layout=w.Layout(width="150px"),
                                      tooltip="centre given as a g-value or as a field")
        self.center_mode.observe(self._convert_center, names="value")
        self.center_g = w.FloatText(value=2.0023, description="=", style={"description_width": "14px"},
                                    layout=w.Layout(width="190px"))
        self.method = w.Dropdown(options=list(METHODS), value="template", description="Intensity", style=STYLE)
        self.n_strong = w.BoundedIntText(value=3, min=1, max=50, description="template: strongest", style=STYLE,
                                         layout=w.Layout(width="190px"))
        self.ref_power = w.FloatText(value=0.0, description="Power to normalise (mW, 0 = first)", style=STYLE,
                                     layout=w.Layout(width="290px"))
        self.max_shift = w.FloatText(value=0.6, description="Shift search (mT)", style=STYLE,
                                     layout=w.Layout(width="190px"))
        self.keep_mT = w.FloatText(value=30.0, description="Sweep kept around the centre (± mT, 0 = all)",
                                   style=STYLE, layout=w.Layout(width="340px"))
        self.freq = w.FloatText(value=0.0, description="Frequency for files without one (GHz)", style=STYLE,
                                layout=w.Layout(width="320px"))
        self.method.observe(lambda ch: setattr(self.n_strong, "disabled", ch["new"] != "template"), names="value")
        self.compute_btn = button("Compute intensities", "play", primary=True, width="200px")
        self.compute_btn.on_click(lambda _: self._guard(self.compute))

        self.points_box = w.VBox()
        self.spectra_out = w.Output()
        self.model = w.Dropdown(options=MODELS, value="LIN", description="Model", style=STYLE,
                                layout=w.Layout(width="160px"))
        self.max_dose = w.FloatText(value=0.0, description="Max dose (Gy, 0 = all)", style=STYLE,
                                    layout=w.Layout(width="220px"))
        self.negative = checkbox("allow De < 0 (shows where the data extrapolate)")
        self.fit_btn = button("Fit", "line-chart", primary=True, width="120px")
        self.fit_btn.on_click(lambda _: self._guard(self.fit))
        self.fit_out = w.Output()
        self.fit_text = w.HTML()
        self.export_btn = button("Prepare downloads", "download")
        self.export_btn.on_click(lambda _: self._guard(self._export))
        self.export_html = w.HTML()
        self.status = w.HTML()

        h = lambda t: w.HTML(f"<h4 style='margin:10px 0 4px'>{t}</h4>")
        self.widget = w.VBox([
            h("1. Spectra"),
            w.HBox([self.upload, self.folder, self.load_folder]),
            w.HTML(colab_hint()),
            w.HTML("<small>.dat + .par, Bruker .DSC + .DTA or .spc + .par, text or CSV. "
                   "Give the added dose of each file; files of the same aliquot get the same label "
                   "and are averaged.</small>"),
            self.files_box, self.group_btn,
            h("2. Intensities"),
            w.HBox([self.width, self.unit, self.center_mode, self.center_g]),
            w.HBox([self.method, self.n_strong]),
            w.HBox([self.ref_power, self.max_shift]),
            w.HBox([self.freq, self.keep_mT]),
            self.compute_btn, self.status, self.spectra_out,
            h("3. Dose response"),
            self.points_box,
            w.HBox([self.model, self.max_dose]), self.negative, self.fit_btn,
            self.fit_text, self.fit_out,
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
        except Exception as e:  # noqa: BLE001 - show the problem in the panel, not as a traceback
            self.status.value = message(f"{type(e).__name__}: {e}", "error")

    # ---- loading ---------------------------------------------------------
    def _on_upload(self, change) -> None:
        files = uploaded(self.upload)
        if not files:
            return
        paths = [self.pool.add_bytes(n, c) for n, c in files]
        clear_upload(self.upload)
        self._guard(lambda: self.load(paths))

    def _load_folder(self) -> None:
        paths = self.pool.add_folder(self.folder.value.strip(), ALL_SUFFIXES)
        self.load(paths)

    def load(self, paths) -> list[str]:
        """Read the spectra among ``paths`` (companion files are skipped) and
        add a row for each. Returns the names read."""
        read, problems = [], []
        for p in paths:
            p = self.pool.add_path(p)
            if p.suffix.lower() not in MAIN_SUFFIXES:
                continue  # companion files are read with their main file
            try:
                s = read_epr(p)
            except Exception as e:  # noqa: BLE001 - reported in the panel
                problems.append(f"{p.name}: {e}")
                continue
            self.spectra[p.name] = s
            if p.name not in self.rows:
                self.rows[p.name] = {
                    "dose": w.Text(value=dose_from_name(p.name), layout=w.Layout(width="80px")),
                    "aliquot": w.Text(value=p.stem, layout=w.Layout(width="200px")),
                    "mass": w.FloatText(value=0.0, layout=w.Layout(width="80px")),
                    "use": w.Checkbox(value=True, indent=False, layout=w.Layout(width="40px")),
                }
            read.append(p.name)
        self._render_files()
        text = f"{len(self.spectra)} spectra loaded."
        if problems:
            text += " Not read: " + "; ".join(problems)
        self.status.value = message(text, "warn" if problems else "ok")
        return read

    def _render_files(self) -> None:
        head = w.HBox([w.HTML(f"<b>{t}</b>", layout=w.Layout(width=wd)) for t, wd in
                       (("use", "40px"), ("file", "260px"), ("dose (Gy)", "80px"), ("aliquot", "200px"),
                        ("mass (mg, 0 = none)", "120px"), ("points · scans · GHz · mW", "220px"))])
        rows = [head]
        for name, s in self.spectra.items():
            r = self.rows[name]
            info = f"{s.B.size} · {s.n_scans} · {s.freq_GHz or '?'} · {s.power_mW or '?'}"
            rows.append(w.HBox([r["use"], w.Label(name, layout=w.Layout(width="260px")), r["dose"],
                                r["aliquot"], r["mass"], w.Label(info, layout=w.Layout(width="220px"))]))
        self.files_box.children = rows

    def set_file(self, name: str, dose: float | None = None, aliquot: str | None = None,
                 mass_mg: float | None = None, use: bool | None = None) -> None:
        """Fill a file's row from code."""
        r = self.rows[name]
        if dose is not None:
            r["dose"].value = f"{dose:g}"
        if aliquot is not None:
            r["aliquot"].value = aliquot
        if mass_mg is not None:
            r["mass"].value = float(mass_mg)
        if use is not None:
            r["use"].value = bool(use)

    def group_by_dose(self) -> None:
        for r in self.rows.values():
            if r["dose"].value.strip():
                r["aliquot"].value = f"{float(r['dose'].value):g} Gy"

    # ---- intensities -----------------------------------------------------
    def window(self) -> IntensityWindow:
        mode, x = self.center_mode.value, float(self.center_g.value)
        if mode == "g":
            return IntensityWindow(float(self.width.value), self.unit.value, center_g=x)
        return IntensityWindow(float(self.width.value), self.unit.value, center_mT=x / 10 if mode == "G" else x)

    def _frequency(self) -> float | None:
        f = next((s.freq_GHz for s in self.spectra.values() if s.freq_GHz), None)
        return f or (float(self.freq.value) or None)

    def _convert_center(self, change) -> None:
        """Keep the same centre when its unit changes (needs a frequency for g ↔ field)."""
        old, new, x = change["old"], change["new"], float(self.center_g.value)
        mT = {"G": lambda v: v / 10, "mT": lambda v: v}
        f = self._frequency()
        if old == "g":
            if f is None:
                return
            field = field_for_g(x, f)
        else:
            field = mT[old](x)
        if new == "g":
            if f is None:
                return
            value = round(float(g_for_field(field, f)), 5)
        else:
            value = round(field * 10 if new == "G" else field, 3)
        self.center_g.value = value

    def _groups(self) -> dict[str, dict]:
        groups: dict[str, dict] = {}
        for name, s in self.spectra.items():
            r = self.rows[name]
            if not r["use"].value:
                continue
            txt = r["dose"].value.strip().replace(",", ".")
            if not txt:
                raise ValueError(f"{name}: no dose given")
            dose = float(txt)
            g = groups.setdefault(r["aliquot"].value.strip() or name,
                                  {"dose": dose, "files": [], "spectra": [], "mass": 0.0})
            if not math.isclose(g["dose"], dose):
                raise ValueError(f"aliquot {r['aliquot'].value!r} has files with different doses")
            if s.freq_GHz is None:
                if not self.freq.value > 0:
                    raise ValueError(f"{name} has no microwave frequency; give it in 'Frequency for files "
                                     "without one'")
                s = dataclasses.replace(s, freq_GHz=float(self.freq.value))
            half = float(self.keep_mT.value)
            if half > 0:  # a wide survey sweep is cut to the region of the dating signal
                c = self.window().center(s.freq_GHz)
                if s.B.min() < c - half or s.B.max() > c + half:
                    s = s.window(c - half, c + half)
            g["files"].append(name)
            g["spectra"].append(s)
            g["mass"] = g["mass"] or float(r["mass"].value)
        if len(groups) < 2:
            raise ValueError("at least two aliquots are needed")
        return groups

    @staticmethod
    def _sweep(s) -> tuple[int, float]:
        return s.B.size, round(float(np.median(np.diff(s.B))), 4)

    def compute(self) -> list[dict]:
        """Intensity of every aliquot (repeats averaged), then a fit.

        Spectra recorded with another sweep than most (e.g. a wide survey
        sweep of the natural) are measured on their own. An aliquot measured
        both ways bridges the two: the ratio of its intensities puts the
        other-sweep aliquots on the scale of the main sweep.
        """
        groups = self._groups()
        window = self.window()
        used = [s for g in groups.values() for s in g["spectra"]]
        ref = float(self.ref_power.value) or next((s.power_mW for s in used if s.power_mW), None)
        method = self.method.value
        sweeps = [self._sweep(s) for s in used]
        main = max(set(sweeps), key=sweeps.count)
        main_spectra = [s for s in used if self._sweep(s) == main]
        template = empirical_template(main_spectra, window, int(self.n_strong.value), float(self.max_shift.value)) \
            if method == "template" else None
        shift = float(self.max_shift.value)

        def measure(spectra, mass):
            kw = {"mass_mg": mass} if mass > 0 else {}
            return combined_intensity(spectra, method, template, window, ref_power_mW=ref, max_shift=shift, **kw)

        old = {p["aliquot"]: p["fit"].value for p in self.points}
        points, bridges = [], {}
        for label, g in sorted(groups.items(), key=lambda kv: kv[1]["dose"]):
            by_sweep: dict = {}
            for name, sp in zip(g["files"], g["spectra"], strict=True):
                by_sweep.setdefault(self._sweep(sp), []).append((name, sp))
            if len(by_sweep) > 2 or (len(by_sweep) == 2 and main not in by_sweep):
                raise ValueError(f"aliquot {label!r} mixes sweeps that cannot be related")
            key = main if main in by_sweep else next(iter(by_sweep))
            files = [n for n, _ in by_sweep[key]]
            spectra = [sp for _, sp in by_sweep[key]]
            r = measure(spectra, g["mass"])
            for other, items in by_sweep.items():
                if other != key:  # measured both ways: a bridge between the sweeps
                    ro = measure([sp for _, sp in items], g["mass"])
                    k = r.value / ro.value
                    sk = abs(k) * float(np.hypot(r.sigma / r.value, ro.sigma / ro.value))
                    bridges.setdefault(other, []).append((k, sk, label))
            points.append({"aliquot": label, "dose": g["dose"], "files": files, "spectra": spectra,
                           "sweep": key, "result": r, "note": "",
                           "fit": w.Checkbox(value=old.get(label, True), indent=False,
                                             layout=w.Layout(width="40px"))})
        notes, kind = [], "ok"
        for p in points:
            if p["sweep"] == main:
                continue
            if p["sweep"] not in bridges:
                notes.append(f"{p['aliquot']} was measured with another sweep and no aliquot links the two; "
                             "its intensity may be on another scale.")
                kind = "warn"
                continue
            ks = np.array([b[0] for b in bridges[p["sweep"]]])
            sks = np.array([b[1] for b in bridges[p["sweep"]]])
            wts = 1 / sks**2
            k, sk = float(np.sum(wts * ks) / wts.sum()), float(1 / np.sqrt(wts.sum()))
            r = p["result"]
            value = r.value * k
            sigma = float(np.hypot(r.sigma * k, r.value * sk))
            p["result"] = dataclasses.replace(r, value=value, sigma=sigma)
            p["note"] = f"×{k:.3f}±{sk:.3f}"
            notes.append(f"{p['aliquot']} put on the main sweep with the factor {k:.3f} ± {sk:.3f} from "
                         f"{', '.join(b[2] for b in bridges[p['sweep']])}.")
        self.template = template
        self.points = points
        self._render_points()
        self._plot_spectra(window, ref)
        note = f"{len(points)} aliquots, window {window.describe(used[0].freq_GHz)}. " + " ".join(notes)
        self.status.value = message(note, kind)
        self.fit()
        return self.results()

    def results(self) -> list[dict]:
        out = []
        for p in self.points:
            r = p["result"]
            out.append({"aliquot": p["aliquot"], "dose_Gy": p["dose"], "files": " + ".join(p["files"]),
                        "scale": p["note"],
                        "intensity": r.value, "sigma": r.sigma, "p_noise": r.p_noise, "n_repeats": r.n_repeats,
                        "chi2_red": r.chi2_red, "shift_mT": r.shift_mT, "in_fit": bool(p["fit"].value)})
        return out

    def set_point(self, aliquot: str, in_fit: bool) -> None:
        next(p for p in self.points if p["aliquot"] == aliquot)["fit"].value = bool(in_fit)

    def _render_points(self) -> None:
        f, unit = scale_of([p["result"].value for p in self.points])
        cols = (("fit", "40px"), ("aliquot", "200px"), ("dose (Gy)", "80px"), (f"intensity{unit}", "170px"),
                ("p_noise", "80px"), ("repeats", "130px"))
        rows = [w.HBox([w.HTML(f"<b>{t}</b>", layout=w.Layout(width=wd)) for t, wd in cols])]
        for p in self.points:
            r = p["result"]
            flag = " ⚠" if r.p_noise is not None and r.p_noise >= 0.01 else ""
            rep = f"{r.n_repeats}, χ²ν {r.chi2_red:.2f}" if r.n_repeats > 1 else "1"
            if p["note"]:
                rep += f"  {p['note']}"
            cells = [p["fit"], w.Label(p["aliquot"]), w.Label(f"{p['dose']:g}"),
                     w.Label(f"{r.value * f:.4g} ± {r.sigma * f:.2g}"),
                     w.Label("—" if r.p_noise is None else f"{r.p_noise:.3f}{flag}"), w.Label(rep)]
            rows.append(w.HBox([c if i == 0 else w.Box([c], layout=w.Layout(width=cols[i][1]))
                                for i, c in enumerate(cells)]))
        rows.append(w.HTML("<small>Untick a point to leave it out of the fit. ⚠ p_noise ≥ 0.01: the "
                           "signal is not clearly above the noise; leaving such points out biases De "
                           "upwards, so keep them unless the measurement failed.</small>"))
        self.points_box.children = rows

    def _plot_spectra(self, window, ref) -> None:
        from .. import plot
        from ..spectra import combine_spectra

        spectra, fits, labels = [], [], []
        for p in self.points:
            c = combine_spectra(p["spectra"], window, ref_power_mW=ref).spectrum
            m = window.mask(c.B, c.freq_GHz)
            spectra.append((c.B[m], c.y[m]))
            r = p["result"]
            fits.append((c.B[m], r.fit.fitted) if r.fit is not None else None)
            labels.append(f"{p['dose']:g} Gy")
        fits = fits if all(f is not None for f in fits) else None
        ax = plot.plot_spectra(spectra, labels, fits=fits, title="Spectra in the intensity window")
        show_figure(self.spectra_out, ax.figure)

    # ---- fit ---------------------------------------------------------------
    def fit(self):
        """Dose-response fit of the ticked points."""
        if not self.points:
            raise ValueError("compute the intensities first")
        use = [p for p in self.points if p["fit"].value]
        if len(use) < 3:
            raise ValueError("at least three points are needed for the fit")
        D = np.array([p["dose"] for p in use])
        I = np.array([p["result"].value for p in use])
        S = np.array([p["result"].sigma for p in use])
        md = float(self.max_dose.value) or None
        self.drc = fit_dose_response(D, I, self.model.value, sigma=S if np.all(np.isfinite(S) & (S > 0)) else None,
                                     max_dose=md, De_min=-np.inf if self.negative.value else 0.0)
        excluded = [(p["dose"], p["result"].value, p["result"].sigma) for p in self.points if not p["fit"].value]
        from .. import plot

        ax = plot.plot_dose_response(self.drc, excluded=excluded or None)
        show_figure(self.fit_out, ax.figure)
        d = self.drc
        self.fit_text.value = (f"<b>De = {d.De:.4g} ± {d.De_sigma:.2g} Gy</b> &nbsp; ({d.model}, "
                               f"{len(use)} points, χ²ν = {d.chi2_red:.2f})")
        if self.on_de is not None and np.isfinite(d.De):
            self.on_de(float(d.De), float(d.De_sigma))
        return self.drc

    # ---- export ------------------------------------------------------------
    def settings(self) -> dict:
        """Everything needed to repeat the analysis."""
        return {
            "eprdating": __version__,
            "files": {n: {"dose_Gy": r["dose"].value, "aliquot": r["aliquot"].value,
                          "mass_mg": r["mass"].value, "use": r["use"].value} for n, r in self.rows.items()},
            "window": {"width": self.width.value, "unit": self.unit.value,
                       "center": self.center_g.value, "center_as": self.center_mode.value},
            "method": self.method.value, "template_strongest": self.n_strong.value,
            "ref_power_mW": self.ref_power.value, "max_shift_mT": self.max_shift.value,
            "frequency_fallback_GHz": self.freq.value, "sweep_kept_mT": self.keep_mT.value,
            "fit": {"model": self.model.value, "max_dose_Gy": self.max_dose.value,
                    "De_may_be_negative": self.negative.value,
                    "points": {p["aliquot"]: bool(p["fit"].value) for p in self.points}},
            "De_Gy": None if self.drc is None else [self.drc.De, self.drc.De_sigma],
        }

    def results_csv(self) -> str:
        buf = io.StringIO()
        rows = self.results()
        if not rows:
            return ""
        wr = csv.DictWriter(buf, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
        return buf.getvalue()

    def _export(self) -> None:
        self.export_html.value = (download_link("eprdating_intensities.csv", self.results_csv(), "intensities (CSV)")
                                  + download_link("eprdating_De_settings.json",
                                                  json.dumps(self.settings(), indent=2), "settings (JSON)"))
