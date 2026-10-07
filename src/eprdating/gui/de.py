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
    table_html,
    uploaded,
)

#: files read as spectra; companions (.par of a .dat, .DTA of a .DSC) are found by the readers
MAIN_SUFFIXES = (".dat", ".dsc", ".spc", ".xml", ".txt", ".csv")
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
        self.comparison: list[dict] = []
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
        self.weights = w.Dropdown(
            options=[("measured errors (noise, repeats)", "errors"), ("1/I² (relative)", "1/I^2"),
                     ("equal", "none")],
            value="errors", description="Weights", style=STYLE, layout=w.Layout(width="330px"),
            tooltip="1/I²: every point with the same relative error, common practice in ESR dating")
        self.negative = checkbox("allow De < 0 (shows where the data extrapolate)")
        self.fit_btn = button("Fit", "line-chart", primary=True, width="120px")
        self.fit_btn.on_click(lambda _: self._guard(self.fit))
        self.fit_out = w.Output()
        self.fit_text = w.HTML()
        self.compare_btn = button("Compare methods", "balance-scale", width="170px")
        self.compare_btn.on_click(lambda _: self._guard(self.compare_methods))
        self.compare_html = w.HTML()
        self.compare_out = w.Output()
        self.export_btn = button("Prepare downloads", "download")
        self.export_btn.on_click(lambda _: self._guard(self._export))
        self.export_html = w.HTML()
        self.status = w.HTML()

        h = lambda t: w.HTML(f"<h4 style='margin:10px 0 4px'>{t}</h4>")
        self.widget = w.VBox([
            h("1. Spectra"),
            w.HBox([self.upload, self.folder, self.load_folder]),
            w.HTML(colab_hint()),
            w.HTML("<small>.dat + .par, Bruker .DSC + .DTA or .spc + .par, Freiberg MS5000 .xml, text or CSV. "
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
            w.HBox([self.model, self.max_dose, self.weights]), self.negative, self.fit_btn,
            self.fit_text, self.fit_out,
            w.HTML("<small>The same De with the four intensity methods (same files, window, ticked points "
                   "and model):</small>"),
            self.compare_btn, self.compare_html, self.compare_out,
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

    def _sweep(self, s) -> tuple[float, float, float]:
        """Sweep of a spectrum (first and last field, step), the same for
        sweeps that differ by a few points (instruments that sample in time)."""
        lo, hi, step = float(s.B[0]), float(s.B[-1]), float(np.median(np.diff(s.B)))
        for k in self._sweeps:
            tol = 0.02 * (k[1] - k[0])
            if abs(lo - k[0]) <= tol and abs(hi - k[1]) <= tol and abs(step / k[2] - 1) < 0.02:
                return k
        key = (round(lo, 3), round(hi, 3), round(step, 5))
        self._sweeps.append(key)
        return key

    def _measure(self, groups: dict, window: IntensityWindow, method: str) -> tuple:
        """Intensity of every aliquot with ``method``.

        Spectra recorded with another sweep than most (e.g. a wide survey
        sweep of the natural) are measured on their own. An aliquot measured
        both ways bridges the two: the ratio of its intensities puts the
        other-sweep aliquots on the scale of the main sweep.

        Returns ``(points, notes, kind, template, ref_power)``.
        """
        used = [s for g in groups.values() for s in g["spectra"]]
        self._sweeps: list[tuple[float, float, float]] = []
        ref = float(self.ref_power.value) or next((s.power_mW for s in used if s.power_mW), None)
        sweeps = [self._sweep(s) for s in used]
        main = max(set(sweeps), key=sweeps.count)
        main_spectra = [s for s in used if self._sweep(s) == main]
        template = empirical_template(main_spectra, window, int(self.n_strong.value), float(self.max_shift.value)) \
            if method == "template" else None
        shift = float(self.max_shift.value)

        def measure(spectra, mass):
            kw = {"mass_mg": mass} if mass > 0 else {}
            return combined_intensity(spectra, method, template, window, ref_power_mW=ref, max_shift=shift, **kw)

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
                           "sweep": key, "result": r, "note": ""})
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
            if not (k > 0 and sk < 0.5 * k):
                notes.append(f"With the {method} method that factor is poorly determined, so "
                             f"{p['aliquot']} carries little weight.")
                kind = "warn"
        return points, notes, kind, template, ref

    def compute(self) -> list[dict]:
        """Intensity of every aliquot (repeats averaged; see :meth:`_measure`
        for spectra recorded with different sweeps), then a fit."""
        groups = self._groups()
        window = self.window()
        points, notes, kind, template, ref = self._measure(groups, window, self.method.value)
        old = {p["aliquot"]: p["fit"].value for p in self.points}
        for p in points:
            p["fit"] = w.Checkbox(value=old.get(p["aliquot"], True), indent=False, layout=w.Layout(width="40px"))
        self.template = template
        self.points = points
        self.comparison = []  # from other intensities
        self.compare_html.value = ""
        self.compare_out.clear_output()
        self._render_points()
        self._plot_spectra(window, ref)
        freq = next(s.freq_GHz for g in groups.values() for s in g["spectra"])
        note = f"{len(points)} aliquots, window {window.describe(freq)}. " + " ".join(notes)
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
        self.drc = self._fit(use)
        excluded = [(p["dose"], p["result"].value, p["result"].sigma) for p in self.points if not p["fit"].value]
        from .. import plot

        ax = plot.plot_dose_response(self.drc, excluded=excluded or None)
        show_figure(self.fit_out, ax.figure)
        d = self.drc
        self.fit_text.value = (f"<b>De = {d.De:.4g} ± {d.De_sigma:.2g} Gy</b> &nbsp; ({d.model}, "
                               f"{len(use)} points, weights {self.weights.label}, {self._fit_quality(d)})")
        if self.on_de is not None and np.isfinite(d.De):
            self.on_de(float(d.De), float(d.De_sigma))
        return self.drc

    def _fit_quality(self, d) -> str:
        """χ²ν with measured errors; with 1/I² the relative scatter about the
        curve (√χ²ν, since the errors are |I|); nothing for equal weights."""
        if self.weights.value == "1/I^2" or (self.weights.value == "errors" and d.sigma is not None
                                               and np.allclose(d.sigma, np.abs(d.intensity))):
            return f"scatter {100 * np.sqrt(d.chi2_red):.2g} %"
        if self.weights.value == "none" or d.sigma is None:
            return "equal weights"
        return f"χ²ν = {d.chi2_red:.2f}"

    def _fit(self, use: list[dict]):
        if len(use) < 3:
            raise ValueError("at least three points are needed for the fit")
        D = np.array([p["dose"] for p in use])
        I = np.array([p["result"].value for p in use])
        S = np.array([p["result"].sigma for p in use])
        md = float(self.max_dose.value) or None
        kw = {"max_dose": md, "De_min": -np.inf if self.negative.value else 0.0}
        if self.weights.value == "errors" and np.all(np.isfinite(S) & (S > 0)):
            return fit_dose_response(D, I, self.model.value, sigma=S, **kw)
        # 1/I² (relative errors, scaled by the scatter) or equal weights; also
        # the fallback when some measured errors are missing
        weighting = "1/I^2" if self.weights.value == "1/I^2" else "none"
        return fit_dose_response(D, I, self.model.value, weighting=weighting, **kw)

    # ---- comparison of methods ---------------------------------------------
    def compare_methods(self) -> list[dict]:
        """De with every intensity method, with the same files, window, ticked
        points and model. The methods agree within their errors when the
        signal is strong and alone in the window; a systematic difference
        points to noise bias (peak-to-peak, T1-B2), baseline problems
        (double integral) or other signals in the window."""
        groups = self._groups()
        window = self.window()
        ticked = {p["aliquot"]: bool(p["fit"].value) for p in self.points}
        rows = []
        for method in METHODS:
            try:
                points, _notes, _kind, _t, _ref = self._measure(groups, window, method)
                use = [p for p in points if ticked.get(p["aliquot"], True)]
                d = self._fit(use)
                rows.append({"method": method, "De_Gy": float(d.De), "sigma_Gy": float(d.De_sigma),
                             "chi2_red": float(d.chi2_red), "points": len(use), "problem": ""})
            except Exception as e:  # noqa: BLE001 - one method failing does not stop the others
                rows.append({"method": method, "De_Gy": float("nan"), "sigma_Gy": float("nan"),
                             "chi2_red": float("nan"), "points": 0, "problem": f"{type(e).__name__}: {e}"})
        self.comparison = rows
        self._show_comparison(rows)
        return rows

    def _show_comparison(self, rows: list[dict]) -> None:
        import matplotlib.pyplot as plt

        ok = [r for r in rows if np.isfinite(r["De_Gy"])]
        shown = []
        for r in rows:
            current = " ◀" if r["method"] == self.method.value else ""
            shown.append({"method": r["method"] + current,
                          "De": "—" if not np.isfinite(r["De_Gy"]) else f"{r['De_Gy']:.4g} ± {r['sigma_Gy']:.2g}",
                          "chi2": "—" if not np.isfinite(r["chi2_red"]) else f"{r['chi2_red']:.2f}",
                          "n": r["points"] or "—", "problem": r["problem"]})
        text = ""
        if len(ok) > 1:
            De = np.array([r["De_Gy"] for r in ok])
            sig = np.array([r["sigma_Gy"] for r in ok])
            spread = float(np.ptp(De))
            typical = float(np.median(sig))
            verdict = ("the differences are within the errors" if spread <= 2 * typical else
                       "the methods disagree by more than their errors: look for noise bias, baseline "
                       "problems or other signals in the window")
            text = (f"<p>Range of De between methods: {spread:.3g} Gy; typical error {typical:.2g} Gy: "
                    f"{verdict}. The methods are not independent (same spectra), so agreement is "
                    "expected only within the errors.</p>")
        self.compare_html.value = (table_html(shown, [("method", "method"), ("De", "De (Gy)"),
                                                      ("chi2", "χ²ν"), ("n", "points"), ("problem", "")])
                                   + text)
        if not ok:
            self.compare_out.clear_output()
            return
        fig, ax = plt.subplots(figsize=(6, 0.5 + 0.55 * len(ok)))
        y = np.arange(len(ok))[::-1]
        ax.errorbar([r["De_Gy"] for r in ok], y, xerr=[r["sigma_Gy"] for r in ok], fmt="o", color="#1f6feb",
                    ecolor="#57606a", elinewidth=1.5, capsize=3)
        ax.set_yticks(y, [r["method"] for r in ok])
        ax.axvline(0, color="#d0d7de", lw=1, zorder=0)
        ax.set_xlabel("De (Gy)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylim(-0.6, len(ok) - 0.4)
        fig.tight_layout()
        show_figure(self.compare_out, fig)

    def comparison_csv(self) -> str:
        if not self.comparison:
            return ""
        buf = io.StringIO()
        wr = csv.DictWriter(buf, fieldnames=list(self.comparison[0]))
        wr.writeheader()
        wr.writerows(self.comparison)
        return buf.getvalue()

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
            "fit": {"model": self.model.value, "max_dose_Gy": self.max_dose.value, "weights": self.weights.value,
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
        html = (download_link("eprdating_intensities.csv", self.results_csv(), "intensities (CSV)")
                + download_link("eprdating_De_settings.json", json.dumps(self.settings(), indent=2),
                                "settings (JSON)"))
        if self.comparison:
            html += download_link("eprdating_De_methods.csv", self.comparison_csv(), "methods (CSV)")
        self.export_html.value = html
