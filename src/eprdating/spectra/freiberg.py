"""Freiberg Instruments MS5000 spectra (ESRStudio ``.xml`` and ``.csv``).

The ``.xml`` file holds everything: the measurement attributes (microwave
frequency ``MwFreq`` in GHz, Q factor, temperature, time stamp), the recipe
(sweep, ``MicrowavePower`` in mW, ``Modulation`` in mT, ``Accumulations``)
and the curves as Base64 little-endian doubles, each sampled in time
(``XOffset + i * XSlope`` s): the field ``BField`` (mT) and the signal
``MWAbsorption``. The spectrum is the signal against the field interpolated
at the signal's sampling times; this reproduces the ``.csv`` export of
ESRStudio (which adds a constant offset to the signal and has no frequency).

The ``.csv`` export has ``key;value;description`` recipe lines and, after a
``Meas`` line, ``BField [mT];MW_Absorption []`` columns. When an ``.xml``
with the same name sits next to it, the frequency is taken from there.

A measurement repeated with a recipe of several runs is saved as
``Name_1.xml`` … ``Name_10.xml`` (each the average of ``Accumulations``
sweeps) and, often, their average ``Name_result.xml``. :func:`read_ms5000`
on a folder reads the runs (not the ``_result``) as the scans of one
spectrum.
"""

from __future__ import annotations

import base64
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from .io import Spectrum

_RUN_RE = re.compile(r"_(\d+)$")


def _floats(text: str) -> np.ndarray:
    """Base64 doubles; ESRStudio writes one padded chunk per value."""
    text = (text or "").strip()
    if not text:
        return np.empty(0)
    raw = b"".join(base64.b64decode(chunk + "=" * (-len(chunk) % 4)) for chunk in text.split("=") if chunk)
    return np.frombuffer(raw, "<f8").astype(float)


def _number(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


def is_ms5000(path: str | Path) -> bool:
    """Whether ``path`` is an ESRStudio ``.xml`` or ``.csv`` file."""
    path = Path(path)
    try:
        head = path.read_bytes()[:400].decode("utf-8-sig", errors="replace")
    except OSError:
        return False
    if path.suffix.lower() == ".xml":
        return "<ESRXmlFile" in head
    return head.lstrip().startswith("Name;") and "\nRecipe" in head.replace("\r", "")


def _from_xml(path: Path) -> tuple[list[tuple[np.ndarray, np.ndarray]], dict, str]:
    root = ET.parse(path).getroot()
    measurements = root.findall("./Data/Measurement")
    if not measurements:
        raise ValueError(f"{path}: no <Measurement> in the file")
    params: dict[str, object] = {}
    curves = []
    for meas in measurements:
        if not params:
            params = {k: _number(v) for k, v in meas.attrib.items()}
            for p in meas.iterfind("./Recipe/Parameters/Param"):
                params[p.get("Name")] = _number(p.text)
            recipe = meas.find("Recipe")
            if recipe is not None:
                params["Recipe"] = recipe.get("Name", "")
        data = {}
        for c in meas.iterfind("./DataCurves/Curve"):
            v = _floats(c.text)
            if v.size:
                t = float(c.get("XOffset", 0)) + float(c.get("XSlope", 1)) * np.arange(v.size)
                data[c.get("YType") or c.get("Name")] = (t, v)
        if "BField" not in data or "MW_Absorption" not in data:
            raise ValueError(f"{path}: the BField or MW_Absorption curve is missing")
        tB, B = data["BField"]
        tY, Y = data["MW_Absorption"]
        keep = (tY >= tB[0]) & (tY <= tB[-1])
        Bi = np.interp(tY[keep], tB, B)
        order = np.argsort(Bi)
        curves.append((Bi[order], Y[keep][order]))
    return curves, params, str(params.get("Name", path.stem))


def _from_csv(path: Path) -> tuple[list[tuple[np.ndarray, np.ndarray]], dict, str]:
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    params: dict[str, object] = {}
    name = path.stem
    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == "Meas")
    except StopIteration:
        raise ValueError(f"{path}: no 'Meas' section; not an ESRStudio export") from None
    for line in lines[:start]:
        parts = line.split(";")
        if len(parts) >= 2 and parts[0]:
            if parts[0] == "Name":
                name = parts[1]
            params[parts[0]] = _number(parts[1])
    rows = []
    for line in lines[start + 2:]:
        parts = line.split(";")
        if len(parts) >= 2:
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"{path}: no data")
    a = np.array(rows)
    a = a[np.argsort(a[:, 0])]
    xml = path.with_suffix(".xml")
    if xml.exists() and "MwFreq" not in params:
        try:
            root = ET.parse(xml).getroot()
            meas = root.find("./Data/Measurement")
            if meas is not None:
                params.update({k: _number(v) for k, v in meas.attrib.items() if k != "Name"})
        except ET.ParseError:
            pass
    return [(a[:, 0], a[:, 1])], params, name


def _spectrum(curves, params: dict, name: str, freq_GHz: float | None) -> Spectrum:
    B = curves[0][0]
    scans = np.array([np.interp(B, b, y) for b, y in curves])
    f = params.get("MwFreq")
    power = params.get("MicrowavePower")
    mod = params.get("Modulation")
    return Spectrum(
        B=B, scans=scans, name=name, params=params,
        freq_GHz=freq_GHz if freq_GHz is not None else (float(f) if isinstance(f, float) and f > 0 else None),
        power_mW=float(power) if isinstance(power, float) else None,
        mod_amp_mT=float(mod) if isinstance(mod, float) else None,
    )


def read_ms5000(path: str | Path, freq_GHz: float | None = None, include_result: bool = False) -> Spectrum:
    """Read an MS5000 (ESRStudio) ``.xml`` or ``.csv`` spectrum, or a folder of runs.

    path      : an ``.xml`` file (preferred: it has the frequency), its
                ``.csv`` export, or a folder whose ``Name_1.xml`` …
                ``Name_n.xml`` runs become the scans of one spectrum (put on
                the field grid of the first). ``Name_result.xml`` files in
                the folder are skipped unless ``include_result``.
    freq_GHz  : overrides the frequency in the file (needed for a ``.csv``
                without its ``.xml``).

    The power (mW) and modulation amplitude (mT) come from the recipe; there
    is no receiver gain (``gain`` is None).
    """
    path = Path(path)
    if path.is_dir():
        files = sorted(p for p in path.glob("*.xml") if is_ms5000(p))
        if not files:
            raise ValueError(f"{path}: no ESRStudio .xml files")
        runs = sorted((p for p in files if _RUN_RE.search(p.stem)),
                      key=lambda p: int(_RUN_RE.search(p.stem).group(1)))
        results = [p for p in files if p.stem.endswith("_result")]
        others = [p for p in files if p not in runs and p not in results]
        files = (runs or others or results) + (results if include_result and runs else [])
        curves, params = [], {}
        for f in files:
            c, p, _ = _from_xml(f)
            curves += c
            params = params or p
        params["files"] = [f.name for f in files]
        return _spectrum(curves, params, path.name, freq_GHz)
    if path.suffix.lower() == ".xml":
        curves, params, name = _from_xml(path)
    else:
        curves, params, name = _from_csv(path)
    return _spectrum(curves, params, name, freq_GHz)


__all__ = ["is_ms5000", "read_ms5000"]
