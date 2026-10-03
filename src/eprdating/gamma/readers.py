"""Readers for gamma-ray spectrum files.

:func:`read_gamma` picks the reader from the extension (and content):

==============  ===============================================  =========
extension       format                                           reader
==============  ===============================================  =========
``.Spe``        ORTEC / IAEA ASCII (``$DATA:``, ``$MEAS_TIM:``)  native
``.Chn``        ORTEC binary (MAESTRO)                           native
``.n42/.xml``   ANSI N42.42 (2006 and 2012 schemas)              native
``.txt``        ``# key: value`` header + channel columns        native
``.csv/.txt``   plain columns (counts, or channel/energy+counts) native
``.cnf``        Canberra Genie 2000                              becquerel
``.iec/.spc``   IEC 61455 / ORTEC binary SPC                     becquerel
==============  ===============================================  =========

The formats marked *becquerel* are read through the optional package
`becquerel <https://github.com/lbl-anp/becquerel>`_ (``pip install eprdating[gamma-formats]``),
whose parsers for these binary formats have been tested on real files.

Every reader returns a :class:`~eprdating.gamma.GammaSpectrum` with counts,
live and real time and, when the file has one, the stored energy
calibration (used only as a first guess; see
:func:`~eprdating.gamma.auto_calibrate`). Plain column files carry no times:
pass ``live_time`` (s).
"""

from __future__ import annotations

import re
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from .spectrum import Calibration, GammaSpectrum, read_spectrum_txt

# ---------------------------------------------------------------- ORTEC/IAEA SPE


def read_spe(path: str | Path) -> GammaSpectrum:
    """ORTEC / IAEA ``.Spe`` ASCII spectrum."""
    path = Path(path)
    lines = [ln.strip() for ln in path.read_text(encoding="latin-1").splitlines()]
    blocks: dict[str, list[str]] = {}
    key = None
    for ln in lines:
        if ln.startswith("$"):
            key = ln[1:].rstrip(":").strip()
            blocks[key] = []
        elif key is not None:
            blocks[key].append(ln)
    if "DATA" not in blocks:
        raise ValueError(f"{path}: no $DATA: block")
    first, last = (int(x) for x in blocks["DATA"][0].split()[:2])
    vals = " ".join(blocks["DATA"][1:]).split()
    counts = np.array([float(v) for v in vals[: last - first + 1]])
    live = real = None
    if blocks.get("MEAS_TIM"):
        t = blocks["MEAS_TIM"][0].split()
        live, real = float(t[0]), float(t[1]) if len(t) > 1 else None
    if live is None:
        raise ValueError(f"{path}: no live time ($MEAS_TIM:)")
    cal = None
    if "MCA_CAL" in blocks and len(blocks["MCA_CAL"]) >= 2:
        n = int(blocks["MCA_CAL"][0].split()[0])
        coef = [float(x) for x in blocks["MCA_CAL"][1].split()[:n]]
        if len(coef) > 1 and coef[1] != 0:
            cal = Calibration(tuple(coef))
    elif blocks.get("ENER_FIT"):
        coef = [float(x) for x in blocks["ENER_FIT"][0].split()]
        if len(coef) > 1 and coef[1] != 0:
            cal = Calibration(tuple(coef))
    start = blocks.get("DATE_MEA", [None])[0]
    header = {k: " / ".join(v) for k, v in blocks.items() if k != "DATA"}
    return GammaSpectrum(counts, live, real, start, path.stem, np.arange(first, last + 1, dtype=float), cal, header)


# --------------------------------------------------------------------- ORTEC CHN

_MONTHS = {m: i + 1 for i, m in enumerate(
    ("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"))}


def read_chn(path: str | Path) -> GammaSpectrum:
    """ORTEC ``.Chn`` binary spectrum (MAESTRO).

    Layout (little endian): int16 -1, int16 MCA number, int16 segment,
    char[2] start seconds, int32 real time and int32 live time in 20 ms
    ticks, char[8] date ``DDMMMYY*`` (``*`` = '1' after 2000), char[4] time
    ``HHMM``, int16 first channel, int16 number of channels, then uint32
    counts; an optional trailer (int16 -101 or -102) holds the energy
    calibration as float32 (offset, slope[, quadratic]).
    """
    path = Path(path)
    b = path.read_bytes()
    if len(b) < 32 or struct.unpack_from("<h", b, 0)[0] != -1:
        raise ValueError(f"{path}: not an ORTEC .Chn file")
    _, _mca, _seg, sec, rt, lt, date, hhmm, first, n = struct.unpack_from("<hhh2sii8s4shh", b, 0)
    counts = np.frombuffer(b, dtype="<u4", count=n, offset=32).astype(float)
    start = None
    try:
        d = date.decode("ascii")
        yy = int(d[5:7]) + (2000 if d[7] == "1" else 1900)
        t = hhmm.decode("ascii")
        start = f"{yy:04d}-{_MONTHS[d[2:5].upper()]:02d}-{int(d[:2]):02d} {t[:2]}:{t[2:]}:{sec.decode('ascii')}"
    except (UnicodeDecodeError, ValueError, KeyError):
        pass
    cal = None
    tr = 32 + 4 * n
    if len(b) >= tr + 16 and struct.unpack_from("<h", b, tr)[0] in (-101, -102):
        nc = 3 if struct.unpack_from("<h", b, tr)[0] == -102 else 2
        coef = struct.unpack_from(f"<{nc}f", b, tr + 4)
        if coef[1] != 0:
            cal = Calibration(tuple(float(c) for c in coef))
    return GammaSpectrum(counts, lt * 0.02, rt * 0.02, start, path.stem,
                         np.arange(first, first + n, dtype=float), cal)


# ----------------------------------------------------------------------- N42


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _duration(text: str) -> float:
    """ISO 8601 duration (``PT86400.5S``, ``PT1H2M3S``) in seconds."""
    m = re.fullmatch(r"\s*P(?:(\d+(?:\.\d*)?)D)?T?(?:(\d+(?:\.\d*)?)H)?(?:(\d+(?:\.\d*)?)M)?(?:(\d+(?:\.\d*)?)S)?\s*",
                     text)
    if not m:
        return float(text)
    d, h, mi, s = (float(x) if x else 0.0 for x in m.groups())
    return ((d * 24 + h) * 60 + mi) * 60 + s


def _channel_data(el) -> np.ndarray:
    vals = np.array(el.text.split(), float)
    if el.attrib.get("compressionCode", el.attrib.get("Compression", "")).lower().startswith("countedzero"):
        out = []
        i = 0
        while i < vals.size:
            if vals[i] == 0 and i + 1 < vals.size:
                out.extend([0.0] * int(vals[i + 1]))
                i += 2
            else:
                out.append(vals[i])
                i += 1
        vals = np.array(out)
    return vals


def read_n42(path: str | Path, index: int = 0) -> GammaSpectrum:
    """ANSI N42.42 XML spectrum (2006 or 2012). ``index`` selects the
    spectrum when the file holds several."""
    path = Path(path)
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:  # tolerate undeclared namespace prefixes (seen in real files)
        text = path.read_text(encoding="utf-8", errors="replace")
        text = re.sub(r"<(/?)[A-Za-z_][\w.-]*:(?=[A-Za-z_])", r"<\1", text)
        root = ET.fromstring(text.encode())
    spectra = [el for el in root.iter() if _local(el.tag) == "Spectrum"]
    if not spectra:
        raise ValueError(f"{path}: no <Spectrum> element")
    sp = spectra[index]
    ch = next((e for e in sp.iter() if _local(e.tag) == "ChannelData"), None)
    if ch is None:
        raise ValueError(f"{path}: no <ChannelData>")
    counts = _channel_data(ch)

    def find_time(names, scopes):
        for scope in scopes:
            for e in scope.iter():
                if _local(e.tag) in names and e.text:
                    return _duration(e.text)
        return None

    parent = next((p for p in root.iter() if sp in list(p)), root)
    live = find_time(("LiveTimeDuration", "LiveTime"), (sp, parent))
    real = find_time(("RealTimeDuration", "RealTime"), (sp, parent))
    if live is None:
        raise ValueError(f"{path}: no live time")
    cal = None
    for e in root.iter():
        if _local(e.tag) in ("CoefficientValues", "Coefficients") and e.text:
            ok = True
            if _local(e.tag) == "Coefficients":  # 2006: check it is an energy calibration
                cal_el = next((p for p in root.iter() if _local(p.tag) == "Calibration" and e in list(p.iter())), None)
                ok = cal_el is None or cal_el.attrib.get("Type", "Energy").lower() == "energy"
            coef = [float(x) for x in e.text.split()]
            if ok and len(coef) > 1 and coef[1] != 0:
                cal = Calibration(tuple(coef))
                break
    start = next((e.text for e in root.iter() if _local(e.tag) in ("StartDateTime", "StartTime") and e.text), None)
    return GammaSpectrum(counts, live, real, start, path.stem, np.arange(counts.size, dtype=float), cal)


# ----------------------------------------------------------------------- columns


def read_columns(path: str | Path, live_time: float, real_time: float | None = None,
                 counts_column: int = -1, energy_column: int | None = None,
                 delimiter: str | None = None, skiprows: int = 0) -> GammaSpectrum:
    """Plain text/CSV spectrum: one column of counts, or several columns
    with the counts in ``counts_column``. If ``energy_column`` is given its
    values are fitted with a linear calibration (first guess only)."""
    path = Path(path)
    if delimiter is None and path.suffix.lower() == ".csv":
        delimiter = ","
    a = np.genfromtxt(path, delimiter=delimiter, skip_header=skiprows, comments="#", invalid_raise=False)
    a = np.atleast_2d(a.T).T if a.ndim == 1 else a
    a = a[np.isfinite(a).all(axis=1)]
    counts = a[:, counts_column]
    cal = None
    if energy_column is not None:
        n = np.arange(counts.size, dtype=float)
        coef = np.polynomial.polynomial.polyfit(n, a[:, energy_column], 1)
        cal = Calibration(tuple(coef))
    return GammaSpectrum(counts, live_time, real_time, None, path.stem, np.arange(counts.size, dtype=float), cal)


# ----------------------------------------------------------------- becquerel


def read_with_becquerel(path: str | Path) -> GammaSpectrum:
    """Read any format supported by becquerel (CNF, SPC, IEC 61455, SPE)."""
    try:
        import becquerel as bq
    except ModuleNotFoundError as exc:
        raise ImportError(
            f"reading {Path(path).suffix} files needs the optional package becquerel: pip install becquerel"
        ) from exc
    import contextlib
    import io

    with contextlib.redirect_stdout(io.StringIO()):  # becquerel prints progress
        s = bq.Spectrum.from_file(str(path))
    counts = np.asarray(s.counts_vals, float)
    cal = None
    c = getattr(s, "energy_cal", None)
    expr = str(getattr(c, "expression", ""))
    coef = getattr(c, "params", None) if expr.startswith("p[0] + p[1] * x") else None  # polynomial only
    if coef is not None and len(coef) > 1:
        cal = Calibration(tuple(float(x) for x in coef))
    return GammaSpectrum(counts, float(s.livetime), float(s.realtime) if s.realtime else None,
                         str(s.start_time) if s.start_time else None, Path(path).stem,
                         np.arange(counts.size, dtype=float), cal)


# ----------------------------------------------------------------- dispatch


def read_gamma(path: str | Path, **kw) -> GammaSpectrum:
    """Read a gamma spectrum, choosing the reader from the file."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".spe":
        return read_spe(path)
    if ext == ".chn":
        return read_chn(path)
    if ext in (".n42", ".xml"):
        return read_n42(path, **kw)
    if ext in (".cnf", ".iec", ".spc"):
        return read_with_becquerel(path)
    head = path.read_text(encoding="latin-1", errors="replace")[:4000]
    if "# Channel data" in head or "Live time (s)" in head:
        return read_spectrum_txt(path)
    if head.lstrip().startswith("$SPEC_ID") or "$DATA:" in head:
        return read_spe(path)
    if "<" in head[:200] and "Spectrum" in head:
        return read_n42(path, **kw)
    if "live_time" not in kw:
        raise ValueError(f"{path}: plain column file, pass live_time=... (s)")
    return read_columns(path, **kw)


__all__ = ["read_chn", "read_columns", "read_gamma", "read_n42", "read_spe", "read_with_becquerel"]
