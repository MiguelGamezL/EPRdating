"""Bruker cw-EPR files.

* BES3T (Xepr, E500/E580/EMXplus/EMXmicro…): ``.DSC`` descriptor + ``.DTA``
  binary data.
* ESP / WinEPR (ESP300, EMX with WinEPR, ECS106…): ``.par`` parameters +
  ``.spc`` binary data.

The format rules follow EasySpin's ``eprload`` (Stoll & Schweiger, MIT
licence), against whose test files these readers are checked. Only 1-D
field sweeps and 2-D data whose second dimension is kept as separate rows
(e.g. repeated scans) are supported. The data are returned as stored: no
scaling by scans, gain or power is applied, but the parameters needed for it
are put into the :class:`~eprdating.spectra.Spectrum` fields.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .io import Spectrum

_FIELD_SCALE = {"G": 0.1, "MT": 1.0, "T": 1000.0}


def _num(text) -> float | None:
    if text is None:
        return None
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(text))
    return float(m.group()) if m else None


def _with_unit(text) -> tuple[float | None, str]:
    if text is None:
        return None, ""
    m = re.match(r"\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*([A-Za-z]*)", str(text))
    return (float(m.group(1)), m.group(2)) if m else (None, "")


def _companion(path: Path, ext: str) -> Path:
    for e in (ext.lower(), ext.upper()):
        p = path.with_suffix(e)
        if p.exists():
            return p
    raise FileNotFoundError(f"{path.with_suffix(ext)} not found")


# ---------------------------------------------------------------- BES3T


def read_dsc(path: str | Path) -> dict[str, str]:
    """Key/value pairs of a BES3T descriptor (comments and device blocks
    skipped, continuation lines joined, quotes removed)."""
    params: dict[str, str] = {}
    lines = Path(path).read_text(encoding="latin-1", errors="replace").splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i].rstrip()
        i += 1
        if not ln or ln[0] in "*#." or ln.startswith("\t*"):
            continue
        while ln.endswith("\\") and i < len(lines):  # continuation
            ln = ln[:-1] + "\n" + lines[i].rstrip()
            i += 1
        parts = ln.split(None, 1)
        key = parts[0]
        val = parts[1].strip() if len(parts) > 1 else ""
        if len(val) >= 2 and val[0] == val[-1] == "'":
            val = val[1:-1]
        params.setdefault(key, val)
    return params


_BES3T_FMT = {"C": "i1", "S": "i2", "I": "i4", "F": "f4", "D": "f8"}


def read_bes3t(path: str | Path, freq_GHz: float | None = None) -> Spectrum:
    """Bruker BES3T ``.DSC``/``.DTA`` cw spectrum."""
    path = Path(path)
    dsc = _companion(path, ".dsc")
    dta = _companion(path, ".dta")
    p = read_dsc(dsc)
    nx = int(_num(p["XPTS"]))
    ny = int(_num(p.get("YPTS", 1)) or 1)
    order = ">" if p.get("BSEQ", "BIG").upper() == "BIG" else "<"
    fmt = p.get("IRFMT", "D").split(",")[0].strip().upper()
    if fmt not in _BES3T_FMT:
        raise ValueError(f"{dsc}: unsupported IRFMT {fmt!r}")
    cplx = p.get("IKKF", "REAL").split(",")[0].strip().upper() == "CPLX"
    data = np.fromfile(dta, dtype=order + _BES3T_FMT[fmt]).astype(float)
    if cplx:
        data = data[0::2]  # real part (in-phase)
    if data.size < nx * ny:
        raise ValueError(f"{dta}: expected {nx * ny} values, found {data.size}")
    rows = data[: nx * ny].reshape(ny, nx)
    if p.get("XTYP", "IDX") == "IGD":
        xfmt = _BES3T_FMT[p.get("XFMT", "D").upper()]
        x = np.fromfile(_companion(path, ".xgf"), dtype=order + xfmt, count=nx).astype(float)
    else:
        x = _num(p["XMIN"]) + np.linspace(0.0, _num(p["XWID"]), nx)
    unit = p.get("XUNI", "G").strip("'").upper()
    if unit not in _FIELD_SCALE:
        raise ValueError(f"{dsc}: x axis in {unit!r} is not a magnetic field (only field sweeps are supported)")
    B = x * _FIELD_SCALE[unit]
    mwfq = _num(p.get("MWFQ"))
    mwpw = _num(p.get("MWPW"))
    rcag = _num(p.get("RCAG"))
    b0ma = _num(p.get("B0MA"))  # T
    rctc = _num(p.get("RCTC"))  # s
    if b0ma is None:  # device layer, e.g. "ModAmp  1.000 G"
        v, u = _with_unit(p.get("ModAmp"))
        b0ma = None if v is None else v * _FIELD_SCALE.get(u.upper(), 0.1) * 1e-3
    if rctc is None:
        v, u = _with_unit(p.get("TimeConst"))
        rctc = None if v is None else v * {"ms": 1e-3, "s": 1.0, "us": 1e-6}.get(u, 1e-3)
    return Spectrum(
        B=B, scans=rows, name=path.stem, params=p,
        freq_GHz=freq_GHz if freq_GHz is not None else (mwfq / 1e9 if mwfq else None),
        power_mW=mwpw * 1e3 if mwpw is not None else None,
        gain=10 ** (rcag / 20) if rcag is not None else None,
        mod_amp_mT=b0ma * 1e3 if b0ma is not None else None,
        time_constant_ms=rctc * 1e3 if rctc is not None else None,
    )


# ---------------------------------------------------------------- ESP / WinEPR


def read_esp_par(path: str | Path) -> dict[str, str]:
    """``KEY value`` pairs of an ESP/WinEPR ``.par`` file."""
    out: dict[str, str] = {}
    for ln in Path(path).read_text(encoding="latin-1", errors="replace").splitlines():
        parts = ln.strip().split(None, 1)
        if parts:
            out.setdefault(parts[0], parts[1].strip() if len(parts) > 1 else "")
    return out


def read_esp(path: str | Path, freq_GHz: float | None = None) -> Spectrum:
    """Bruker ESP / WinEPR ``.par``/``.spc`` spectrum.

    WinEPR files (``DOS`` key present) hold little-endian float32; ESP files
    big-endian int32 (EasySpin's rules).
    """
    path = Path(path)
    par = _companion(path, ".par")
    spc = _companion(path, ".spc")
    p = read_esp_par(par)
    winepr = "DOS" in p
    flags = int(_num(p.get("JSS", 0)) or 0)
    cplx = bool(flags & (1 << 4))
    two_d = bool(flags & (1 << 12))
    nx = 1024
    ny = 1
    if "SSX" in p and two_d:
        nx = int(_num(p["SSX"])) // (2 if cplx else 1)
    if "SSY" in p and two_d:
        ny = int(_num(p["SSY"]))
    if "ANZ" in p and not two_d:
        nx = int(_num(p["ANZ"])) // (2 if cplx else 1)
    if "RES" in p:
        nx = int(_num(p["RES"]))
    if "REY" in p and two_d:
        ny = int(_num(p["REY"]))
    if "XPLS" in p:
        nx = int(_num(p["XPLS"]))
    dtype = "<f4" if winepr else ">i4"
    data = np.fromfile(spc, dtype=dtype).astype(float)
    if cplx:
        data = data[0::2]
    if data.size < nx * ny:
        raise ValueError(f"{spc}: expected {nx * ny} values, found {data.size}")
    rows = data[: nx * ny].reshape(ny, nx)
    hcf, hsw, gst, gsi = (_num(p.get(k)) for k in ("HCF", "HSW", "GST", "GSI"))
    if p.get("JEX", "field-sweep") not in ("field-sweep", ""):
        raise ValueError(f"{par}: JEX = {p['JEX']!r}; only field sweeps are supported")
    if hcf is not None and hsw is not None and gst is not None and gsi is not None:
        x = gst + gsi * np.linspace(0, 1, nx)
    elif hcf is not None and hsw is not None:
        x = hcf + hsw / 2 * np.linspace(-1, 1, nx)
    elif gst is not None and gsi is not None:
        x = gst + gsi * np.linspace(0, 1, nx)
    else:
        raise ValueError(f"{par}: no field axis (HCF/HSW or GST/GSI)")
    mf = _num(p.get("MF"))
    mp = _num(p.get("MP"))
    rrg = _num(p.get("RRG"))
    rma = _num(p.get("RMA"))  # G
    rtc = _num(p.get("RTC"))  # ms
    return Spectrum(
        B=x * 0.1, scans=rows, name=path.stem, params=p,
        freq_GHz=freq_GHz if freq_GHz is not None else mf,
        power_mW=mp, gain=rrg, mod_amp_mT=rma * 0.1 if rma is not None else None,
        time_constant_ms=rtc,
    )


__all__ = ["read_bes3t", "read_dsc", "read_esp", "read_esp_par"]
