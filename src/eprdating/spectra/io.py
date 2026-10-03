"""Reading cw-EPR spectra from disk.

:func:`read_epr` reads any supported file: Bruker BES3T and ESP/WinEPR
(:mod:`eprdating.spectra.bruker`), the ``.dat``/``.par`` pairs described
below, and plain field/signal columns. All return a :class:`Spectrum` with
the field in mT.

``.dat``/``.par`` format
------------------------
``.dat`` + ``.par`` ASCII pairs (one spectrum, possibly several scans). The
``.par`` file holds ``KEY : value`` lines::

    N : 512            points per scan
    CF : 3350.0        centre field set by the operator (G)
    CF_ : 3354.0       actual centre field (G)
    SW : 499.1453      sweep width (G)
    SW_ : 500.0        second sweep-width value (kept, not used)
    Nscans : 4
    Freq : 9.43        microwave frequency (GHz)
    TC, MA, OF, PH, RG, CT   time constant, modulation amplitude, offset,
                             phase, receiver gain, conversion time (raw units)

and the ``.dat`` file has five whitespace-separated columns: running index,
point index within the scan, field (G), signal, and a flag. Scans follow one
another; a line with ``NaN`` signal may separate them, and an interrupted
acquisition leaves a trailing incomplete scan, which is dropped.

The field column is CF ± SW/2, i.e. it is built on the *set* centre field.
:func:`read_dat` moves it by ``CF_ - CF`` so that the axis is the actual
field (``actual_field=False`` keeps the raw axis). For the M18 series this
correction (+0.4 mT) brings the CO2- signal to within 0.05 mT of the
simulated one (literature g-values, nominal 9.43 GHz). Fields are in mT, as
everywhere in :mod:`eprdating`.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass, field, replace
from itertools import pairwise
from pathlib import Path

import numpy as np

_POWER_RE = re.compile(r"(\d+(?:\.\d+)?)\s*mW", re.IGNORECASE)


@dataclass
class Spectrum:
    """A cw-EPR (first-derivative) spectrum with its individual scans.

    B       : field grid in mT (increasing).
    scans   : array (n_scans, n_points) of raw signal.
    params  : acquisition parameters as read from the file (raw units).
    freq_GHz, power_mW : microwave frequency and power (for the .dat/.par
              format the power comes from the file name, e.g.
              ``M18_3_19mW_4SCAN.dat``).
    gain    : receiver gain as a linear factor (None if unknown).
    mod_amp_mT, time_constant_ms : field modulation (peak to peak) and
              lock-in time constant, when the file gives them in known units.
    """

    B: np.ndarray
    scans: np.ndarray
    name: str = ""
    params: dict[str, object] = field(default_factory=dict)
    freq_GHz: float | None = None
    power_mW: float | None = None
    dropped_points: int = 0
    gain: float | None = None
    mod_amp_mT: float | None = None
    time_constant_ms: float | None = None

    @property
    def n_scans(self) -> int:
        return self.scans.shape[0]

    @property
    def y(self) -> np.ndarray:
        """Average of the scans."""
        return self.scans.mean(axis=0)

    def window(self, lo_mT: float, hi_mT: float) -> Spectrum:
        """Copy restricted to ``lo_mT <= B <= hi_mT``."""
        m = (self.B >= lo_mT) & (self.B <= hi_mT)
        return replace(self, B=self.B[m], scans=self.scans[:, m], params=dict(self.params))

    def __repr__(self) -> str:
        return (f"Spectrum({self.name!r}, {self.B.size} points {self.B[0]:.2f}–{self.B[-1]:.2f} mT, "
                f"{self.n_scans} scan(s), {self.freq_GHz} GHz, {self.power_mW} mW)")


def read_par(path: str | Path) -> dict[str, float]:
    """Parse a ``KEY : value`` parameter file."""
    out: dict[str, float] = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        if ":" not in line:
            continue
        k, v = (s.strip() for s in line.split(":", 1))
        try:
            out[k] = float(v)
        except ValueError:
            continue
    return out


def read_dat(path: str | Path, par: str | Path | None = None, power_mW: float | None = None,
             field_unit: str = "G", actual_field: bool = True) -> Spectrum:
    """Read a ``.dat`` spectrum (and its ``.par`` file if present).

    ``par`` defaults to the file with the same stem. ``power_mW`` defaults to
    the number before ``mW`` in the file name. ``field_unit`` is the unit of
    the field column ("G" or "mT"). With ``actual_field`` the axis is shifted
    by ``CF_ - CF`` (actual minus set centre field) when both are present.
    """
    path = Path(path)
    par = Path(par) if par is not None else path.with_suffix(".par")
    params = read_par(par) if par.exists() else {}

    a = np.genfromtxt(path, usecols=(1, 2, 3))
    a = a[np.isfinite(a).all(axis=1)]
    if a.size == 0:
        raise ValueError(f"{path}: no data")
    idx = a[:, 0].astype(int)
    starts = np.r_[0, np.where(np.diff(idx) <= 0)[0] + 1, len(a)]
    chunks = [a[s:e] for s, e in pairwise(starts)]
    n = int(params.get("N", max(len(c) for c in chunks)))
    full = [c for c in chunks if len(c) == n]
    dropped = sum(len(c) for c in chunks if len(c) != n)
    if not full:
        raise ValueError(f"{path}: no complete scan of {n} points")
    B = full[0][:, 1].copy()
    for c in full[1:]:
        if not np.allclose(c[:, 1], B):
            raise ValueError(f"{path}: scans have different field axes")
    scale = {"G": 0.1, "mT": 1.0}[field_unit]
    if actual_field and "CF" in params and "CF_" in params:
        B = B + (params["CF_"] - params["CF"])
    if power_mW is None:
        m = _POWER_RE.search(path.stem)
        power_mW = float(m.group(1)) if m else None
    return Spectrum(
        B=B * scale,
        scans=np.array([c[:, 2] for c in full]),
        name=path.stem,
        params=params,
        freq_GHz=params.get("Freq"),
        power_mW=power_mW,
        dropped_points=dropped,
        gain=params.get("RG"),
    )


def read_columns(path: str | Path, field_unit: str = "G", freq_GHz: float | None = None,
                  power_mW: float | None = None, field_column: int = 0, signal_columns=None,
                  delimiter: str | None = None) -> Spectrum:
    """Plain text/CSV spectrum: a field column and one or more signal
    columns (each taken as a scan). Lines starting with ``#`` are skipped."""
    path = Path(path)
    if delimiter is None and path.suffix.lower() == ".csv":
        delimiter = ","
    a = np.genfromtxt(path, delimiter=delimiter, comments="#", invalid_raise=False)
    a = a[np.isfinite(a).all(axis=1)]
    cols = [c for c in range(a.shape[1]) if c != field_column] if signal_columns is None else list(signal_columns)
    B = a[:, field_column] * {"G": 0.1, "mT": 1.0, "T": 1000.0}[field_unit]
    order = np.argsort(B)
    return Spectrum(B[order], a[:, cols].T[:, order], path.stem, {}, freq_GHz, power_mW)


def read_epr(path: str | Path, **kw) -> Spectrum:
    """Read a cw-EPR spectrum, choosing the reader from the file.

    ======================  =============================================
    ``.DSC`` / ``.DTA``     Bruker BES3T (Xepr)
    ``.par`` + ``.spc``     Bruker ESP / WinEPR
    ``.dat`` + ``.par``     ``KEY : value`` + five-column ASCII (see above)
    ``.txt`` / ``.csv``     field and signal columns
    ======================  =============================================
    """
    from .bruker import read_bes3t, read_esp

    path = Path(path)
    ext = path.suffix.lower()
    if ext in (".dsc", ".dta"):
        return read_bes3t(path, **kw)
    if ext == ".spc" or (ext == ".par" and _sibling(path, ".spc") is not None):
        return read_esp(path, **kw)
    if ext == ".dat" or (ext == ".par" and _sibling(path, ".dat") is not None):
        return read_dat(path.with_suffix(".dat") if ext == ".par" else path, **kw)
    return read_columns(path, **kw)


def _sibling(path: Path, ext: str) -> Path | None:
    for e in (ext, ext.upper()):
        p = path.with_suffix(e)
        if p.exists():
            return p
    return None


def read_series(paths: Iterable[str | Path], **kw) -> list[Spectrum]:
    """Read several spectra (e.g. a dose series), any supported format."""
    return [read_epr(p, **kw) for p in paths]


__all__ = ["Spectrum", "read_columns", "read_dat", "read_epr", "read_par", "read_series"]
