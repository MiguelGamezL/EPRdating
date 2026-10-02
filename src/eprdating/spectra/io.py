"""Reading cw-EPR spectra from disk.

Supported format
----------------
``.dat`` + ``.par`` ASCII pairs (one spectrum, possibly several scans). The
``.par`` file holds ``KEY : value`` lines::

    N : 512            points per scan
    CF : 3350.0        actual centre field (G); the field axis is CF ± SW/2
    CF_ : 3354.0       centre field set by the operator (kept, not used)
    SW : 499.1453      actual sweep width (G)
    SW_ : 500.0        sweep width set by the operator
    Nscans : 4
    Freq : 9.43        microwave frequency (GHz)
    TC, MA, OF, PH, RG, CT   time constant, modulation amplitude, offset,
                             phase, receiver gain, conversion time (raw units)

and the ``.dat`` file has five whitespace-separated columns: running index,
point index within the scan, field (G), signal, and a flag. Scans follow one
another; a line with ``NaN`` signal may separate them, and an interrupted
acquisition leaves a trailing incomplete scan, which is dropped.

Fields are converted to mT, the unit used everywhere in :mod:`eprdating`.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass, field
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
    freq_GHz, power_mW : microwave frequency and power (power usually comes
              from the file name, e.g. ``M18_3_19mW_4SCAN.dat``).
    """

    B: np.ndarray
    scans: np.ndarray
    name: str = ""
    params: dict[str, float] = field(default_factory=dict)
    freq_GHz: float | None = None
    power_mW: float | None = None
    dropped_points: int = 0

    @property
    def n_scans(self) -> int:
        return self.scans.shape[0]

    @property
    def y(self) -> np.ndarray:
        """Average of the scans."""
        return self.scans.mean(axis=0)

    @property
    def gain(self) -> float | None:
        return self.params.get("RG")

    def window(self, lo_mT: float, hi_mT: float) -> Spectrum:
        """Copy restricted to ``lo_mT <= B <= hi_mT``."""
        m = (self.B >= lo_mT) & (self.B <= hi_mT)
        return Spectrum(self.B[m], self.scans[:, m], self.name, dict(self.params), self.freq_GHz,
                        self.power_mW, self.dropped_points)

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
             field_unit: str = "G") -> Spectrum:
    """Read a ``.dat`` spectrum (and its ``.par`` file if present).

    ``par`` defaults to the file with the same stem. ``power_mW`` defaults to
    the number before ``mW`` in the file name. ``field_unit`` is the unit of
    the field column ("G" or "mT").
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
    )


def read_series(paths: Iterable[str | Path], **kw) -> list[Spectrum]:
    """Read several spectra (e.g. a dose series)."""
    return [read_dat(p, **kw) for p in paths]


__all__ = ["Spectrum", "read_dat", "read_par", "read_series"]
