"""U, Th and K contents by the comparative (relative) method.

The sample and reference materials of known content are measured in the same
container and geometry. For every gamma line the background-corrected count
rate per unit mass of the sample is divided by that of the reference:

    C_sample = C_ref · (R_s / m_s) / (R_ref / m_ref),   R = net area / live time

Efficiency, emission probability and coincidence summing cancel; what does
not cancel is a difference in fill height or density (self-absorption),
which matters most below ~200 keV and when the masses differ.

Line groups
-----------
``K``      40K, 1460.8 keV.
``Ra226``  214Pb and 214Bi lines: 226Ra through its radon daughters. Gives
           the U content if the chain is in equilibrium (sample sealed for
           ~3–4 weeks so that 222Rn and its daughters grow back).
``U238``   234Th (63.3 keV) and 234mPa (1001.0 keV): the top of the chain,
           i.e. 238U itself. Weak; compared with ``Ra226`` it tests the
           equilibrium the dose-rate model assumes.
``Th232``  228Ac, 212Pb and 208Tl lines.

Peak areas come from Gaussians at the calibrated positions and widths of each
spectrum (calibrated separately, so gain drift between days is absorbed) on a
linear background; only amplitudes and background are fitted, which keeps
weak peaks unbiased. Lines of one group share the reference uncertainty, so
it is added after averaging.

Reference values (dry mass, 1σ)
-------------------------------
IAEA-RGU-1: 238U 4941 ± 99 Bq/kg → U = 400 ± 8 µg/g (certificate);
IAEA-RGTh-1: Th = 800 ± 16 µg/g (IAEA/RL/148);
IAEA-RGK-1: K = 448 ± 3 g/kg (certificate, 2016).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .._types import Value, ValueLike, as_value
from ..dose_rate import Sediment
from .spectrum import NATURAL_LINES, Calibration, GammaSpectrum


@dataclass(frozen=True)
class Line:
    energy: float  # keV
    emitter: str
    group: str
    neighbours: tuple[float, ...] = ()  # other peaks fitted in the same window
    window: float | None = None  # half width (keV); default 8σ beyond the outermost peak


LINES: tuple[Line, ...] = (
    Line(1460.820, "40K", "K"),
    Line(295.224, "214Pb", "Ra226", (300.087,)),
    Line(351.932, "214Pb", "Ra226"),
    Line(609.312, "214Bi", "Ra226"),
    Line(1120.287, "214Bi", "Ra226"),
    Line(1764.494, "214Bi", "Ra226"),
    Line(63.29, "234Th", "U238"),
    Line(1001.03, "234mPa", "U238"),
    Line(238.632, "212Pb", "Th232", (241.997,)),
    Line(338.320, "228Ac", "Th232"),
    Line(583.187, "208Tl", "Th232"),
    Line(911.204, "228Ac", "Th232"),
    Line(968.971, "228Ac", "Th232", (964.766,)),
    Line(2614.511, "208Tl", "Th232"),
)

#: which element each group measures
GROUP_ELEMENT = {"K": "K", "Ra226": "U", "U238": "U", "Th232": "Th"}


@dataclass
class Reference:
    """A reference material: content per element (U, Th in µg/g, K in %)."""

    name: str
    content: Mapping[str, ValueLike]
    mass_g: float = 500.0


IAEA_RGU_1 = Reference("IAEA-RGU-1", {"U": (400.0, 8.0)})
IAEA_RGTH_1 = Reference("IAEA-RGTh-1", {"Th": (800.0, 16.0)})
IAEA_RGK_1 = Reference("IAEA-RGK-1", {"K": (44.8, 0.3)})


def _gauss(x, mu, s):
    return np.exp(-0.5 * ((x - mu) / s) ** 2)


def line_area(spec: GammaSpectrum, line: Line, cal: Calibration | None = None) -> tuple[float, float]:
    """Area (counts) of ``line`` and its 1σ, with neighbours fitted jointly.

    Gaussian positions and widths come from the spectrum calibration; the
    amplitudes and a linear background are linear parameters, fitted with
    Poisson weights from the model (two iterations).
    """
    cal = cal or spec.calibration
    if cal is None:
        raise ValueError(f"{spec.name}: spectrum not calibrated")
    E = cal.energy(spec.channels)
    half = line.window
    if half is None:
        span = max((abs(n - line.energy) for n in line.neighbours), default=0.0)
        half = span + 8.0 * float(cal.sigma_keV(line.energy))
    m = np.abs(E - line.energy) <= half
    x, y = E[m], spec.counts[m]
    peaks = (line.energy, *line.neighbours)
    cols = [_gauss(x, e, cal.sigma_keV(e)) for e in peaks]
    A = np.column_stack([*cols, np.ones_like(x), x - line.energy])
    w = 1.0 / np.maximum(y, 1.0)
    for _ in range(3):
        Aw = A * np.sqrt(w)[:, None]
        p, *_ = np.linalg.lstsq(Aw, y * np.sqrt(w), rcond=None)
        w = 1.0 / np.maximum(A @ p, 1.0)
    cov = np.linalg.inv((A * w[:, None]).T @ A)
    gain = float(np.median(np.diff(E)))  # keV per channel
    k = cal.sigma_keV(line.energy) * np.sqrt(2 * np.pi) / gain
    return float(p[0] * k), float(np.sqrt(cov[0, 0]) * k)


@dataclass
class LineResult:
    line: Line
    sample_rate: float  # net counts/s
    sample_rate_sigma: float
    ref_rate: float
    ref_rate_sigma: float
    ratio: float  # (R_s/m_s)/(R_ref/m_ref)
    ratio_sigma: float
    content: float  # ratio × reference value
    content_sigma: float  # counting only


@dataclass
class GroupResult:
    group: str
    element: str
    content: Value  # including the reference uncertainty
    ratio: float
    ratio_sigma: float
    chi2_red: float
    lines: list[LineResult] = field(default_factory=list)


@dataclass
class GammaResult:
    sample: str
    mass_g: float
    groups: dict[str, GroupResult]
    activity_ratio_ra226_u238: Value | None = None

    @property
    def U(self) -> Value:
        """U (µg/g) from the 226Ra daughters (equilibrium assumed)."""
        return self.groups["Ra226"].content

    @property
    def Th(self) -> Value:
        return self.groups["Th232"].content

    @property
    def K(self) -> Value:
        return self.groups["K"].content

    def sediment(self, water: ValueLike = 0.0, equilibrium: bool = True) -> Sediment:
        """A :class:`~eprdating.Sediment` with these (dry-mass) contents.

        With ``equilibrium=False`` the U chain is split: 238U from the
        234Th/234mPa lines and 226Ra (+ daughters) from 214Pb/214Bi, which are
        independent measurements. Use it when the samples were sealed long
        enough for radon to grow back and the 226Ra/238U ratio departs from 1.
        """
        if equilibrium:
            return Sediment(U=self.U, Th=self.Th, K=self.K, water=water)
        if "U238" not in self.groups:
            raise ValueError("no 238U (234Th/234mPa) lines measured")
        return Sediment(U=self.groups["U238"].content, Th=self.Th, K=self.K, water=water, U_ra226=self.U)

    def summary(self) -> str:
        unit = {"K": "%", "U": "µg/g", "Th": "µg/g"}
        out = [f"{self.sample} ({self.mass_g:g} g)"]
        for g in self.groups.values():
            c = g.content
            out.append(f"  {g.group:6s} {g.element:2s} = {c.value:.4g} ± {c.sigma:.2g} {unit[g.element]}"
                       f"   ({len(g.lines)} lines, chi2_red = {g.chi2_red:.2g})")
            for lr in g.lines:
                out.append(f"      {lr.line.emitter:6s} {lr.line.energy:8.2f} keV  "
                           f"{lr.content:.4g} ± {lr.content_sigma:.2g}")
        if self.activity_ratio_ra226_u238 is not None:
            r = self.activity_ratio_ra226_u238
            out.append(f"  226Ra/238U activity ratio = {r.value:.3g} ± {r.sigma:.2g}")
        return "\n".join(out)


def net_rate(spec: GammaSpectrum, line: Line, background: GammaSpectrum | None) -> tuple[float, float]:
    a, sa = line_area(spec, line)
    r, v = a / spec.live_time, (sa / spec.live_time) ** 2
    if background is not None:
        b, sb = line_area(background, line)
        r -= b / background.live_time
        v += (sb / background.live_time) ** 2
    return r, float(np.sqrt(v))


def _weighted(x, s):
    x, s = np.asarray(x, float), np.asarray(s, float)
    w = 1 / s**2
    m = float(np.sum(w * x) / np.sum(w))
    e = float(np.sqrt(1 / np.sum(w)))
    chi = float(np.sum(w * (x - m) ** 2) / (x.size - 1)) if x.size > 1 else 0.0
    return m, e * np.sqrt(max(chi, 1.0)), chi


def analyse(
    sample: GammaSpectrum,
    mass_g: float,
    references: Mapping[str, tuple[GammaSpectrum, Reference]],
    background: GammaSpectrum | None = None,
    lines: Sequence[Line] = LINES,
    groups: Sequence[str] | None = None,
) -> GammaResult:
    """Contents of the sample by the comparative method.

    ``references`` maps an element ("U", "Th", "K") to (spectrum, Reference).
    All spectra must be calibrated (see :func:`calibrate_natural`). Lines of
    groups whose element has no reference are skipped.
    """
    by_group: dict[str, list[LineResult]] = {}
    for line in lines:
        el = GROUP_ELEMENT[line.group]
        if el not in references or (groups is not None and line.group not in groups):
            continue
        ref_spec, ref = references[el]
        rs, srs = net_rate(sample, line, background)
        rr, srr = net_rate(ref_spec, line, background)
        if rr <= 0:
            continue
        k = ref.mass_g / mass_g
        ratio = k * rs / rr
        sratio = k * np.hypot(srs / rr, rs * srr / rr**2)
        c = as_value(ref.content[el]).value
        by_group.setdefault(line.group, []).append(
            LineResult(line, rs, srs, rr, srr, ratio, float(sratio), ratio * c, float(sratio) * c))
    out = {}
    for g, lrs in by_group.items():
        el = GROUP_ELEMENT[g]
        m, e, chi = _weighted([lr.ratio for lr in lrs], [lr.ratio_sigma for lr in lrs])
        ref_val = as_value(references[el][1].content[el])
        c = m * ref_val.value
        sc = abs(c) * np.hypot(e / m if m else 0.0, ref_val.sigma / ref_val.value)
        out[g] = GroupResult(g, el, Value(c, float(sc)), m, e, chi, lrs)
    ar = None
    if "Ra226" in out and "U238" in out and out["U238"].ratio > 0:
        a, b = out["Ra226"], out["U238"]
        r = a.ratio / b.ratio  # same (equilibrium) reference -> activity ratio
        ar = Value(r, abs(r) * float(np.hypot(a.ratio_sigma / a.ratio, b.ratio_sigma / b.ratio)))
    return GammaResult(sample.name, mass_g, out, ar)


def analyse_files(
    sample: str | Path,
    mass_g: float,
    references: Mapping[str, tuple[str | Path, Reference]],
    background: str | Path | None = None,
    lines: Sequence[Line] = LINES,
    calibration_lines: Sequence[float] = NATURAL_LINES,
    **read_kw,
) -> GammaResult:
    """One call from files to contents.

    Every file is read with :func:`~eprdating.gamma.read_gamma` (any
    supported format) and calibrated on its own natural lines with
    :func:`~eprdating.gamma.auto_calibrate` (no first guess needed), then
    :func:`analyse` is applied. Example::

        r = analyse_files("soil.Spe", 600.0,
                          {"U": ("RGU1.Spe", IAEA_RGU_1), "Th": ("RGTh1.Spe", IAEA_RGTH_1),
                           "K": ("RGK1.Spe", IAEA_RGK_1)},
                          background="empty.Spe")
    """
    from .readers import read_gamma
    from .spectrum import auto_calibrate

    def load(p):
        sp = read_gamma(p, **read_kw)
        auto_calibrate(sp, calibration_lines)
        return sp

    refs = {el: (load(p), ref) for el, (p, ref) in references.items()}
    bg = load(background) if background is not None else None
    return analyse(load(sample), mass_g, refs, background=bg, lines=lines)


def calibrate_natural(spec: GammaSpectrum, guess: Calibration, lines: Sequence[float] = NATURAL_LINES,
                      min_significance: float = 8.0) -> Calibration:
    """Calibrate an environmental spectrum on its own natural lines (absorbs
    gain drift between measurements)."""
    return spec.calibrate(lines, guess=guess, min_significance=min_significance)


def water_content(wet_mass_g: ValueLike, dry_mass_g: ValueLike) -> Value:
    """Water content as mass of water / dry mass, the convention of
    :class:`~eprdating.Sediment`."""
    w, d = as_value(wet_mass_g), as_value(dry_mass_g)
    r = (w.value - d.value) / d.value
    s = np.hypot(w.sigma / d.value, w.value * d.sigma / d.value**2)
    return Value(r, float(s))


#: specific activities, Bq/kg per unit content (natural isotopic composition)
BQ_PER_KG = {"U": 12.35, "Th": 4.057, "K": 316.5}  # per µg/g U, µg/g Th, % K


__all__ = [
    "BQ_PER_KG",
    "GROUP_ELEMENT",
    "IAEA_RGK_1",
    "IAEA_RGTH_1",
    "IAEA_RGU_1",
    "LINES",
    "GammaResult",
    "GroupResult",
    "Line",
    "LineResult",
    "Reference",
    "analyse",
    "analyse_files",
    "calibrate_natural",
    "line_area",
    "net_rate",
    "water_content",
]
