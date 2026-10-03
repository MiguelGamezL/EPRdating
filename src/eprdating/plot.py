"""Figures for every step: spectra, dose-response, dose rates, ages, gamma.

Needs matplotlib (``pip install "eprdating[plot]"``). Every function draws
on the axes passed as ``ax`` (or a new figure) and returns the axes, so the
figures can be combined and restyled with ordinary matplotlib calls.

==============================  ==========================================
:func:`plot_spectra`            a series of cw-EPR spectra, stacked
:func:`plot_spectrum`           one spectrum with its individual scans
:func:`plot_dose_response`      intensity vs added dose, fit and De
:func:`plot_dose_rate`          contribution of each dose-rate component
:func:`plot_age_distribution`   Monte Carlo age distribution
:func:`plot_gamma_spectrum`     HPGe spectrum with the analysed lines
:func:`plot_gamma_lines`        content from each gamma line, per group
==============================  ==========================================

Colours: series follow a fixed categorical order, ordered series (doses) a
single-hue ramp, text stays in neutral inks, and grids are hairlines.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

# palette (validated for colour-vision deficiencies; see docs/plotting.md)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#8a8984"
GRID = "#e4e3df"
CATEGORICAL = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
BLUE_RAMP = ("#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95",
             "#104281", "#0d366b")  # ordinal steps 250 -> 700
RAW = "#c9c8c3"


def _plt():
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:  # pragma: no cover
        raise ImportError('plotting needs matplotlib: pip install "eprdating[plot]"') from exc
    return plt


def _axes(ax, figsize=(6.4, 4.2)):
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=figsize, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_3)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK_2, labelcolor=INK_2, width=0.8)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.title.set_color(INK)
    return ax


def _grid(ax, axis="y"):
    ax.grid(True, axis=axis, color=GRID, linewidth=0.6, linestyle="-")
    ax.set_axisbelow(True)


def ramp(n: int) -> list[str]:
    """``n`` colours from the single-hue ramp, light to dark (for doses)."""
    if n <= 1:
        return [BLUE_RAMP[5]]
    idx = np.linspace(0, len(BLUE_RAMP) - 1, n).round().astype(int)
    return [BLUE_RAMP[i] for i in idx]


# --------------------------------------------------------------------- EPR


def _xy(s):
    if isinstance(s, tuple):
        return np.asarray(s[0], float), np.asarray(s[1], float), None
    return np.asarray(s.B, float), np.asarray(s.y, float), np.asarray(s.scans, float)


def plot_spectra(spectra: Sequence, labels: Sequence[str] | None = None, *, window: tuple[float, float] | None = None,
                 offset: float | None = None, fits: Sequence | None = None, baseline: bool = True, ax=None,
                 title: str | None = None):
    """Stacked spectra of a dose series, light to dark with the dose.

    spectra : :class:`~eprdating.spectra.Spectrum` objects or ``(B, y)`` pairs.
    labels  : text at the right of each trace (e.g. ``"40 Gy"``).
    window  : field range (mT) to show.
    offset  : vertical spacing; default 1.2 × the largest peak-to-peak.
    fits    : optional fitted curves, ``(B, y)`` pairs or fit results with a
              ``fitted`` array on the same window (drawn over the data).
    baseline: subtract each trace's mean inside the window.
    """
    ax = _axes(ax, (5.6, 1.2 + 0.55 * len(spectra)))
    data = []
    for s in spectra:
        B, y, _ = _xy(s)
        m = np.ones_like(B, bool) if window is None else (B >= window[0]) & (B <= window[1])
        B, y = B[m], y[m]
        data.append((B, y, y.mean() if baseline else 0.0))
    if offset is None:
        offset = 1.2 * max(np.ptp(y) for _, y, _ in data)
    colors = ramp(len(data))
    for i, ((B, y, mean), c) in enumerate(zip(data, colors)):
        base = i * offset - mean
        f = None if fits is None else fits[i]
        if f is not None:
            fb, fy = (np.asarray(a, float) for a in f) if isinstance(f, tuple) else (B, np.asarray(f.fitted, float))
            ax.plot(B, y + base, color=RAW, lw=0.9)
            ax.plot(fb, fy + base, color=c, lw=1.6, solid_capstyle="round")
        else:
            ax.plot(B, y + base, color=c, lw=1.4, solid_capstyle="round")
        base = i * offset
        if labels is not None:
            ax.annotate(labels[i], (B[-1], base), xytext=(6, 0), textcoords="offset points", va="center",
                        fontsize=9, color=INK_2, annotation_clip=False)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Magnetic field (mT)")
    ax.set_xlim(data[0][0][0], data[0][0][-1])
    if title:
        ax.set_title(title, loc="left", fontsize=11)
    return ax


def plot_spectrum(spectrum, *, window: tuple[float, float] | None = None, show_scans: bool = True,
                  fit=None, ax=None, title: str | None = None):
    """One spectrum: individual scans in grey, their average, and an
    optional fitted curve (``(B, y)`` or a fit result on the same window)."""
    ax = _axes(ax)
    B, y, scans = _xy(spectrum)
    m = np.ones_like(B, bool) if window is None else (B >= window[0]) & (B <= window[1])
    if show_scans and scans is not None and scans.shape[0] > 1:
        for k, sc in enumerate(scans):
            ax.plot(B[m], sc[m], color=RAW, lw=0.7, label="scans" if k == 0 else None)
    ax.plot(B[m], y[m], color=CATEGORICAL[0], lw=1.5, label="average" if show_scans else None)
    if fit is not None:
        fb, fy = (np.asarray(a, float) for a in fit) if isinstance(fit, tuple) else (B[m], np.asarray(fit.fitted))
        ax.plot(fb, fy, color=CATEGORICAL[1], lw=1.6, label="fit")
    _grid(ax, "both")
    ax.set_xlabel("Magnetic field (mT)")
    ax.set_ylabel("dχ″/dB (a.u.)")
    handles, _ = ax.get_legend_handles_labels()
    if len(handles) > 1:
        ax.legend(frameon=False, fontsize=9, labelcolor=INK_2)
    name = title if title is not None else getattr(spectrum, "name", "")
    if name:
        ax.set_title(name, loc="left", fontsize=11)
    return ax


def plot_dose_response(drc, *, excluded: Sequence[tuple] | None = None, ax=None, dose_unit: str = "Gy",
                       intensity_label: str = "EPR intensity (a.u.)", show_De: bool = True,
                       title: str | None = None):
    """Additive-dose points, fitted curve extrapolated to zero and De.

    drc      : :class:`~eprdating.DoseResponseResult`.
    excluded : points left out of the fit, as ``(dose, intensity[, sigma])``;
               drawn hollow.
    """
    ax = _axes(ax)
    D, I = drc.dose, drc.intensity
    s = drc.sigma
    lo = min(-drc.De - (drc.De_sigma if np.isfinite(drc.De_sigma) else 0.0), 0.0)
    span = max(D.max() - lo, 1.0)
    x = np.linspace(lo - 0.08 * span, D.max() + 0.05 * span, 300)
    ax.axhline(0, color=INK_3, lw=0.8)
    ax.axvline(0, color=GRID, lw=0.8)
    ax.plot(x, drc.predict(x), color=CATEGORICAL[0], lw=1.6, label=f"{drc.model} fit")
    ax.errorbar(D, I, yerr=s, fmt="o", ms=5.5, color=INK, mfc=INK, mec=SURFACE, mew=1.2, ecolor=INK_2,
                elinewidth=0.9, capsize=0, label="fitted points", zorder=3)
    if excluded:
        ex = np.array([tuple(e) + (np.nan,) * (3 - len(e)) for e in excluded], float)
        ax.errorbar(ex[:, 0], ex[:, 1], yerr=None if np.all(np.isnan(ex[:, 2])) else ex[:, 2], fmt="o", ms=5.5,
                    mfc=SURFACE, mec=INK_2, mew=1.1, ecolor=INK_3, elinewidth=0.9, capsize=0,
                    label="excluded", zorder=3)
    if show_De:
        ax.errorbar([-drc.De], [0.0], xerr=[[drc.De_sigma], [drc.De_sigma]], fmt="D", ms=6,
                    color=CATEGORICAL[1], mec=SURFACE, mew=1.2, ecolor=CATEGORICAL[1], elinewidth=1.4, capsize=0,
                    label="−De", zorder=4)
        ax.annotate(f"De = {drc.De:.3g} ± {drc.De_sigma:.2g} {dose_unit}", (-drc.De, 0.0), xytext=(-4, 8),
                    textcoords="offset points", ha="right", va="bottom", fontsize=9, color=INK,
                    bbox={"boxstyle": "square,pad=0.15", "fc": SURFACE, "ec": "none"})
    _grid(ax, "y")
    ax.set_xlabel(f"Added dose ({dose_unit})")
    ax.set_ylabel(intensity_label)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK_2, loc="upper left")
    if title:
        ax.set_title(title, loc="left", fontsize=11)
    return ax


def plot_dose_rate(result, *, ax=None, unit: str = "Gy/ka", title: str | None = None):
    """Time-averaged dose rate of each component at the age found.

    result : :class:`~eprdating.age.AgeResult` (``sample.age()``), or a
             ``{name: rate}`` mapping.
    """
    ax = _axes(ax, (6.0, 3.2))
    if isinstance(result, dict):
        items = {k: v for k, v in result.items() if v}
    else:
        items = {k: result.accumulated[k] / result.age for k in result.components if result.accumulated[k]}
    names = list(items)[::-1]
    vals = [items[k] for k in names]
    y = np.arange(len(names))
    ax.barh(y, vals, height=0.55, color=CATEGORICAL[0])
    total = sum(vals)
    for yi, v in zip(y, vals):
        ax.annotate(f"{v:.3g}  ({100 * v / total:.0f} %)", (v, yi), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=9, color=INK_2)
    ax.set_yticks(y, names)
    ax.tick_params(axis="y", length=0)
    _grid(ax, "x")
    ax.set_xlabel(f"Dose rate ({unit})")
    ax.set_xlim(0, max(vals) * 1.35)
    ax.set_title(title if title is not None else f"Total {total:.3g} {unit}", loc="left", fontsize=11)
    return ax


def plot_age_distribution(mc, *, ax=None, bins: int = 50, unit: str = "ka", reference: tuple | None = None,
                          title: str | None = None):
    """Histogram of Monte Carlo ages with the 68 % interval.

    mc        : :class:`~eprdating.AgeMC` or a US-ESR Monte Carlo result.
    reference : optional ``(low, high)`` range to compare with (e.g. the
                archaeological context), drawn as a band.
    """
    ax = _axes(ax)
    ages = np.asarray(getattr(mc, "samples", getattr(mc, "ages", mc)), float)
    ages = ages[np.isfinite(ages)]
    lo, hi = np.quantile(ages, [0.16, 0.84])
    if reference is not None:
        ax.axvspan(*reference, color=CATEGORICAL[2], alpha=0.12, lw=0, label="reference range")
    ax.hist(ages, bins=bins, color=CATEGORICAL[0], alpha=0.85, rwidth=0.9, label="Monte Carlo")
    ax.axvspan(lo, hi, color=CATEGORICAL[0], alpha=0.10, lw=0, label="68 % interval")
    med = float(np.median(ages))
    ax.axvline(med, color=INK, lw=1.0)
    ax.text(0.02, 0.97, f"median {med:.3g} (+{hi - med:.2g} / −{med - lo:.2g}) {unit}", transform=ax.transAxes,
            ha="left", va="top", fontsize=9, color=INK)
    _grid(ax, "y")
    ax.set_xlabel(f"Age ({unit})")
    ax.set_ylabel("Draws")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK_2, loc="upper right")
    if title:
        ax.set_title(title, loc="left", fontsize=11)
    return ax


# -------------------------------------------------------------------- gamma

_GROUP_COLOR = {"K": CATEGORICAL[0], "Ra226": CATEGORICAL[1], "U238": CATEGORICAL[3], "Th232": CATEGORICAL[2]}
_GROUP_NAME = {"K": "⁴⁰K", "Ra226": "²²⁶Ra (U)", "U238": "²³⁸U", "Th232": "²³²Th"}


def plot_gamma_spectrum(spec, *, ax=None, lines=None, energy_range: tuple[float, float] | None = (30, 2700),
                        log: bool = True, label_lines: bool = True, title: str | None = None):
    """Calibrated HPGe spectrum with the analysed lines marked by group."""
    from .gamma.comparative import LINES

    ax = _axes(ax, (7.5, 3.6))
    E = spec.energy
    m = np.ones_like(E, bool) if energy_range is None else (E >= energy_range[0]) & (E <= energy_range[1])
    ax.plot(E[m], np.maximum(spec.counts[m], 0.5), color=INK_2, lw=0.6, drawstyle="steps-mid")
    if log:
        ax.set_yscale("log")
    seen = set()
    shown = sorted((ln for ln in (LINES if lines is None else lines)
                    if not energy_range or energy_range[0] <= ln.energy <= energy_range[1]), key=lambda q: q.energy)
    span = (energy_range[1] - energy_range[0]) if energy_range else float(np.ptp(E))
    last = {}  # label row -> energy of its last label
    for ln in shown:
        c = _GROUP_COLOR.get(ln.group, CATEGORICAL[4])
        ax.axvline(ln.energy, color=c, lw=0.8, alpha=0.7,
                   label=_GROUP_NAME.get(ln.group, ln.group) if ln.group not in seen else None)
        seen.add(ln.group)
        if label_lines:  # stagger labels of close lines on up to three rows
            row = next((r for r in range(3) if ln.energy - last.get(r, -np.inf) > 0.025 * span), None)
            if row is None:
                continue
            last[row] = ln.energy
            ax.annotate(f"{ln.energy:.0f}", (ln.energy, 1.0), xycoords=("data", "axes fraction"),
                        xytext=(2, -2 - 22 * row), textcoords="offset points", rotation=90, va="top",
                        fontsize=7, color=INK_3)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts per channel")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK_2, ncol=4, loc="lower left")
    ax.set_title(title if title is not None else spec.name, loc="left", fontsize=11)
    return ax


def plot_gamma_lines(result, *, ax=None, title: str | None = None):
    """Content obtained from each line, grouped, with the group mean ± 1σ.

    Lines that disagree with their group point to interferences or
    self-absorption; the 238U group against the 226Ra group shows
    disequilibrium.
    """
    plt = _plt()
    groups = [g for g in ("K", "Ra226", "U238", "Th232") if g in result.groups]
    if ax is None:
        _fig, axs = plt.subplots(1, len(groups), figsize=(2.6 * len(groups) + 0.6, 3.4), facecolor=SURFACE,
                                gridspec_kw={"width_ratios": [max(len(result.groups[g].lines), 1.5) for g in groups]})
        axs = np.atleast_1d(axs)
    else:
        axs = np.atleast_1d(ax)
    unit = {"K": "%", "U": "µg/g", "Th": "µg/g"}
    for a, g in zip(axs, groups):
        a = _axes(a)
        gr = result.groups[g]
        x = np.arange(len(gr.lines))
        c = _GROUP_COLOR[g]
        a.axhspan(gr.content.value - gr.content.sigma, gr.content.value + gr.content.sigma, color=c, alpha=0.12,
                  lw=0)
        a.axhline(gr.content.value, color=c, lw=1.4)
        a.errorbar(x, [lr.content for lr in gr.lines], yerr=[lr.content_sigma for lr in gr.lines], fmt="o", ms=5,
                   color=INK, mec=SURFACE, mew=1.0, ecolor=INK_2, elinewidth=0.9, capsize=0, zorder=3)
        a.set_xticks(x, [f"{lr.line.energy:.0f}" for lr in gr.lines], rotation=0, fontsize=8)
        a.set_xlim(-0.6, len(x) - 0.4)
        a.set_xlabel("Line (keV)", fontsize=9)
        a.set_ylabel(f"{gr.element} ({unit[gr.element]})", fontsize=9)
        a.set_title(f"{_GROUP_NAME[g]}: {gr.content.value:.3g} ± {gr.content.sigma:.2g}", loc="left", fontsize=10)
        _grid(a, "y")
    if title:
        axs[0].figure.suptitle(title, x=0.02, ha="left", fontsize=11, color=INK)
    axs[0].figure.tight_layout()
    return axs


__all__ = [
    "plot_age_distribution",
    "plot_dose_rate",
    "plot_dose_response",
    "plot_gamma_lines",
    "plot_gamma_spectrum",
    "plot_spectra",
    "plot_spectrum",
    "ramp",
]
