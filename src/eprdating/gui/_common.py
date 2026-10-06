"""Pieces shared by the panels: files, value ± sigma fields, downloads, plots."""

from __future__ import annotations

import base64
import html
import tempfile
from pathlib import Path

import ipywidgets as w

STYLE = {"description_width": "initial"}
WIDE = w.Layout(width="auto")


def uploaded(upload: w.FileUpload) -> list[tuple[str, bytes]]:
    """``(name, content)`` of the files in a FileUpload, ipywidgets 7 or 8."""
    v = upload.value
    if isinstance(v, dict):  # ipywidgets 7: {name: {"metadata": ..., "content": bytes}}
        return [(name, bytes(d["content"])) for name, d in v.items()]
    return [(d["name"], bytes(d["content"])) for d in v]  # ipywidgets 8: tuple of dicts


def clear_upload(upload: w.FileUpload) -> None:
    try:
        upload.value = () if not isinstance(upload.value, dict) else {}
    except Exception:  # noqa: BLE001, S110 - some front ends keep the value read-only
        pass


class FilePool:
    """Files given by upload or from a folder, kept on disk so that readers
    that need a companion file (``.dat`` + ``.par``, ``.DSC`` + ``.DTA``)
    find it next to the main one."""

    def __init__(self) -> None:
        self.dir = Path(tempfile.mkdtemp(prefix="eprdating_"))
        self.paths: dict[str, Path] = {}

    def add_bytes(self, name: str, content: bytes) -> Path:
        p = self.dir / Path(name).name
        p.write_bytes(content)
        self.paths[p.name] = p
        return p

    def add_path(self, path: str | Path) -> Path:
        p = Path(path).expanduser()
        if not p.is_file():
            raise FileNotFoundError(p)
        self.paths[p.name] = p
        return p

    def add_folder(self, folder: str | Path, suffixes: tuple[str, ...]) -> list[Path]:
        d = Path(folder).expanduser()
        if not d.is_dir():
            raise NotADirectoryError(d)
        found = sorted(p for p in d.iterdir() if p.is_file() and p.suffix.lower() in suffixes)
        for p in found:
            self.paths[p.name] = p
        return found

    def names(self) -> list[str]:
        return sorted(self.paths)


def value_field(description: str, value: float = 0.0, sigma: float = 0.0, unit: str = "",
                width: str = "300px", label_width: str = "170px") -> w.HBox:
    """``value ± sigma unit`` input; read it with :func:`read_value`."""
    v = w.FloatText(value=value, description=description, style={"description_width": label_width},
                    layout=w.Layout(width=width))
    s = w.FloatText(value=sigma, description="±", style={"description_width": "14px"},
                    layout=w.Layout(width="110px"))
    box = w.HBox([v, s, w.Label(unit)])
    box.value_widget, box.sigma_widget = v, s
    return box


def read_value(box: w.HBox) -> tuple[float, float]:
    return float(box.value_widget.value), float(box.sigma_widget.value)


def _round(x: float, digits: int = 4) -> float:
    return float(f"{x:.{digits}g}") if x else 0.0


def set_value(box: w.HBox, value: float, sigma: float) -> None:
    """Fill a value ± sigma field, rounded to four significant digits."""
    box.value_widget.value = _round(float(value))
    box.sigma_widget.value = _round(float(sigma), 2)


def checkbox(description: str, value: bool = True) -> w.Checkbox:
    return w.Checkbox(value=value, description=description, indent=False, layout=w.Layout(width="auto"))


def button(description: str, icon: str = "", primary: bool = False, width: str = "190px") -> w.Button:
    return w.Button(description=description, icon=icon, button_style="primary" if primary else "",
                    layout=w.Layout(width=width))


def scale_of(values) -> tuple[float, str]:
    """Power of ten that brings the largest ``values`` to 1-1000, and its label."""
    import math

    m = max((abs(v) for v in values if math.isfinite(v)), default=0.0)
    if not m:
        return 1.0, ""
    k = 3 * math.floor(math.log10(m) / 3)
    return (10.0 ** -k, f" (×10^{k})") if k else (1.0, "")


def download_link(filename: str, content: str | bytes, label: str) -> str:
    """HTML link that saves ``content`` as ``filename`` (works in Jupyter and Colab)."""
    data = content.encode() if isinstance(content, str) else content
    b64 = base64.b64encode(data).decode()
    return (f'<a download="{html.escape(filename)}" href="data:application/octet-stream;base64,{b64}" '
            f'style="margin-right:1em">{html.escape(label)}</a>')


def message(text: str, kind: str = "info") -> str:
    color = {"info": "#555", "ok": "#1a7f37", "warn": "#9a6700", "error": "#cf222e"}[kind]
    return f'<div style="color:{color};margin:4px 0">{html.escape(text)}</div>'


def table_html(rows: list[dict], columns: list[tuple[str, str]]) -> str:
    """Small HTML table; ``columns`` = ``(key, header)``."""
    head = "".join(f"<th style='text-align:right;padding:2px 8px'>{html.escape(h)}</th>" for _, h in columns)
    body = []
    for r in rows:
        cells = "".join(f"<td style='text-align:right;padding:2px 8px'>{html.escape(str(r.get(k, '')))}</td>"
                        for k, _ in columns)
        body.append(f"<tr>{cells}</tr>")
    return f"<table style='border-collapse:collapse;font-size:13px'><tr>{head}</tr>{''.join(body)}</table>"


def show_figure(out: w.Output, *figs) -> None:
    """Draw matplotlib figures in an Output widget, replacing what it held."""
    import matplotlib.pyplot as plt
    from IPython.display import display

    out.clear_output(wait=True)
    with out:
        for fig in figs:
            display(fig)
    for fig in figs:
        plt.close(fig)


def colab_hint() -> str:
    """A line for Google Colab users, empty elsewhere."""
    import sys

    if "google.colab" not in sys.modules:
        return ""
    return ("<small>In Colab you can also drop the files in the Files panel (left) or mount Google Drive, "
            "and load the folder (e.g. <code>/content</code>).</small>")
