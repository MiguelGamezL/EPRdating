"""Write the API reference (docs/api/*.md) from the docstrings.

    python tools/gen_api_docs.py

One page per module, listing what the module exports in ``__all__`` with
its signature and docstring. The pages are plain Markdown, readable on GitHub
and by MkDocs. Re-run after changing a public function.
"""

from __future__ import annotations

import dataclasses
import importlib
import inspect
import re
import textwrap
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "docs" / "api"

PAGES = [
    ("eprdating", "Top level", "Everything most scripts need, importable as `from eprdating import ...`."),
    ("eprdating.doseresponse", "Dose-response and De", None),
    ("eprdating.age", "Ages", None),
    ("eprdating.usesr", "US-ESR", None),
    ("eprdating.dose_rate", "Dose rates", None),
    ("eprdating.history", "Environmental histories", None),
    ("eprdating.uptake", "Uranium uptake", None),
    ("eprdating.series", "U-series ingrowth", None),
    ("eprdating.onegroup", "One-group beta attenuation", None),
    ("eprdating.beta", "Fixed beta factors", None),
    ("eprdating.alpha", "Alpha efficiency", None),
    ("eprdating.spectra.io", "EPR spectra: reading", None),
    ("eprdating.spectra.bruker", "EPR spectra: Bruker files", None),
    ("eprdating.spectra.preprocess", "EPR spectra: pre-processing", None),
    ("eprdating.spectra.intensity", "EPR spectra: scalar intensities", None),
    ("eprdating.spectra.deconvolution", "EPR spectra: template and component fits", None),
    ("eprdating.spectra.epraya_backend", "EPR spectra: EPRAYA simulations", None),
    ("eprdating.gamma.readers", "Gamma: reading spectra", None),
    ("eprdating.gamma.spectrum", "Gamma: calibration", None),
    ("eprdating.gamma.comparative", "Gamma: comparative method", None),
    ("eprdating.plot", "Plotting", None),
]

_ROLE = re.compile(r":(?:func|class|mod|meth|attr|data):`~?([^`]+)`")


def _clean(doc: str) -> str:
    doc = inspect.cleandoc(doc or "")
    doc = _ROLE.sub(lambda m: f"`{m.group(1).split('.')[-1]}`", doc)
    doc = re.sub(r"``([^`]+)``", r"`\1`", doc)
    # reST section underlines -> bold headings
    doc = re.sub(r"^(\S[^\n]*)\n[-=~]{3,}\s*$", r"**\1**\n", doc, flags=re.MULTILINE)
    # literal blocks introduced by "::" and reST simple tables -> fenced code
    out, lines, i = [], doc.splitlines(), 0
    rule = re.compile(r"^\s*=+(\s+=+)+\s*$")
    while i < len(lines):
        ln = lines[i]
        if rule.match(ln):
            block = [ln]
            i += 1
            while i < len(lines):
                block.append(lines[i])
                i += 1
                if rule.match(block[-1]) and (i >= len(lines) or not lines[i].strip()):
                    break
            body = "\n".join(b for b in block if not rule.match(b)).rstrip().replace("`", "")
            out += ["", "```text", body, "```", ""]
            continue
        if ln.rstrip().endswith("::"):
            out.append(ln.rstrip()[:-1] if not ln.rstrip().endswith(" ::") else ln.rstrip()[:-3])
            i += 1
            block = []
            while i < len(lines) and (not lines[i].strip() or lines[i].startswith((" ", "\t"))):
                block.append(lines[i])
                i += 1
            while block and not block[-1].strip():
                block.pop()
            out += ["", "```", textwrap.dedent("\n".join(block)).strip("\n"), "```", ""]
            continue
        out.append(ln)
        i += 1
    return _params_to_list("\n".join(out))


def _ann(a) -> str:
    if a is inspect.Parameter.empty:
        return ""
    return a if isinstance(a, str) else inspect.formatannotation(a)


def _sig(obj) -> str:
    try:
        sig = inspect.signature(obj)
    except (TypeError, ValueError):
        return ""
    parts, star = [], False
    for p in sig.parameters.values():
        if p.name.startswith("_"):
            continue
        if p.kind is p.KEYWORD_ONLY and not star:
            parts.append("*")
            star = True
        txt = ("*" if p.kind is p.VAR_POSITIONAL else "**" if p.kind is p.VAR_KEYWORD else "") + p.name
        if p.annotation is not p.empty:
            txt += f": {_ann(p.annotation)}"
        if p.default is not p.empty:
            d = re.sub(r"<function (\w+)[^>]*>", r"\1", repr(p.default))
            txt += f" = {d}"
        parts.append(txt)
    ret = f" -> {_ann(sig.return_annotation)}" if sig.return_annotation not in (sig.empty, "None", None) else ""
    return f"({', '.join(parts)}){ret}"


_PARAM = re.compile(r"^([A-Za-z_][\w]*(?:\s*,\s*[A-Za-z_][\w]*)*)\s*:\s+(\S.*)$")
_NAME_ONLY = re.compile(r"^([A-Za-z_][\w]*(?:\s*,\s*[A-Za-z_][\w]*)*)\s*$")


def _params_to_list(doc: str) -> str:
    """'name : text' entries and numpy-style 'name' + indented text -> bullets."""
    lines, out, i, fence = doc.splitlines(), [], 0, False
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("```"):
            fence = not fence
        if fence:
            out.append(ln)
            i += 1
            continue
        m = _PARAM.match(ln)
        nm = _NAME_ONLY.match(ln)
        nxt = lines[i + 1] if i + 1 < len(lines) else ""
        if m or (nm and nxt.startswith((" ", "\t")) and nxt.strip()):
            names, text = (m.group(1), m.group(2)) if m else (nm.group(1), "")
            i += 1
            while i < len(lines) and lines[i].startswith((" ", "\t")) and lines[i].strip():
                text += " " + lines[i].strip()
                i += 1
            label = ", ".join(f"`{n.strip()}`" for n in names.split(","))
            out.append(f"- {label}: {text.strip()}")
            continue
        out.append(ln)
        i += 1
    return "\n".join(out)


def _entry(name: str, obj) -> str:
    parts = []
    if inspect.isclass(obj):
        kind = "dataclass" if dataclasses.is_dataclass(obj) else "class"
        parts.append(f"### `{name}`\n\n*{kind}* `{name}{_sig(obj)}`\n")
        parts.append(_clean(obj.__doc__))
        methods = []
        for m, f in inspect.getmembers(obj):
            if m.startswith("_") or m in getattr(obj, "__dataclass_fields__", {}):
                continue
            if inspect.isfunction(f) and f.__qualname__.startswith(obj.__name__):
                doc = _clean(f.__doc__).split("\n\n")[0] if f.__doc__ else ""
                methods.append(f"- `{m}{_sig(f)}` — {doc}" if doc else f"- `{m}{_sig(f)}`")
            elif isinstance(inspect.getattr_static(obj, m), property):
                doc = _clean(f.__doc__).split("\n\n")[0] if getattr(f, "__doc__", None) else ""
                methods.append(f"- `{m}` *(property)*" + (f" — {doc}" if doc else ""))
        if methods:
            parts.append("\n**Members**\n\n" + "\n".join(methods))
    elif callable(obj):
        parts.append(f"### `{name}`\n\n`{name}{_sig(obj)}`\n")
        parts.append(_clean(obj.__doc__))
    else:
        r = repr(obj)
        r = r if len(r) < 300 else r[:300] + " …"
        parts.append(f"### `{name}`\n\n```python\n{name} = {r}\n```")
    return "\n".join(p for p in parts if p) + "\n"


def page(modname: str, title: str, intro: str | None) -> str:
    mod = importlib.import_module(modname)
    names = list(getattr(mod, "__all__", []))
    if not names:  # public objects defined in the module itself
        names = sorted(n for n, o in vars(mod).items()
                       if not n.startswith("_") and getattr(o, "__module__", None) == modname
                       and (inspect.isclass(o) or inspect.isfunction(o)))
    out = [f"# {title}\n", f"`{modname}`\n"]
    if intro:
        out.append(intro + "\n")
    if mod.__doc__ and modname != "eprdating":
        out.append(_clean(mod.__doc__) + "\n")
    out.append("**Contents:** " + ", ".join(f"[`{n}`](#{n.lower()})" for n in names) + "\n")
    for n in names:
        out.append(_entry(n, getattr(mod, n)))
    return "\n".join(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    index = ["# API reference\n", "Generated from the docstrings by `tools/gen_api_docs.py`.\n"]
    for modname, title, intro in PAGES:
        fname = modname.replace("eprdating.", "").replace("eprdating", "top-level") + ".md"
        (OUT / fname).write_text(page(modname, title, intro), encoding="utf-8")
        index.append(f"- [{title}]({fname}) — `{modname}`")
    (OUT / "index.md").write_text("\n".join(index) + "\n", encoding="utf-8")
    print(f"wrote {len(PAGES) + 1} pages to {OUT}")


if __name__ == "__main__":
    main()
