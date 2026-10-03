"""Batch-run DATA.EXE (Grün 2009) as a black box through DOSBox-X.

Each case is written as an .EPR file, loaded with <F5>, calculated with <F1>
and printed with <F10> (the printer port is captured to a file). A fresh DOS
session is started for every case.

    DATA_DIR=/path/with/DATA.EXE DISPLAY=:99 python run_data.py cases.json results.json
"""
import json
import pathlib
import re
import sys
import zlib

from drv import D, key, printed, shot, start, stop, typ

# the 19 fields of an .EPR file, in order; all are value,error pairs except
# the "beta only" switch (Y: external gamma given directly)
ORDER = [("De", 2), ("U_en", 2), ("r48", 2), ("k", 2), ("U_d1", 2), ("U_d2", 2), ("W_d", 2), ("Rn", 2),
         ("thick", 2), ("rem1", 2), ("rem2", 2), ("dens", 2), ("U_s", 2), ("Th_s", 2), ("K_s", 2), ("W_s", 2),
         ("beta_only", 1), ("ext_gamma", 2), ("depth", 2)]


def write_epr(name, c):
    lines = ['"SITE"', '"E",19']
    for k, n in ORDER:
        v = c[k]
        if n == 1:
            lines.append(f'"{v}"')
        else:
            v = v if isinstance(v, (list, tuple)) else (v, 0.0)
            lines.append(f"{float(v[0])!r},{float(v[1])!r}")
    (D / "esr" / "epr" / f"{name}.EPR").write_bytes(("\r\n".join(lines) + "\r\n").encode())


NUM = r"(-?\d*\.?\d+)±(\d*\.?\d+)"
ROWS = (("τ-SED", "gamma_sed"), ("int.", "internal"), ("ß-DE1", "beta_de1"), ("ß-DE2", "beta_de2"),
        ("ß-SED", "beta_sed"), ("TOTAL", "total"), ("AGE", "age"))


def parse(txt):
    """Dose rates (µGy/a) and ages (ka) for EU and LU, and the beta factors."""
    out = {}
    for ln in txt.splitlines():
        for row, key_ in ROWS:
            m = re.search(r"║\s*" + re.escape(row) + r"\s*│\s*" + NUM + r"\s+" + NUM, ln)
            if m:
                out[key_] = {"EU": [float(m.group(1)), float(m.group(2))],
                             "LU": [float(m.group(3)), float(m.group(4))]}
    lines = txt.splitlines()
    head = next((ln for ln in lines if "CORRECTION" in ln and "SED" in ln), None)
    if head:
        cols = {c: head.index(c) for c in ("SED", "DE-1", "DE-2")}
        for ln in lines:
            for lab, key_ in (("URANIUM", "U"), ("THORIUM", "Th"), ("K-40", "K")):
                if lab in ln:
                    s0 = ln.index("│", ln.index(lab)) + 1
                    seg = ln[s0: ln.index("│", s0)]
                    vals = {}
                    for m in re.finditer(r"\.\d+|\d+\.\d*", seg):
                        pos = s0 + m.start()
                        col = min(cols, key=lambda c, pos=pos: abs(cols[c] - pos))
                        vals[col] = float(m.group())
                    out[f"beta_corr_{key_}"] = vals
    return out


def run_case(name, c):
    short = "C" + format(zlib.crc32(name.encode()) % 10**6, "06d")  # DOS 8.3 name, no dots
    write_epr(short, c)
    key("F5", pause=1.2)
    typ("C:\\esr\\epr\\", 0.3)
    key("Return", pause=1.2)
    typ(short, 0.3)
    key("Return", pause=1.2)
    key("Right", pause=0.6)
    key("Return", pause=1.5)
    key("Return", pause=0.8)  # leave the De edit field unchanged
    key("F1", pause=2.5)
    txt = printed()
    res = parse(txt)
    res["raw"] = txt
    return res


if __name__ == "__main__":
    cases = json.loads(pathlib.Path(sys.argv[1]).read_text())
    outp = pathlib.Path(sys.argv[2])
    done = json.loads(outp.read_text()) if outp.exists() else {}
    for name, c in cases.items():
        if name in done and "age" in done[name]:
            continue
        for _ in range(2):
            start()
            try:
                r = run_case(name, c)
            except Exception as e:  # noqa: BLE001
                r = {"error": repr(e)}
            if "age" in r:
                break
            shot(f"fail_{name}")
        done[name] = r
        outp.write_text(json.dumps(done, indent=1, ensure_ascii=False))
        print(name, {k: v for k, v in r.items() if k != "raw"}, flush=True)
    stop()
