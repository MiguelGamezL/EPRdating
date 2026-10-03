"""Batch-run the U-series/ESR part of DATA (Grün 2009) through DOSBox-X.

<F6> "READ ENAMEL INTO U-SERIES" loads an .EPR file into a screen that adds
234U/238U and 230Th/234U for the enamel and each dentine side. These cannot
be stored in the .EPR file, so they are typed in: a field is edited with
<Enter> value <Enter> <Enter> error <Enter>, after which the cursor moves to
the next row. <F1> calculates; <F10> prints EU, LU, US-ESR (with the p
values) and CS-US results. The printout echoes every input, which is checked
against the case so that a mistyped field is caught.

    DATA_DIR=/path/with/DATA.EXE DISPLAY=:99 python run_useries.py cases.json results.json

A case is a run_data case plus "useries": {"enamel"|"de1"|"de2": [[234U/238U,
err], [230Th/234U, err]]}.
"""
import hashlib
import json
import pathlib
import re
import sys
import time
import zlib

from drv import ENV, key, printed, shot, start, stop, typ
from run_data import NUM, ROWS, write_epr

# screen rows: 0 De, 1 U, 2 234/238 (enamel), 3 230/234, 4 k, 5 U DE1, 6 234/238,
# 7 230/234, 8 U DE2, 9 234/238, 10 230/234, ...
TISSUE_ROWS = {"enamel": 2, "de1": 6, "de2": 9}


def edit(value, error, pause=0.5):
    key("Return", pause=pause)
    typ(f"{value:g}", 0.2)
    key("Return", pause=pause)
    key("Return", pause=pause)
    typ(f"{error:g}", 0.2)
    key("Return", pause=pause)


def type_useries(us):
    row = 0
    for tissue in ("enamel", "de1", "de2"):
        target = TISSUE_ROWS[tissue]
        for _ in range(target - row):
            key("Down", pause=0.3)
        (r48, r48e), (th, the) = us[tissue]
        edit(r48, r48e)
        edit(th, the)
        row = target + 2


def parse_useries(txt):
    """EU/LU block (as run_data.parse), then the US-ESR and CS-US block."""
    head, _, tail = txt.partition("U-SERIES/ESR")
    out = {}
    for ln in head.splitlines():
        for row, key_ in ROWS:
            m = re.search(r"║\s*" + re.escape(row) + r"\s*│\s*" + NUM + r"\s+" + NUM, ln)
            if m:
                out[key_] = {"EU": [float(m.group(1)), float(m.group(2))],
                             "LU": [float(m.group(3)), float(m.group(4))]}
    us = {}
    for ln in tail.splitlines():
        for row, key_ in (("int.", "internal"), ("ß-DE1", "beta_de1"), ("ß-DE2", "beta_de2"), ("TOTAL", "total")):
            m = re.search(r"║\s*" + re.escape(row) + r"\s*│\s*" + NUM + r"(?:\s+" + NUM + r")?", ln)
            if m:
                us[key_] = [float(m.group(1)), float(m.group(2))]
                if m.group(3):
                    us[f"p_{key_}"] = [float(m.group(3)), float(m.group(4))]
        m = re.search(r"║\s*AGE\s*│\s*(\d*\.?\d+)\s*\+(\d*\.?\d+)\s*-(\d*\.?\d+)", ln)
        if m:
            us["age"] = [float(m.group(1)), float(m.group(2)), float(m.group(3))]
        m = re.search(r"║\s*CS-US\s*│\s*" + NUM + r"\s+" + NUM, ln)
        if m:
            out["csus"] = {"age": [float(m.group(1)), float(m.group(2))],
                           "col2": [float(m.group(3)), float(m.group(4))]}
        if "NO RESULT" in ln:
            us["no_result"] = ln.split("│")[1].strip(" ║")
    out["usesr"] = us
    # echoed inputs: U-234/U-238 and Th-230/U-234 lines, in screen order
    out["echo"] = [[float(a), float(b)] for a, b in
                   re.findall(r"(?:U-234/U-238|Th-230/U-234)\s+(-?\d*\.?\d+)±\s*(\d*\.?\d+)", head)]
    return out


def screen_hash():
    import subprocess

    png = subprocess.run(["import", "-window", "root", "png:-"], env=ENV, capture_output=True, check=False).stdout
    return hashlib.md5(png).hexdigest()


def wait_until_still(timeout=300.0, still=6.0, every=1.5):
    """US-ESR is iterative (a counter runs on screen): wait until the screen
    has not changed for ``still`` seconds."""
    t0, last, since = time.time(), None, time.time()
    while time.time() - t0 < timeout:
        h = screen_hash()
        if h != last:
            last, since = h, time.time()
        elif time.time() - since >= still:
            return True
        time.sleep(every)
    return False


def expected_echo(us):
    return [v for t in ("enamel", "de1", "de2") for v in us[t]]


def run_case(name, c):
    short = "U" + format(zlib.crc32(name.encode()) % 10**6, "06d")
    write_epr(short, c)
    key("F6", pause=1.2)
    typ("C:\\esr\\epr\\", 0.3)
    key("Return", pause=1.2)
    typ(short, 0.3)
    key("Return", pause=1.2)
    key("Right", pause=0.6)
    key("Return", pause=1.5)
    key("Return", pause=2.0)
    type_useries(c["useries"])
    key("F1", pause=2.0)
    wait_until_still()
    txt = printed()
    res = parse_useries(txt)
    res["raw"] = txt
    want = expected_echo(c["useries"])
    ok = len(res["echo"]) == 6 and all(
        abs(a - x) < 5e-5 and abs(b - y) < 5e-5 for (a, b), (x, y) in zip(res["echo"], want, strict=True))
    res["inputs_ok"] = ok
    return res


if __name__ == "__main__":
    cases = json.loads(pathlib.Path(sys.argv[1]).read_text())
    outp = pathlib.Path(sys.argv[2])
    done = json.loads(outp.read_text()) if outp.exists() else {}
    for name, c in cases.items():
        if name in done and done[name].get("inputs_ok") and "age" in done[name]:
            continue
        for _ in range(3):
            start()
            try:
                r = run_case(name, c)
            except Exception as e:  # noqa: BLE001
                r = {"error": repr(e)}
            if r.get("inputs_ok") and "age" in r:
                break
            shot(f"fail_{name}")
        done[name] = r
        outp.write_text(json.dumps(done, indent=1, ensure_ascii=False))
        print(name, r.get("inputs_ok"), r.get("usesr"), r.get("csus"), flush=True)
    stop()
