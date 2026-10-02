"""Drive the ROSY 2.0 GUI under Wine + Xvfb to compute reference ages in batch.

ROSY is NOT distributed with EPRdating. Point ``--prefix`` to a 32-bit Wine
prefix whose ``drive_c/ROSY`` folder holds your own copy of ROSY
(``rosy.bat`` and friends). Results go to ``drive_c/out``.

Usage::

    Xvfb :99 -screen 0 1280x900x24 &
    python run_rosy.py cases.json results.json --prefix ~/.wine-rosy --restart

``cases.json`` maps a case id to a dict of inputs (see ``rosy_io.DEFAULTS``).
The click coordinates assume ROSY's default window position on a 1280x900
screen. Only ROSY's inputs and outputs are used (black-box comparison).
"""

import argparse
import json
import os
import subprocess
import time

import rosy_io

HERE = os.path.dirname(os.path.abspath(__file__))
ENV = dict(os.environ, DISPLAY=os.environ.get("DISPLAY", ":99"))


def x(*args):
    subprocess.run(["xdotool", *map(str, args)], env=ENV, check=True)


def click(px, py, wait=1.0):
    x("mousemove", px, py, "click", 1)
    time.sleep(wait)


def typ(text):
    x("key", "ctrl+a")
    x("type", "--delay", 15, text)


class Rosy:
    def __init__(self, prefix):
        self.prefix = os.path.abspath(os.path.expanduser(prefix))
        self.drive_c = os.path.join(self.prefix, "drive_c")
        self.out = os.path.join(self.drive_c, "out")
        os.makedirs(self.out, exist_ok=True)
        self.wenv = dict(ENV, WINEDEBUG="-all", WINEPREFIX=self.prefix, WINEARCH="win32")

    def restart(self):
        subprocess.run(["wineserver", "-k"], env=self.wenv, check=False)
        time.sleep(2)
        subprocess.Popen(["wine", "cmd", "/c", "rosy.bat"], cwd=os.path.join(self.drive_c, "ROSY"),
                         env=self.wenv, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for _ in range(60):
            r = subprocess.run(["xdotool", "search", "--name", "^ROSY$"], env=ENV, capture_output=True, check=False)
            if r.stdout.strip():
                break
            time.sleep(1)
        time.sleep(8)
        click(662, 482, 2.0)  # "Would you like to input data from a file?" -> No

    def run_case(self, case_id, params, timeout=20):
        inp = os.path.join(self.out, f"in_{case_id}")
        res = os.path.join(self.out, f"res_{case_id}")
        rosy_io.write_case(os.path.join(HERE, "template.out"), inp, case_id, params)
        if os.path.exists(res):
            os.remove(res)
        click(822, 732, 2.0)  # Import File
        click(657, 514, 2.0)  # import from a ROSY output file -> OK
        click(664, 375, 0.3)  # file name field
        typ(f"C:\\out\\in_{case_id}")
        click(741, 560, 3.0)  # OK
        click(610, 279, 0.3)  # output file name
        typ(f"C:\\out\\res_{case_id}")
        click(729, 732, 0.5)  # Calculate
        t0 = time.time()
        while time.time() - t0 < timeout:
            if os.path.exists(res) and os.path.getsize(res) > 1000:
                time.sleep(0.5)
                return rosy_io.parse_output(res)
            time.sleep(0.5)
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases")
    ap.add_argument("results")
    ap.add_argument("--prefix", default="~/.wine-rosy")
    ap.add_argument("--restart", action="store_true", help="(re)start ROSY before the first case")
    a = ap.parse_args()
    rosy = Rosy(a.prefix)
    if a.restart:
        rosy.restart()
    results = {}
    with open(a.cases) as fh:
        cases = json.load(fh)
    for cid, params in cases.items():
        r = rosy.run_case(cid, params)
        results[cid] = r
        print(cid, "OK" if r else "FAIL", (r or {}).get("EU", {}).get("age"), flush=True)
        if r is None:  # usually a non-convergence dialog: restart to clear it
            rosy.restart()
        with open(a.results, "w") as fh:
            json.dump(results, fh, indent=1)


if __name__ == "__main__":
    main()
