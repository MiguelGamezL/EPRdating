"""Drive DATA.EXE (Grün 2009) in DOSBox-X: start, keys, screenshots, print capture.

DATA is a DOS program; it is run unchanged as a black box. The working folder
(``DATA_DIR``, default ``./data_work``) must contain your copy of ``DATA.EXE``;
it is mounted as drive C:, cases are written to ``C:\\esr\\epr\\`` (DATA's
default path) and the printer port LPT1 is captured to ``out/lpt1.txt``.
"""
import os
import pathlib
import subprocess
import time

D = pathlib.Path(os.environ.get("DATA_DIR", "data_work")).resolve()
ENV = dict(os.environ, DISPLAY=os.environ.get("DISPLAY", ":99"), SDL_AUDIODRIVER="dummy")
LPT = D / "out" / "lpt1.txt"

CONF = """[sdl]
output=surface
autolock=false
[dosbox]
machine=svga_s3
memsize=16
[cpu]
cycles=max
[parallel]
parallel1=file file:{lpt}
[render]
scaler=none
[autoexec]
mount c {dir}
c:
DATA.EXE
"""


def sh(*a, **kw):
    return subprocess.run(a, env=ENV, capture_output=True, text=True, check=False, **kw)


def win():
    out = sh("xdotool", "search", "--name", "DOSBox").stdout.split()
    return out[0] if out else None


def start():
    """Start a fresh DOSBox-X session running DATA."""
    stop()
    if not (D / "DATA.EXE").exists():
        raise FileNotFoundError(f"DATA.EXE not found in {D} (set DATA_DIR)")
    (D / "esr" / "epr").mkdir(parents=True, exist_ok=True)
    LPT.parent.mkdir(exist_ok=True)
    if LPT.exists():
        LPT.unlink()
    (D / "dosbox.conf").write_text(CONF.format(lpt=LPT, dir=D))
    with open(D / "dosbox.log", "w") as log:
        subprocess.Popen(["dosbox-x", "-conf", str(D / "dosbox.conf"), "-nopromptfolder", "-fastlaunch"],
                         cwd=D, env=ENV, stdout=log, stderr=subprocess.STDOUT)
    for _ in range(50):
        time.sleep(0.2)
        if win():
            break
    time.sleep(5.0)


def stop():
    # match the process name exactly: `pkill -f dosbox-x` would also kill the
    # shell that runs this script if its command line mentions dosbox-x
    for pid in sh("pgrep", "-x", "dosbox-x").stdout.split():
        sh("kill", "-9", pid)
    time.sleep(0.5)


def key(*ks, pause=0.5):
    w = win()
    for k in ks:
        sh("xdotool", "key", "--window", w, k)
        time.sleep(pause)


def typ(text, pause=0.5):
    sh("xdotool", "type", "--window", win(), "--delay", "90", str(text))
    time.sleep(pause)


def field(v, pause=0.6):
    typ(v, 0.2)
    key("Return", pause=pause)


def shot(name):
    """Screenshot of the virtual display (for debugging a stuck case)."""
    sh("import", "-window", "root", str(D / f"{name}.png"))


def printed(wait=4.0):
    """Press F10 (print screen) and return the newly printed text block."""
    before = LPT.read_bytes() if LPT.exists() else b""
    key("F10", pause=0.3)
    t0 = time.time()
    while time.time() - t0 < wait:
        time.sleep(0.3)
        now = LPT.read_bytes() if LPT.exists() else b""
        if len(now) > len(before):
            time.sleep(1.0)
            return LPT.read_bytes()[len(before):].decode("cp437")
    return ""
