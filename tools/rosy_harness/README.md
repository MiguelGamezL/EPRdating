# ROSY harness

Batch-runs the original ROSY 2.0 (Windows, Java compiled with IBM HPJ) under
Wine to produce reference data for EPRdating. ROSY is **not** included: use
your own copy.

## Setup (Ubuntu)

```bash
sudo dpkg --add-architecture i386 && sudo apt update
sudo apt install wine wine32:i386 xvfb xdotool
export WINEPREFIX=~/.wine-rosy WINEARCH=win32
wineboot -i
cp -r /path/to/ROSY $WINEPREFIX/drive_c/ROSY
Xvfb :99 -screen 0 1280x900x24 &
```

## Run

```bash
cd tools/rosy_harness
DISPLAY=:99 python run_rosy.py cases.json results.json --prefix ~/.wine-rosy --restart
```

- `rosy_io.py` writes ROSY's "DATA FOR INPUT" block (the import format of
  `File > Import File > ROSY output file`) and parses ROSY output files.
  Field order is documented in `rosy_io.FIELDS`.
- `template.out` is a ROSY output file used as the skeleton for imports.
- Each case takes ~10 s. On a non-convergence dialog the harness restarts ROSY.

Findings from the first campaign: `tests/validation/ROSY_FINDINGS.md`.
