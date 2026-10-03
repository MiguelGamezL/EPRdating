# DATA harness

Batch-runs the original DATA program (Grün 2009, DOS) in DOSBox-X to produce
reference data for EPRdating. DATA is **not** included: use your own copy.

## Setup (Ubuntu)

```bash
sudo apt install dosbox-x xvfb xdotool imagemagick
mkdir data_work && cp /path/to/DATA.EXE data_work/
Xvfb :99 -screen 0 1280x900x24 &
```

## Run

```bash
cd tools/data_harness
DATA_DIR=data_work DISPLAY=:99 python run_data.py \
    ../../tests/validation/data_reference/cases.json results.json
```

- Each case is written as an `.EPR` file in `C:\esr\epr\` (DATA's default
  path), loaded with <F5>, calculated with <F1> and printed with <F10>;
  DOSBox-X captures the printer port LPT1 to `data_work/out/lpt1.txt`.
- `run_data.ORDER` documents the 19 fields of an `.EPR` file. The "beta
  only" switch (`Y`) takes the external gamma dose rate (µGy/a) directly.
- File names are DOS 8.3 (short hashed names are used).
- A fresh DOS session is started for every case (~15 s each); a failed case
  is retried once and a screenshot is left in `data_work/`.
- `parse` reads the RESULTS block (EU and LU dose rates in µGy/a, ages in
  ka) and DATA's beta-correction table. The raw printouts of the reference
  run are in `tests/validation/data_reference/raw/`.

## U-series/ESR (US-ESR, CS-US)

```bash
DATA_DIR=data_work DISPLAY=:99 python run_useries.py \
    ../../tests/validation/data_reference/useries_cases.json useries_results.json
```

- <F6> loads the `.EPR` file into the U-series screen. The ratios are typed
  in, because they cannot be stored in the file. A field is edited with
  <Enter> value <Enter> <Enter> error <Enter>; `run_useries.TISSUE_ROWS`
  lists the screen rows.
- US-ESR is iterative (a counter runs on the screen), so the harness waits
  until the screen stops changing before printing (~50 s per case).
- The echoed inputs on the printout are compared with the case. A mismatch
  is retried.
- Avoid inputs that crash DATA: a 100 % error on De (overflow, line 116) or
  no U in the dentine (line 114).

Findings: `tests/validation/DATA_FINDINGS.md`.
