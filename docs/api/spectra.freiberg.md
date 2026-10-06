# EPR spectra: Freiberg MS5000 files

`eprdating.spectra.freiberg`

Freiberg Instruments MS5000 spectra (ESRStudio `.xml` and `.csv`).

The `.xml` file holds everything: the measurement attributes (microwave
frequency `MwFreq` in GHz, Q factor, temperature, time stamp), the recipe
(sweep, `MicrowavePower` in mW, `Modulation` in mT, `Accumulations`)
and the curves as Base64 little-endian doubles, each sampled in time
(`XOffset + i * XSlope` s): the field `BField` (mT) and the signal
`MWAbsorption`. The spectrum is the signal against the field interpolated
at the signal's sampling times; this reproduces the `.csv` export of
ESRStudio (which adds a constant offset to the signal and has no frequency).

The `.csv` export has `key;value;description` recipe lines and, after a
`Meas` line, `BField [mT];MW_Absorption []` columns. When an `.xml`
with the same name sits next to it, the frequency is taken from there.

A measurement repeated with a recipe of several runs is saved as
`Name_1.xml` … `Name_10.xml` (each the average of `Accumulations`
sweeps) and, often, their average `Name_result.xml`. `read_ms5000`
on a folder reads the runs (not the `_result`) as the scans of one
spectrum.

**Contents:** [`is_ms5000`](#is_ms5000), [`read_ms5000`](#read_ms5000)

### `is_ms5000`

`is_ms5000(path: str | Path) -> bool`

Whether `path` is an ESRStudio `.xml` or `.csv` file.

### `read_ms5000`

`read_ms5000(path: str | Path, freq_GHz: float | None = None, include_result: bool = False) -> Spectrum`

Read an MS5000 (ESRStudio) `.xml` or `.csv` spectrum, or a folder of runs.

- `path`: an `.xml` file (preferred: it has the frequency), its `.csv` export, or a folder whose `Name_1.xml` … `Name_n.xml` runs become the scans of one spectrum (put on the field grid of the first). `Name_result.xml` files in the folder are skipped unless `include_result`.
- `freq_GHz`: overrides the frequency in the file (needed for a `.csv` without its `.xml`).

The power (mW) and modulation amplitude (mT) come from the recipe; there
is no receiver gain (`gain` is None).
