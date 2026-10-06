# EPR spectra: reading

`eprdating.spectra.io`

Reading cw-EPR spectra from disk.

`read_epr` reads any supported file: Bruker BES3T and ESP/WinEPR
(`bruker`), Freiberg MS5000 (`freiberg`),
the `.dat`/`.par` pairs described below, and plain field/signal columns. All return a `Spectrum` with
the field in mT.

**`.dat`/`.par` format**

`.dat` + `.par` ASCII pairs (one spectrum, possibly several scans). The
`.par` file holds `KEY : value` lines:

```
N : 512            points per scan
CF : 3350.0        centre field set by the operator (G)
CF_ : 3354.0       actual centre field (G)
SW : 499.1453      sweep width (G)
SW_ : 500.0        second sweep-width value (kept, not used)
Nscans : 4
Freq : 9.43        microwave frequency (GHz)
TC, MA, OF, PH, RG, CT   time constant, modulation amplitude, offset,
                         phase, receiver gain, conversion time (raw units)
```

and the `.dat` file has five whitespace-separated columns: running index,
point index within the scan, field (G), signal, and a flag. Scans follow one
another; a line with `NaN` signal may separate them, and an interrupted
acquisition leaves a trailing incomplete scan, which is dropped.

The field column is CF ± SW/2, i.e. it is built on the *set* centre field.
`read_dat` moves it by `CF_ - CF` so that the axis is the actual
field (`actual_field=False` keeps the raw axis). For the M18 series this
correction (+0.4 mT) brings the CO2- signal to within 0.05 mT of the
simulated one (literature g-values, nominal 9.43 GHz). Fields are in mT, as
everywhere in `eprdating`.

**Contents:** [`Spectrum`](#spectrum), [`read_columns`](#read_columns), [`read_dat`](#read_dat), [`read_epr`](#read_epr), [`read_par`](#read_par), [`read_series`](#read_series)

### `Spectrum`

*dataclass* `Spectrum(B: np.ndarray, scans: np.ndarray, name: str = '', params: dict[str, object] = <factory>, freq_GHz: float | None = None, power_mW: float | None = None, dropped_points: int = 0, gain: float | None = None, mod_amp_mT: float | None = None, time_constant_ms: float | None = None)`

A cw-EPR (first-derivative) spectrum with its individual scans.

- `B`: field grid in mT (increasing).
- `scans`: array (n_scans, n_points) of raw signal.
- `params`: acquisition parameters as read from the file (raw units).
- `freq_GHz`, `power_mW`: microwave frequency and power (for the .dat/.par format the power comes from the file name, e.g. `M18_3_19mW_4SCAN.dat`).
- `gain`: receiver gain as a linear factor (None if unknown).
- `mod_amp_mT`, `time_constant_ms`: field modulation (peak to peak) and lock-in time constant, when the file gives them in known units.

**Members**

- `n_scans` *(property)*
- `window(self, lo_mT: float, hi_mT: float) -> Spectrum` — Copy restricted to `lo_mT <= B <= hi_mT`.
- `y` *(property)* — Average of the scans.

### `read_columns`

`read_columns(path: str | Path, field_unit: str = 'G', freq_GHz: float | None = None, power_mW: float | None = None, field_column: int = 0, signal_columns = None, delimiter: str | None = None) -> Spectrum`

Plain text/CSV spectrum: a field column and one or more signal
columns (each taken as a scan). Lines starting with `#` are skipped.

### `read_dat`

`read_dat(path: str | Path, par: str | Path | None = None, power_mW: float | None = None, field_unit: str = 'G', actual_field: bool = True) -> Spectrum`

Read a `.dat` spectrum (and its `.par` file if present).

`par` defaults to the file with the same stem. `power_mW` defaults to
the number before `mW` in the file name. `field_unit` is the unit of
the field column ("G" or "mT"). With `actual_field` the axis is shifted
by `CF_ - CF` (actual minus set centre field) when both are present.

### `read_epr`

`read_epr(path: str | Path, **kw) -> Spectrum`

Read a cw-EPR spectrum, choosing the reader from the file.


```text
.DSC / .DTA     Bruker BES3T (Xepr)
.par + .spc     Bruker ESP / WinEPR
.dat + .par     KEY : value + five-column ASCII (see above)
.xml / .csv     Freiberg MS5000 (ESRStudio); a folder of MS5000
                        runs is read as the scans of one spectrum
.txt / .csv     field and signal columns
```

### `read_par`

`read_par(path: str | Path) -> dict[str, float]`

Parse a `KEY : value` parameter file.

### `read_series`

`read_series(paths: Iterable[str | Path], **kw) -> list[Spectrum]`

Read several spectra (e.g. a dose series), any supported format.
