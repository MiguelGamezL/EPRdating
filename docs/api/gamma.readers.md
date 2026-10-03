# Gamma: reading spectra

`eprdating.gamma.readers`

Readers for gamma-ray spectrum files.

`read_gamma` picks the reader from the extension (and content):


```text
extension       format                                           reader
.Spe        ORTEC / IAEA ASCII ($DATA:, $MEAS_TIM:)  native
.Chn        ORTEC binary (MAESTRO)                           native
.n42/.xml   ANSI N42.42 (2006 and 2012 schemas)              native
.txt        # key: value header + channel columns        native
.csv/.txt   plain columns (counts, or channel/energy+counts) native
.cnf        Canberra Genie 2000                              becquerel
.iec/.spc   IEC 61455 / ORTEC binary SPC                     becquerel
```


The formats marked *becquerel* are read through the optional package
`becquerel <https://github.com/lbl-anp/becquerel>`_ (`pip install eprdating[gamma-formats]`),
whose parsers for these binary formats have been tested on real files.

Every reader returns a `GammaSpectrum` with counts,
live and real time and, when the file has one, the stored energy
calibration (used only as a first guess; see
`auto_calibrate`). Plain column files carry no times:
pass `live_time` (s).

**Contents:** [`read_chn`](#read_chn), [`read_columns`](#read_columns), [`read_gamma`](#read_gamma), [`read_n42`](#read_n42), [`read_spe`](#read_spe), [`read_with_becquerel`](#read_with_becquerel)

### `read_chn`

`read_chn(path: str | Path) -> GammaSpectrum`

ORTEC `.Chn` binary spectrum (MAESTRO).

Layout (little endian): int16 -1, int16 MCA number, int16 segment,
char[2] start seconds, int32 real time and int32 live time in 20 ms
ticks, char[8] date `DDMMMYY*` (`*` = '1' after 2000), char[4] time
`HHMM`, int16 first channel, int16 number of channels, then uint32
counts; an optional trailer (int16 -101 or -102) holds the energy
calibration as float32 (offset, slope[, quadratic]).

### `read_columns`

`read_columns(path: str | Path, live_time: float, real_time: float | None = None, counts_column: int = -1, energy_column: int | None = None, delimiter: str | None = None, skiprows: int = 0) -> GammaSpectrum`

Plain text/CSV spectrum: one column of counts, or several columns
with the counts in `counts_column`. If `energy_column` is given its
values are fitted with a linear calibration (first guess only).

### `read_gamma`

`read_gamma(path: str | Path, **kw) -> GammaSpectrum`

Read a gamma spectrum, choosing the reader from the file.

### `read_n42`

`read_n42(path: str | Path, index: int = 0) -> GammaSpectrum`

ANSI N42.42 XML spectrum (2006 or 2012). `index` selects the
spectrum when the file holds several.

### `read_spe`

`read_spe(path: str | Path) -> GammaSpectrum`

ORTEC / IAEA `.Spe` ASCII spectrum.

### `read_with_becquerel`

`read_with_becquerel(path: str | Path) -> GammaSpectrum`

Read any format supported by becquerel (CNF, SPC, IEC 61455, SPE).
