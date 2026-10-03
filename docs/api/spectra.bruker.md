# EPR spectra: Bruker files

`eprdating.spectra.bruker`

Bruker cw-EPR files.

* BES3T (Xepr, E500/E580/EMXplus/EMXmicro…): `.DSC` descriptor + `.DTA`
  binary data.
* ESP / WinEPR (ESP300, EMX with WinEPR, ECS106…): `.par` parameters +
  `.spc` binary data.

The format rules follow EasySpin's `eprload` (Stoll & Schweiger, MIT
licence), against whose test files these readers are checked. Only 1-D
field sweeps and 2-D data whose second dimension is kept as separate rows
(e.g. repeated scans) are supported. The data are returned as stored: no
scaling by scans, gain or power is applied, but the parameters needed for it
are put into the `Spectrum` fields.

**Contents:** [`read_bes3t`](#read_bes3t), [`read_dsc`](#read_dsc), [`read_esp`](#read_esp), [`read_esp_par`](#read_esp_par)

### `read_bes3t`

`read_bes3t(path: str | Path, freq_GHz: float | None = None) -> Spectrum`

Bruker BES3T `.DSC`/`.DTA` cw spectrum.

### `read_dsc`

`read_dsc(path: str | Path) -> dict[str, str]`

Key/value pairs of a BES3T descriptor (comments and device blocks
skipped, continuation lines joined, quotes removed).

### `read_esp`

`read_esp(path: str | Path, freq_GHz: float | None = None) -> Spectrum`

Bruker ESP / WinEPR `.par`/`.spc` spectrum.

WinEPR files (`DOS` key present) hold little-endian float32; ESP files
big-endian int32 (EasySpin's rules).

### `read_esp_par`

`read_esp_par(path: str | Path) -> dict[str, str]`

`KEY value` pairs of an ESP/WinEPR `.par` file.
