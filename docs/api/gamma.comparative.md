# Gamma: comparative method

`eprdating.gamma.comparative`

U, Th and K contents by the comparative (relative) method.

The sample and reference materials of known content are measured in the same
container and geometry. For every gamma line the background-corrected count
rate per unit mass of the sample is divided by that of the reference:

    C_sample = C_ref · (R_s / m_s) / (R_ref / m_ref),   R = net area / live time

Efficiency, emission probability and coincidence summing cancel; what does
not cancel is a difference in fill height or density (self-absorption),
which matters most below ~200 keV and when the masses differ.

**Line groups**

`K`      40K, 1460.8 keV.
`Ra226`  214Pb and 214Bi lines: 226Ra through its radon daughters. Gives
           the U content if the chain is in equilibrium (sample sealed for
           ~3–4 weeks so that 222Rn and its daughters grow back).
`U238`   234Th (63.3 keV) and 234mPa (1001.0 keV): the top of the chain,
           i.e. 238U itself. Weak; compared with `Ra226` it tests the
           equilibrium the dose-rate model assumes.
`Th232`  228Ac, 212Pb and 208Tl lines.

Peak areas come from Gaussians at the calibrated positions and widths of each
spectrum (calibrated separately, so gain drift between days is absorbed) on a
linear background; only amplitudes and background are fitted, which keeps
weak peaks unbiased. Lines of one group share the reference uncertainty, so
it is added after averaging.

**Reference values (dry mass, 1σ)**

IAEA-RGU-1: 238U 4941 ± 99 Bq/kg → U = 400 ± 8 µg/g (certificate);
IAEA-RGTh-1: Th = 800 ± 16 µg/g (IAEA/RL/148);
IAEA-RGK-1: K = 448 ± 3 g/kg (certificate, 2016).

**Contents:** [`BQ_PER_KG`](#bq_per_kg), [`GROUP_ELEMENT`](#group_element), [`IAEA_RGK_1`](#iaea_rgk_1), [`IAEA_RGTH_1`](#iaea_rgth_1), [`IAEA_RGU_1`](#iaea_rgu_1), [`LINES`](#lines), [`GammaResult`](#gammaresult), [`GroupResult`](#groupresult), [`Line`](#line), [`LineResult`](#lineresult), [`Reference`](#reference), [`analyse`](#analyse), [`analyse_files`](#analyse_files), [`calibrate_natural`](#calibrate_natural), [`line_area`](#line_area), [`net_rate`](#net_rate), [`water_content`](#water_content)

### `BQ_PER_KG`

```python
BQ_PER_KG = {'U': 12.35, 'Th': 4.057, 'K': 316.5}
```

### `GROUP_ELEMENT`

```python
GROUP_ELEMENT = {'K': 'K', 'Ra226': 'U', 'U238': 'U', 'Th232': 'Th'}
```

### `IAEA_RGK_1`

```python
IAEA_RGK_1 = Reference(name='IAEA-RGK-1', content={'K': (44.8, 0.3)}, mass_g=500.0)
```

### `IAEA_RGTH_1`

```python
IAEA_RGTH_1 = Reference(name='IAEA-RGTh-1', content={'Th': (800.0, 16.0)}, mass_g=500.0)
```

### `IAEA_RGU_1`

```python
IAEA_RGU_1 = Reference(name='IAEA-RGU-1', content={'U': (400.0, 8.0)}, mass_g=500.0)
```

### `LINES`

```python
LINES = (Line(energy=1460.82, emitter='40K', group='K', neighbours=(), window=None), Line(energy=295.224, emitter='214Pb', group='Ra226', neighbours=(300.087,), window=None), Line(energy=351.932, emitter='214Pb', group='Ra226', neighbours=(), window=None), Line(energy=609.312, emitter='214Bi', group='Ra226' …
```

### `GammaResult`

*dataclass* `GammaResult(sample: str, mass_g: float, groups: dict[str, GroupResult], activity_ratio_ra226_u238: Value | None = None)`

GammaResult(sample: 'str', mass_g: 'float', groups: 'dict[str, GroupResult]', activity_ratio_ra226_u238: 'Value | None' = None)

**Members**

- `K` *(property)*
- `Th` *(property)*
- `U` *(property)* — U (µg/g) from the 226Ra daughters (equilibrium assumed).
- `sediment(self, water: ValueLike = 0.0, equilibrium: bool = True) -> Sediment` — A `Sediment` with these (dry-mass) contents.
- `summary(self) -> str`

### `GroupResult`

*dataclass* `GroupResult(group: str, element: str, content: Value, ratio: float, ratio_sigma: float, chi2_red: float, lines: list[LineResult] = <factory>)`

GroupResult(group: 'str', element: 'str', content: 'Value', ratio: 'float', ratio_sigma: 'float', chi2_red: 'float', lines: 'list[LineResult]' = <factory>)

### `Line`

*dataclass* `Line(energy: float, emitter: str, group: str, neighbours: tuple[float, ...] = (), window: float | None = None)`

Line(energy: 'float', emitter: 'str', group: 'str', neighbours: 'tuple[float, ...]' = (), window: 'float | None' = None)

### `LineResult`

*dataclass* `LineResult(line: Line, sample_rate: float, sample_rate_sigma: float, ref_rate: float, ref_rate_sigma: float, ratio: float, ratio_sigma: float, content: float, content_sigma: float)`

LineResult(line: 'Line', sample_rate: 'float', sample_rate_sigma: 'float', ref_rate: 'float', ref_rate_sigma: 'float', ratio: 'float', ratio_sigma: 'float', content: 'float', content_sigma: 'float')

### `Reference`

*dataclass* `Reference(name: str, content: Mapping[str, ValueLike], mass_g: float = 500.0)`

A reference material: content per element (U, Th in µg/g, K in %).

### `analyse`

`analyse(sample: GammaSpectrum, mass_g: float, references: Mapping[str, tuple[GammaSpectrum, Reference]], background: GammaSpectrum | None = None, lines: Sequence[Line] = (Line(energy=1460.82, emitter='40K', group='K', neighbours=(), window=None), Line(energy=295.224, emitter='214Pb', group='Ra226', neighbours=(300.087,), window=None), Line(energy=351.932, emitter='214Pb', group='Ra226', neighbours=(), window=None), Line(energy=609.312, emitter='214Bi', group='Ra226', neighbours=(), window=None), Line(energy=1120.287, emitter='214Bi', group='Ra226', neighbours=(), window=None), Line(energy=1764.494, emitter='214Bi', group='Ra226', neighbours=(), window=None), Line(energy=63.29, emitter='234Th', group='U238', neighbours=(), window=None), Line(energy=1001.03, emitter='234mPa', group='U238', neighbours=(), window=None), Line(energy=238.632, emitter='212Pb', group='Th232', neighbours=(241.997,), window=None), Line(energy=338.32, emitter='228Ac', group='Th232', neighbours=(), window=None), Line(energy=583.187, emitter='208Tl', group='Th232', neighbours=(), window=None), Line(energy=911.204, emitter='228Ac', group='Th232', neighbours=(), window=None), Line(energy=968.971, emitter='228Ac', group='Th232', neighbours=(964.766,), window=None), Line(energy=2614.511, emitter='208Tl', group='Th232', neighbours=(), window=None)), groups: Sequence[str] | None = None) -> GammaResult`

Contents of the sample by the comparative method.

`references` maps an element ("U", "Th", "K") to (spectrum, Reference).
All spectra must be calibrated (see `calibrate_natural`). Lines of
groups whose element has no reference are skipped.

### `analyse_files`

`analyse_files(sample: str | Path, mass_g: float, references: Mapping[str, tuple[str | Path, Reference]], background: str | Path | None = None, lines: Sequence[Line] = (Line(energy=1460.82, emitter='40K', group='K', neighbours=(), window=None), Line(energy=295.224, emitter='214Pb', group='Ra226', neighbours=(300.087,), window=None), Line(energy=351.932, emitter='214Pb', group='Ra226', neighbours=(), window=None), Line(energy=609.312, emitter='214Bi', group='Ra226', neighbours=(), window=None), Line(energy=1120.287, emitter='214Bi', group='Ra226', neighbours=(), window=None), Line(energy=1764.494, emitter='214Bi', group='Ra226', neighbours=(), window=None), Line(energy=63.29, emitter='234Th', group='U238', neighbours=(), window=None), Line(energy=1001.03, emitter='234mPa', group='U238', neighbours=(), window=None), Line(energy=238.632, emitter='212Pb', group='Th232', neighbours=(241.997,), window=None), Line(energy=338.32, emitter='228Ac', group='Th232', neighbours=(), window=None), Line(energy=583.187, emitter='208Tl', group='Th232', neighbours=(), window=None), Line(energy=911.204, emitter='228Ac', group='Th232', neighbours=(), window=None), Line(energy=968.971, emitter='228Ac', group='Th232', neighbours=(964.766,), window=None), Line(energy=2614.511, emitter='208Tl', group='Th232', neighbours=(), window=None)), calibration_lines: Sequence[float] = (238.632, 351.932, 583.187, 609.312, 911.204, 1460.82, 1764.494, 2614.511), **read_kw) -> GammaResult`

One call from files to contents.

Every file is read with `read_gamma` (any
supported format) and calibrated on its own natural lines with
`auto_calibrate` (no first guess needed), then
`analyse` is applied. Example:

```
r = analyse_files("soil.Spe", 600.0,
                  {"U": ("RGU1.Spe", IAEA_RGU_1), "Th": ("RGTh1.Spe", IAEA_RGTH_1),
                   "K": ("RGK1.Spe", IAEA_RGK_1)},
                  background="empty.Spe")
```

### `calibrate_natural`

`calibrate_natural(spec: GammaSpectrum, guess: Calibration, lines: Sequence[float] = (238.632, 351.932, 583.187, 609.312, 911.204, 1460.82, 1764.494, 2614.511), min_significance: float = 8.0) -> Calibration`

Calibrate an environmental spectrum on its own natural lines (absorbs
gain drift between measurements).

### `line_area`

`line_area(spec: GammaSpectrum, line: Line, cal: Calibration | None = None) -> tuple[float, float]`

Area (counts) of `line` and its 1σ, with neighbours fitted jointly.

Gaussian positions and widths come from the spectrum calibration; the
amplitudes and a linear background are linear parameters, fitted with
Poisson weights from the model (two iterations).

### `net_rate`

`net_rate(spec: GammaSpectrum, line: Line, background: GammaSpectrum | None) -> tuple[float, float]`


### `water_content`

`water_content(wet_mass_g: ValueLike, dry_mass_g: ValueLike) -> Value`

Water content as mass of water / dry mass, the convention of
`Sediment`.
