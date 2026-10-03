# Dose rates

`eprdating.dose_rate`

Environmental dose-rate building blocks.

* Infinite-matrix dose rates from U, Th, K concentrations (conversion factors).
* Water correction (Zimmerman 1971; Aitken 1985).
* Alpha efficiency (k-value).
* Cosmic dose rate (Prescott & Hutton 1994).

All dose rates are in Gy/ka.

**Contents:** [`ConversionFactors`](#conversionfactors), [`Sediment`](#sediment), [`available_factor_sets`](#available_factor_sets), [`conversion_factors`](#conversion_factors), [`cosmic_dose_rate`](#cosmic_dose_rate), [`cosmic_dose_rate_sea_level`](#cosmic_dose_rate_sea_level), [`geomagnetic_latitude`](#geomagnetic_latitude), [`matrix_dose_rates`](#matrix_dose_rates), [`samples_from`](#samples_from), [`u_disequilibrium_factor`](#u_disequilibrium_factor), [`u_series_split`](#u_series_split), [`water_correction`](#water_correction)

### `ConversionFactors`

*dataclass* `ConversionFactors(key: str, reference: str, table: dict[str, dict[str, Value]])`

Dose-rate conversion factors in Gy/ka per ppm (U, Th) or per % (K).

**Members**

- `get(self, nuclide: str, radiation: str) -> Value`

### `Sediment`

*dataclass* `Sediment(U: ValueLike = 0.0, Th: ValueLike = 0.0, K: ValueLike = 0.0, water: ValueLike = 0.0, U_ra226: ValueLike | None = None)`

Radionuclide content of a sediment or soil (U, Th in ppm, K in %).

`U_ra226`: 226Ra and its daughters as ppm of U in equilibrium, when it
differs from the 238U content `U` (e.g. from gamma spectrometry, 214Pb
and 214Bi lines vs 234Th and 234mPa). `None` means equilibrium.

**Members**

- `dose_rate(self, radiation: str, factors: ConversionFactors | None = None, values = None) -> float` — Wet dose rate for one radiation type. `values` overrides the inputs
(used by the Monte Carlo engine).

### `available_factor_sets`

`available_factor_sets() -> list`


### `conversion_factors`

`conversion_factors(key: str = 'guerin_2011') -> ConversionFactors`

Load a published set of conversion factors (see `available_factor_sets`).

### `cosmic_dose_rate`

`cosmic_dose_rate(depth_m: float, density: float, lat_deg: float, lon_deg: float, altitude_m: float, rel_sigma: float = 0.1) -> Value`

Cosmic dose rate after Prescott & Hutton (1994), Gy/ka.

Corrected for altitude and geomagnetic latitude with the F, J, H
factors of Prescott & Stefan (1982). A 10 % relative uncertainty is
assigned by default, as is common practice.

### `cosmic_dose_rate_sea_level`

`cosmic_dose_rate_sea_level(depth_m: float, density: float) -> float`

Cosmic dose rate at 55°N geomagnetic latitude and sea level (Gy/ka).

`depth_m` is the overburden thickness (m) and `density` its mean
density (g/cm³); their product is the shielding in hg/cm².

### `geomagnetic_latitude`

`geomagnetic_latitude(lat_deg: float, lon_deg: float) -> float`

Geomagnetic latitude (degrees) from geographic coordinates.

Dipole approximation with the pole at 78.3°N, 291°E, as used by
Prescott & Hutton (1994) and DRAC.

### `matrix_dose_rates`

`matrix_dose_rates(U: float = 0.0, Th: float = 0.0, K: float = 0.0, factors: ConversionFactors | None = None, U_ra226: float | None = None) -> dict[str, float]`

Dry infinite-matrix alpha/beta/gamma dose rates (Gy/ka).

`U` and `Th` in ppm, `K` in %. The Th chain is taken in secular
equilibrium; so is the U chain unless `U_ra226` (226Ra and daughters
expressed as ppm of U in equilibrium, as measured by gamma spectrometry
through 214Pb/214Bi) differs from `U` (238U).

### `samples_from`

`samples_from(obj, names, rng: np.random.Generator, n: int) -> dict[str, np.ndarray]`

Draw Gaussian samples for the listed attributes of `obj`.

### `u_disequilibrium_factor`

`u_disequilibrium_factor(radiation: str, U: float, U_ra226: float | None) -> float`

U-equivalent content giving the dose rate of a chain whose 226Ra and
daughters correspond to `U_ra226` ppm while the rest follows `U`.

### `u_series_split`

`u_series_split(radiation: str) -> tuple[float, float]`

Fractions of the natural-U dose rate emitted before and from 226Ra on.

"Before" is 238U, 234Th, 234Pa, 234U, 230Th and the 235U chain; "from
226Ra" is 226Ra, 222Rn and its daughters. Returns `(pre, post)`.

### `water_correction`

`water_correction(rate: float, water: float, radiation: str) -> float`

Correct a dry dose rate for water content.

`water` is the ratio mass of water / mass of dry material (e.g. 0.15).
