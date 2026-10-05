# Top level

`eprdating`

Everything most scripts need, importable as `from eprdating import ...`.

**Contents:** [`AgeMC`](#agemc), [`AgeResult`](#ageresult), [`BetaGeometry`](#betageometry), [`DelayedUptake`](#delayeduptake), [`DoseRateComponent`](#doseratecomponent), [`DoseResponseResult`](#doseresponseresult), [`EarlyUptake`](#earlyuptake), [`History`](#history), [`LinearUptake`](#linearuptake), [`Material`](#material), [`Sediment`](#sediment), [`ToothLayers`](#toothlayers), [`ToothSample`](#toothsample), [`USESRSample`](#usesrsample), [`USModel`](#usmodel), [`USeries`](#useries), [`UseriesData`](#useriesdata), [`Value`](#value), [`available_factor_sets`](#available_factor_sets), [`bootstrap_De`](#bootstrap_de), [`compound`](#compound), [`conversion_factors`](#conversion_factors), [`cosmic_dose_rate`](#cosmic_dose_rate), [`cosmic_history`](#cosmic_history), [`dentine_material`](#dentine_material), [`fit_dose_response`](#fit_dose_response), [`matrix_dose_rates`](#matrix_dose_rates), [`mixture`](#mixture), [`sediment_material`](#sediment_material), [`solve_age`](#solve_age), [`water_correction`](#water_correction)

### `AgeMC`

*dataclass* `AgeMC(samples: np.ndarray, nominal: AgeResult, n_failed: int = 0)`

Monte Carlo age distribution.

**Members**

- `interval(self, level: float = 0.68) -> tuple`
- `mean` *(property)*
- `std` *(property)*
- `summary(self) -> str`

### `AgeResult`

*dataclass* `AgeResult(age: float, De: float, components: dict[str, float], accumulated: dict[str, float])`

AgeResult(age: 'float', De: 'float', components: 'dict[str, float]', accumulated: 'dict[str, float]')

**Members**

- `mean_dose_rate` *(property)* — Time-averaged dose rate De / T (Gy/ka).
- `summary(self) -> str`

### `BetaGeometry`

*dataclass* `BetaGeometry(internal: ValueLike, dentine: ValueLike, external: ValueLike)`

Fractions of infinite-matrix beta dose rate reaching the dated enamel.

- `internal`: enamel self-dose (U inside the enamel layer).
- `dentine`: from U in the adjacent dentine.
- `external`: from the sediment (or cement) on the outer side.

**Members**

- `values(self) -> dict`

### `DelayedUptake`

*dataclass* `DelayedUptake(t_uptake: float)`

All U taken up at once `t_uptake` ka before present (none before).

Used by the CSUS-ESR model (Grün 2000), where `t_uptake` is the
closed-system U-series age of the tissue. For a sample older than
`t_uptake` the tissue delivers dose only during its last `t_uptake`.

**Members**

- `accumulated(self, T: float, G = None) -> float`
- `fraction(self, t, T: float)` — U(t)/U_m at time `t` after burial for a sample of age `T`.

### `DoseRateComponent`

*dataclass* `DoseRateComponent(name: str, rate: float, uptake: USModel | None = None, G: Callable[[float], float] | None = None, profile: tuple[Sequence[float], Sequence[float]] | None = None)`

One contribution to the total dose rate (Gy/ka).

- `rate`: present-day dose rate the source would deliver in secular equilibrium with its present U content.
- `uptake`: None for a constant source, else a `USModel`.
- `G`: time-integrated activity ratio of incorporated U (`tau -> ∫ D/D_eq`); None means secular equilibrium.
- `profile`: for a constant source whose rate changed in the past, `(breaks, rates)`: piecewise-constant rates in ka before present (see `History`); `rate` is then `rates[0]`.

**Members**

- `accumulated(self, T: float) -> float` — Dose (Gy) delivered by this component over `T` ka.

### `DoseResponseResult`

*dataclass* `DoseResponseResult(model: str, params: dict[str, float], errors: dict[str, float], covariance: np.ndarray, dose: np.ndarray, intensity: np.ndarray, sigma: np.ndarray | None, chi2_red: float, r2: float, De_min: float = 0.0)`

Outcome of `fit_dose_response`.

**Members**

- `De` *(property)* — Equivalent dose (same units as the input doses, usually Gy).
- `De_samples(self, n: int, rng: np.random.Generator | None = None) -> np.ndarray` — Draw De from the fit's multivariate-normal parameter distribution.
- `De_sigma` *(property)*
- `predict(self, dose) -> np.ndarray` — Evaluate the fitted curve at the given added doses.
- `residuals(self) -> np.ndarray`
- `summary(self) -> str`

### `EarlyUptake`

`EarlyUptake() -> USModel`


### `History`

*dataclass* `History(values: Sequence[ValueLike], breaks: Sequence[float] = ())`

Piecewise-constant parameter history, in ka before present.

`values[0]` holds from today back to `breaks[0]`, `values[i]`
between `breaks[i-1]` and `breaks[i]`, and the last value before
the last break.

**Members**

- `at(self, t: float) -> float` — Nominal value at `t` ka before present.
- `mean(self, T: float) -> float` — Time-weighted nominal mean over the last `T` ka.
- `nominal(self) -> list[float]`
- `sample(self, rng: np.random.Generator, n: int) -> list[np.ndarray]` — `n` draws of every segment value, sampled independently.

### `LinearUptake`

`LinearUptake() -> USModel`


### `Material`

*dataclass* `Material(name: str, fractions: dict)`

Mass fractions of elements (normalised on construction).

**Members**

- `coefficients(self, E: float, scatter_factor: float = 1.0) -> tuple[float, float]` — (μa, μs) in cm²/g at energy E (MeV).
- `key` *(property)* — Identity by composition (two materials with the same name but
different compositions never share cached coefficients).
- `with_water(self, water: float) -> Material` — Add `water` grams of water per gram of dry material.

### `Sediment`

*dataclass* `Sediment(U: ValueLike = 0.0, Th: ValueLike = 0.0, K: ValueLike = 0.0, water: ValueLike | History = 0.0, U_ra226: ValueLike | None = None)`

Radionuclide content of a sediment or soil (U, Th in ppm, K in %).

`U_ra226`: 226Ra and its daughters as ppm of U in equilibrium, when it
differs from the 238U content `U` (e.g. from gamma spectrometry, 214Pb
and 214Bi lines vs 234Th and 234mPa). `None` means equilibrium.
`water` (mass of water / dry mass) may be a
`History` when it changed during burial; its
first value is the present-day one.

**Members**

- `dose_rate(self, radiation: str, factors: ConversionFactors | None = None, values = None) -> float` — Wet dose rate for one radiation type, with the present-day water.
`values` overrides the inputs (used by the Monte Carlo engine).

### `ToothLayers`

*dataclass* `ToothLayers(enamel_um: ValueLike, dentine_um: ValueLike = 2000.0, cementum_um: ValueLike = 0.0, strip_outer_um: ValueLike = 0.0, strip_inner_um: ValueLike = 0.0, enamel_density: ValueLike = 3.0, dentine_density: ValueLike = 2.82, cementum_density: ValueLike = 2.54, sediment_density: ValueLike = 2.0, sediment_water: float = 0.0, dentine_water: float = 0.0, cementum_water: float = 0.0, enamel: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397}), dentine: Material = Material(name='dentine', fractions={'Ca': 0.2792575058679254, 'P': 0.12949331783152535, 'O': 0.4306537832343316, 'H': 0.026595393066217624, 'C': 0.10000000000000002, 'N': 0.03400000000000001}), cementum: Material = Material(name='dentine', fractions={'Ca': 0.2792575058679254, 'P': 0.12949331783152535, 'O': 0.4306537832343316, 'H': 0.026595393066217624, 'C': 0.10000000000000002, 'N': 0.03400000000000001}), sediment: Material = Material(name='silica', fractions={'Si': 0.46743671254764246, 'O': 0.5325632874523576}), scatter_factor: float = 1.0)`

Planar tooth geometry for one-group beta attenuation.

Thicknesses in µm, densities in g/cm³, water in g per g of dry material.
- `Order`: outer sediment | cementum | enamel | dentine | inner sediment.
Geometric inputs accept a number, a `(value, sigma)` tuple or a
`Value`; the uncertainties are sampled by
`age_mc`.

**Members**

- `at(self, **values) -> ToothLayers` — Copy with the given fields set (floats), e.g. one Monte Carlo draw.
- `chain_fraction(self, source: str, chain: str) -> float` — Energy-weighted fraction for a whole decay chain or segment.
- `fraction(self, source: str, E: float) -> float` — Mean dose in the dated enamel per unit infinite-matrix dose of
`source` ('enamel', 'dentine', 'cementum', 'sediment' — both sides —,
'sediment_outer', 'sediment_inner') for energy E.
- `nominal_values(self) -> dict`

### `ToothSample`

*dataclass* `ToothSample(De: ValueLike, enamel_U: ValueLike, dentine_U: ValueLike, sediment: Sediment, beta: BetaGeometry | ToothLayers, cosmic: ValueLike | History, gamma: ValueLike | History | None = None, k_alpha: ValueLike = 0.13 ± 0.02, dentine_water: ValueLike = 0.0, uptake_enamel: USModel = <factory>, uptake_dentine: USModel = <factory>, u234_u238_enamel: ValueLike = 1.0, u234_u238_dentine: ValueLike = 1.0, radon_loss_enamel: ValueLike = 0.0, radon_loss_dentine: ValueLike = 0.0, enamel_water: ValueLike = 0.0, cementum_U: ValueLike = 0.0, cementum_water: ValueLike = 0.0, uptake_cementum: USModel = <factory>, u234_u238_cementum: ValueLike = 1.0, radon_loss_cementum: ValueLike = 0.0, ingrowth: bool = True, partition: dict | None = None, factors: str = 'guerin_2011', sample_geometry: bool = True, alpha_efficiency: str = 'constant', u234_u238_is: str = 'present', alpha_eref: float = 5.3, beta_by_segment: bool = True)`

ESR dating of tooth enamel with the classical component model.

- `Concentrations`: U in ppm. Dose in Gy, dose rates in Gy/ka, ages in ka.

**Parameters**

- `De`: equivalent dose of the enamel.
- `enamel_U`, `dentine_U`: present-day U in each tissue.
- `sediment`: U, Th, K and water content of the surrounding sediment. The water may be a `History` (wetter or drier periods); the sediment beta and gamma dose rates then follow it, and the present-day value (its first segment) is used elsewhere.
- `beta`: either fixed geometry factors (`BetaGeometry`) or a layered geometry (`ToothLayers`), in which case beta attenuation is computed with one-group theory for each emitter and U-series segment (as ROSY does). Uncertainties on the layer thicknesses, stripping and densities are sampled in `age_mc` together with the water contents.
- `gamma`: external gamma dose rate. If None, computed from `sediment` as an infinite matrix (use in-situ measurements when available). A measured value is taken as today's; with a water history it is rescaled to the water of each period. A `History` is used as given.
- `cosmic`: cosmic dose rate (see `cosmic_dose_rate`), or a `History` of it, e.g. from a burial depth history with `cosmic_history`.
- `k_alpha`: alpha efficiency of enamel. With `alpha_efficiency="energy"` it is the value at `alpha_eref` MeV.
- `alpha_efficiency`: `"constant"` (default, as in DATA) or `"energy"` (as in ROSY): k varies with alpha energy as R(E)/E, so each U-series segment gets its own efficiency (see `alpha`).
- `alpha_eref`: reference alpha energy for `k_alpha` in MeV (ROSY: 5.3).
- `dentine_water`: water content used for the dentine beta contribution.
- `enamel_water`: water content of the enamel (corrects the internal alpha and beta dose rates, Zimmerman coefficients).
cementum_U, cementum_water, uptake_cementum, u234_u238_cementum,
- `radon_loss_cementum`: the same for a cementum layer on the outer side of the enamel (needs a `ToothLayers` geometry with `cementum_um > 0`).
- `uptake_enamel`, `uptake_dentine`: uptake models (EU, LU, US with p).
- `u234_u238_enamel`, `u234_u238_dentine`: 234U/238U activity ratio of each tissue, measured today (`u234_u238_is="present"`, default) or of the incoming uranium (`"initial"`, used by `usesr`).
- `radon_loss_enamel`, `radon_loss_dentine`: fraction of 222Rn escaping each tissue (0 to 1).
- `ingrowth`: model U-series daughter ingrowth after uptake (default). With `False` the tissues are taken in secular equilibrium (warns).
- `partition`: optional custom U-series segment table (see `series`).
- `factors`: name of the conversion-factor set.
- `sample_geometry`: with a `ToothLayers` geometry, recompute the one-group factors for every Monte Carlo draw (default). With `False` they are kept at their nominal values.
- `beta_by_segment`: with a `ToothLayers` geometry, attenuate the beta dose of each U-series segment with its own factor (default, as ROSY). `False` applies one factor for the whole chain to the dose with ingrowth, as the DATA program (Grün 2009) does; for young teeth this lowers the dentine beta dose by up to ~40 % because the hard 234mPa betas dominate before 226Ra grows in.

**Members**

- `age(self) -> AgeResult` — Nominal age with the central value of every input.
- `age_mc(self, n: int = 2000, seed: int | None = None) -> AgeMC` — Monte Carlo age: every input is sampled from its Gaussian uncertainty.
- `components(self, v: dict[str, float] | None = None, uptake_enamel: USModel | None = None, uptake_dentine: USModel | None = None, uptake_cementum: USModel | None = None) -> list[DoseRateComponent]` — Dose-rate components for the input values `v` (nominal if None).

### `USESRSample`

*dataclass* `USESRSample(tooth: ToothSample, enamel: UseriesData | None = None, dentine: UseriesData | None = None, cementum: UseriesData | None = None)`

A tooth for US-ESR dating.

`tooth` carries everything except the uptake history (its uptake models
and 234U/238U fields are ignored for tissues with U-series data).

**Members**

- `age(self, model: str = 'US') -> USESRResult` — Nominal age. `model="US"`: US-ESR (uptake parameter p of each
tissue solved with the age, Grün et al. 1988). `model="CSUS"`:
CSUS-ESR (Grün 2000), U taken up at once at each tissue's
closed-system U-series age; comparing both shows how much the age
depends on the uptake model.
- `age_mc(self, n: int = 1000, seed: int | None = None, marginal_below: float = 0.8, model: str = 'US') -> USESRMC` — Monte Carlo over all inputs, U-series ratios included.

### `USModel`

*dataclass* `USModel(p: float = 0.0)`

Uptake following `(t/T)**(p+1)`; `p=-1` is EU, `p=0` is LU.

**Members**

- `accumulated(self, T: float, G = None) -> float` — Dose accumulated over `T` per unit present-day dose rate.
- `fraction(self, t, T: float)` — U(t)/U_m at time `t` after burial for a sample of age `T`.
- `is_early` *(property)*

### `USeries`

*dataclass* `USeries(ratio: float = 1.0, ratio_is: str = 'present', radon_loss: float = 0.0, partition: dict[str, dict[str, float]] | None = None)`

Time-integrated dose of U incorporated in a tissue, per unit of its
equilibrium dose rate, accounting for daughter ingrowth.

**Parameters**

- `ratio`: 234U/238U activity ratio.
- `ratio_is`: `"present"` (measured today, default) or `"initial"`.
- `radon_loss`: fraction of 222Rn escaping the tissue (0 to 1).
- `partition`: segment fractions; default is the bundled table.

`G(radiation)(tau)` returns ∫_0^tau (D(s)/D_eq) ds for U that has been
in the tissue for `tau` ka.

**Members**

- `G(self, radiation: str, weights: dict[str, float] | None = None)` — Time-integrated dose function for one radiation type.
- `initial_ratio(self, tau: float) -> float` — Initial 234U/238U of a U parcel incorporated `tau` ka ago.

### `UseriesData`

*dataclass* `UseriesData(th230_u234: ValueLike, u234_u238: ValueLike)`

Measured U-series activity ratios of one dental tissue.

### `Value`

*dataclass* `Value(value: float, sigma: float = 0.0)`

A measured quantity with a 1-sigma (Gaussian) uncertainty.

**Members**

- `rel` *(property)* — Relative uncertainty sigma/|value| (`inf` if value is 0).
- `sample(self, rng: np.random.Generator, n: int) -> np.ndarray` — Draw `n` Gaussian samples (constant array if `sigma == 0`).

### `available_factor_sets`

`available_factor_sets() -> list`


### `bootstrap_De`

`bootstrap_De(result: DoseResponseResult, n: int = 1000, seed: int | None = None) -> np.ndarray`

Residual bootstrap of De. Returns the array of resampled De values.

### `compound`

`compound(name: str, formula: dict) -> Material`

A material from its chemical formula, e.g. `compound("calcite", {"Ca": 1, "C": 1, "O": 3})`.

Elements available: `ELEMENTS` (H, C, N, O, Na, Mg, Al, Si, P, K, Ca, Fe).

### `conversion_factors`

`conversion_factors(key: str = 'guerin_2011') -> ConversionFactors`

Load a published set of conversion factors (see `available_factor_sets`).

### `cosmic_dose_rate`

`cosmic_dose_rate(depth_m: float, density: float, lat_deg: float, lon_deg: float, altitude_m: float, rel_sigma: float = 0.1) -> Value`

Cosmic dose rate after Prescott & Hutton (1994), Gy/ka.

Corrected for altitude and geomagnetic latitude with the F, J, H
factors of Prescott & Stefan (1982). A 10 % relative uncertainty is
assigned by default, as is common practice.

### `cosmic_history`

`cosmic_history(depth_m: History, density: float, lat_deg: float, lon_deg: float, altitude_m: float, rel_sigma: float = 0.1) -> History`

Cosmic dose-rate history from a burial-depth history (m, ka before present).

Each segment gets `cosmic_dose_rate` at its depth, with the
`rel_sigma` uncertainty combined in quadrature with that of the depth
(propagated through the local slope of the depth curve). A gradual
burial is approximated by several short segments.

### `dentine_material`

`dentine_material(mineral: float = 0.7, collagen: float = 0.2, water: float = 0.1, name: str = 'dentine') -> Material`

Dentine (or cementum) as hydroxyapatite + collagen + water by mass.

### `fit_dose_response`

`fit_dose_response(dose: Sequence[float], intensity: Sequence[float], model: str = 'SSE', sigma: Sequence[float] | None = None, weighting: str = 'none', p0: Sequence[float] | None = None, max_dose: float | None = None, De_min: float = 0.0, scale_errors: bool = True) -> DoseResponseResult`

Fit an additive-dose curve and return De with its uncertainty.

**Parameters**

- `dose`, `intensity`: Added laboratory dose (0 for the natural aliquot) and ESR intensity.
- `model`: One of `"SSE"`, `"EXPLIN"`, `"DSE"`, `"LIN"`.
- `sigma`: Per-point 1-sigma uncertainty of the intensity. Overrides `weighting`.
- `weighting`: `"none"` (equal weights) or `"1/I^2"` (relative errors, sigma ∝ I), the latter being common practice in ESR dating.
- `p0`: Optional initial guess, in the order of `PARAM_NAMES[model]`.
- `max_dose`: Discard points with added dose above this value (Dmax test).
- `De_min`: Lower bound of De (default 0). Use `-np.inf` to see where the data really extrapolate: a negative De flags an inconsistent series (e.g. weak low-dose points measured too low), which a bound at 0 would hide.
- `scale_errors`: With explicit `sigma`, inflate the covariance by `chi2_red` when it exceeds 1 (Birge ratio), so scatter not explained by the per-point errors (aliquot inhomogeneity, positioning in the cavity) reaches the De uncertainty.

### `matrix_dose_rates`

`matrix_dose_rates(U: float = 0.0, Th: float = 0.0, K: float = 0.0, factors: ConversionFactors | None = None, U_ra226: float | None = None) -> dict[str, float]`

Dry infinite-matrix alpha/beta/gamma dose rates (Gy/ka).

`U` and `Th` in ppm, `K` in %. The Th chain is taken in secular
equilibrium; so is the U chain unless `U_ra226` (226Ra and daughters
expressed as ppm of U in equilibrium, as measured by gamma spectrometry
through 214Pb/214Bi) differs from `U` (238U).

### `mixture`

`mixture(name: str, parts: list[tuple[Material, float]]) -> Material`

Mix materials by mass fraction: `[(material, mass_fraction), ...]`
(fractions are normalised).

### `sediment_material`

`sediment_material(quartz: float = 1.0, calcite: float = 0.0, dolomite: float = 0.0, kaolinite: float = 0.0, illite: float = 0.0, feldspar: float = 0.0, iron_oxide: float = 0.0, name: str = 'sediment') -> Material`

Dry sediment from its mineral mass fractions (normalised), e.g.
`sediment_material(quartz=0.6, calcite=0.3, kaolinite=0.1)`; water is
added separately (`ToothLayers.sediment_water`).

### `solve_age`

`solve_age(De: float, components: Sequence[DoseRateComponent], t_max: float = 100000.0) -> AgeResult`

Find the age T (ka) such that the accumulated dose equals `De` (Gy).

### `water_correction`

`water_correction(rate: float, water: float, radiation: str) -> float`

Correct a dry dose rate for water content.

`water` is the ratio mass of water / mass of dry material (e.g. 0.15).
