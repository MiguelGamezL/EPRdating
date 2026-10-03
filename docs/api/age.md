# Ages

`eprdating.age`

Age calculation: solve ∫_0^T Ḋ(t) dt = De for the age T.

**Generic layer**

`DoseRateComponent` is one contribution to the dose rate with its
present-day (equilibrium) value, an optional uptake model and an optional
U-series ingrowth function. `solve_age` combines any list of them.

**Tooth-enamel layer**

`ToothSample` assembles the usual components of ESR dating of
enamel (internal alpha and beta, beta from dentine and sediment, gamma,
cosmic) and provides a Monte Carlo age with `age_mc`.

**Contents:** [`AgeMC`](#agemc), [`AgeResult`](#ageresult), [`DoseRateComponent`](#doseratecomponent), [`ToothSample`](#toothsample), [`solve_age`](#solve_age)

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

### `DoseRateComponent`

*dataclass* `DoseRateComponent(name: str, rate: float, uptake: USModel | None = None, G: Callable[[float], float] | None = None)`

One contribution to the total dose rate (Gy/ka).

- `rate`: present-day dose rate the source would deliver in secular equilibrium with its present U content.
- `uptake`: None for a constant source, else a `USModel`.
- `G`: time-integrated activity ratio of incorporated U (`tau -> ∫ D/D_eq`); None means secular equilibrium.

**Members**

- `accumulated(self, T: float) -> float` — Dose (Gy) delivered by this component over `T` ka.

### `ToothSample`

*dataclass* `ToothSample(De: ValueLike, enamel_U: ValueLike, dentine_U: ValueLike, sediment: Sediment, beta: BetaGeometry | ToothLayers, cosmic: ValueLike, gamma: ValueLike | None = None, k_alpha: ValueLike = 0.13 ± 0.02, dentine_water: ValueLike = 0.0, uptake_enamel: USModel = <factory>, uptake_dentine: USModel = <factory>, u234_u238_enamel: ValueLike = 1.0, u234_u238_dentine: ValueLike = 1.0, radon_loss_enamel: ValueLike = 0.0, radon_loss_dentine: ValueLike = 0.0, enamel_water: ValueLike = 0.0, cementum_U: ValueLike = 0.0, cementum_water: ValueLike = 0.0, uptake_cementum: USModel = <factory>, u234_u238_cementum: ValueLike = 1.0, radon_loss_cementum: ValueLike = 0.0, ingrowth: bool = True, partition: dict | None = None, factors: str = 'guerin_2011', sample_geometry: bool = True, alpha_efficiency: str = 'constant', u234_u238_is: str = 'present', alpha_eref: float = 5.3, beta_by_segment: bool = True)`

ESR dating of tooth enamel with the classical component model.

- `Concentrations`: U in ppm. Dose in Gy, dose rates in Gy/ka, ages in ka.

**Parameters**

- `De`: equivalent dose of the enamel.
- `enamel_U`, `dentine_U`: present-day U in each tissue.
- `sediment`: U, Th, K and water content of the surrounding sediment.
- `beta`: either fixed geometry factors (`BetaGeometry`) or a layered geometry (`ToothLayers`), in which case beta attenuation is computed with one-group theory for each emitter and U-series segment (as ROSY does). Uncertainties on the layer thicknesses, stripping and densities are sampled in `age_mc` together with the water contents.
- `gamma`: external gamma dose rate. If None, computed from `sediment` as an infinite matrix (use in-situ measurements when available).
- `cosmic`: cosmic dose rate (see `cosmic_dose_rate`).
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

### `solve_age`

`solve_age(De: float, components: Sequence[DoseRateComponent], t_max: float = 100000.0) -> AgeResult`

Find the age T (ka) such that the accumulated dose equals `De` (Gy).
