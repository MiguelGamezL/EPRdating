# US-ESR

`eprdating.usesr`

Combined U-series / ESR (US-ESR) dating of teeth.

Grün, Schwarcz & Chadam (1988): the uptake parameter `p` of each dental
tissue is not assumed (EU, LU) but solved together with the age from the
tissue's measured U-series activity ratios.

**Model**

Uranium accumulates as `U(t) = U_m (t/T)^(p+1)` (`USModel`).
Every parcel of uranium arrives with the same 234U/238U ratio `r_in` (that of
the groundwater) and without 230Th. A parcel that has spent a time `τ` in the
tissue shows today

    234U/238U = 1 + (r_in - 1) e^{-λ234 τ}
    230Th/238U = 1 - e^{-λ230 τ} + (r_in - 1) λ230/(λ230 - λ234) (e^{-λ234 τ} - e^{-λ230 τ})

and the tissue ratios are averages over parcels weighted by `dU`. For given
`T` and `p` the measured 234U/238U fixes `r_in` (linear relation), and
the measured 230Th/234U then fixes `p`. The age is the `T` at which the
ESR dose equation is also satisfied, with each tissue's `p(T)` and
`r_in(T)`.

The U-series data alone set a lower bound on the age: the closed-system
(early uptake) U-series age. When the ESR dose is already exceeded at that
bound there is no US-ESR solution (typically uranium leaching); this is
reported, not forced.

**Contents:** [`USESRMC`](#usesrmc), [`USESRResult`](#usesrresult), [`USESRSample`](#usesrsample), [`UseriesData`](#useriesdata), [`closed_system_age`](#closed_system_age), [`incoming_ratio`](#incoming_ratio), [`predicted_ratios`](#predicted_ratios), [`solve_p`](#solve_p), [`th230_u234`](#th230_u234)

### `USESRMC`

*dataclass* `USESRMC(ages: np.ndarray, p_enamel: np.ndarray, p_dentine: np.ndarray, nominal: USESRResult, p_cementum: np.ndarray | None = None, n_failed: int = 0, marginal_below: float = 0.8)`

Monte Carlo result of a US-ESR age.

Draws without a US-ESR solution (ESR dose exceeded before the U-series
bound, or numerical failure) are counted in `n_failed` and excluded from
the statistics. When the solved fraction is below `marginal_below` the
result is flagged `marginal`: the age is conditional on the draws that
admit a solution and should be read with care (often a sample close to
the closed-system U-series bound, or affected by uranium leaching).

**Members**

- `interval68` *(property)*
- `marginal` *(property)*
- `mean` *(property)*
- `n` *(property)*
- `solved_fraction` *(property)* — Fraction of draws with a US-ESR solution.
- `std` *(property)*
- `summary(self) -> str`

### `USESRResult`

*dataclass* `USESRResult(age: float | None, p_enamel: float | None, p_dentine: float | None, r_in_enamel: float | None, r_in_dentine: float | None, min_age: float, status: str, detail: AgeResult | None = None, p_cementum: float | None = None)`

USESRResult(age: 'float | None', p_enamel: 'float | None', p_dentine: 'float | None', r_in_enamel: 'float | None', r_in_dentine: 'float | None', min_age: 'float', status: 'str', detail: 'AgeResult | None' = None, p_cementum: 'float | None' = None)

**Members**

- `summary(self) -> str`

### `USESRSample`

*dataclass* `USESRSample(tooth: ToothSample, enamel: UseriesData | None = None, dentine: UseriesData | None = None, cementum: UseriesData | None = None)`

A tooth for US-ESR dating.

`tooth` carries everything except the uptake history (its uptake models
and 234U/238U fields are ignored for tissues with U-series data).

**Members**

- `age(self) -> USESRResult` — Nominal US-ESR age and uptake parameters.
- `age_mc(self, n: int = 1000, seed: int | None = None, marginal_below: float = 0.8) -> USESRMC` — Monte Carlo over all inputs, U-series ratios included.

### `UseriesData`

*dataclass* `UseriesData(th230_u234: ValueLike, u234_u238: ValueLike)`

Measured U-series activity ratios of one dental tissue.

### `closed_system_age`

`closed_system_age(th230_u234_meas: float, u234_u238: float, t_max: float = 5000.0) -> float`

Closed-system (early uptake) U-series age in ka; `inf` if beyond range.

### `incoming_ratio`

`incoming_ratio(T: float, p: float, u234_u238: float) -> float`

234U/238U of the incoming uranium that reproduces the measured ratio.

### `predicted_ratios`

`predicted_ratios(T: float, p: float, r_in: float) -> tuple[float, float]`

Present-day (234U/238U, 230Th/234U) of a tissue for age T, uptake p and
incoming ratio r_in.

### `solve_p`

`solve_p(T: float, th230_u234_meas: float, u234_u238: float) -> float | None`

Uptake parameter p of a tissue for age T, or None if T is younger than
the tissue's closed-system U-series age (no p >= -1 fits).

### `th230_u234`

`th230_u234(T: float, p: float, u234_u238: float) -> float`

230Th/234U predicted for a tissue with the measured 234U/238U.
