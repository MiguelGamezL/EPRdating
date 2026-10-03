# Dose-response and De

`eprdating.doseresponse`

Dose-response curve (DRC) fitting and equivalent dose (De) estimation.

Additive-dose protocol: aliquots receive laboratory doses `D` on top of the
natural (archaeological) dose. The ESR intensity is modelled as a function of
the *total* dose `D + De` and the curve is extrapolated back to zero
intensity, where `D = -De`.

Models (`x = D + De`):


```text
"SSE"     single saturating exponential  Imax * (1 - exp(-x/D0))
"EXPLIN"  exponential + linear  Imax * (1 - exp(-x/D0)) + m*x
"DSE"     double saturating exponential
"LIN"     linear m * x (only for doses far below saturation)
```

**Contents:** [`DoseResponseResult`](#doseresponseresult), [`bootstrap_De`](#bootstrap_de), [`fit_dose_response`](#fit_dose_response)

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

### `bootstrap_De`

`bootstrap_De(result: DoseResponseResult, n: int = 1000, seed: int | None = None) -> np.ndarray`

Residual bootstrap of De. Returns the array of resampled De values.

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
