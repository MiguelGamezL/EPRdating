# EPR spectra: template and component fits

`eprdating.spectra.deconvolution`

Spectral deconvolution into known components.

The measured derivative spectrum is modelled as a non-negative combination
of component line shapes (e.g. the axial and orthorhombic CO2- species, the
native organic signal…) plus a free polynomial baseline:

    y(B) ≈ Σ_i a_i · s_i(B) + Σ_k c_k · B^k ,   a_i >= 0

The component shapes `s_i` can come from any source: reference spectra,
analytic line shapes, or spin-Hamiltonian simulations with EPRAYA
(`epraya_backend`). Each shape is normalised to unit
peak-to-peak so the amplitudes `a_i` are peak-to-peak intensities.

Fitting the amplitudes is a bounded *linear* problem, so it is fast and has
a unique solution; the (slow) simulation of the shapes is done once.

**Contents:** [`ComponentBasis`](#componentbasis), [`DeconvolutionResult`](#deconvolutionresult), [`component_vs_dose`](#component_vs_dose), [`normalise_pp`](#normalise_pp)

### `ComponentBasis`

*class* `ComponentBasis(B, components: Mapping[str, Shape], baseline_order: int = 1)`

A set of normalised component shapes on a fixed field grid.

**Members**

- `fit(self, spectrum, nonnegative: bool = True, max_shift: float = 0.0, shift_step: float | None = None, noise = None, n_noise: int = 300, seed: int | None = 0) -> DeconvolutionResult` — Fit amplitudes (peak-to-peak units) and the polynomial baseline.

### `DeconvolutionResult`

*dataclass* `DeconvolutionResult(amplitudes: dict[str, float], baseline: np.ndarray, fitted: np.ndarray, residuals: np.ndarray, r2: float, errors: dict[str, float] = <factory>, shift: float = 0.0, noise_corr_length: int = 1)`

DeconvolutionResult(amplitudes: 'dict[str, float]', baseline: 'np.ndarray', fitted: 'np.ndarray', residuals: 'np.ndarray', r2: 'float', errors: 'dict[str, float]' = <factory>, shift: 'float' = 0.0, noise_corr_length: 'int' = 1)

**Members**

- `component(self, name: str) -> float`

### `component_vs_dose`

`component_vs_dose(basis: ComponentBasis, spectra: Sequence[np.ndarray], component: str, **fit_kw) -> tuple[np.ndarray, np.ndarray, list]`

Amplitude of `component` in each spectrum of a dose series.

Returns `(amplitudes, errors, results)`; feed amplitudes and errors
with their doses to `fit_dose_response`.
`fit_kw` go to `fit`.

### `normalise_pp`

`normalise_pp(y: np.ndarray) -> np.ndarray`

