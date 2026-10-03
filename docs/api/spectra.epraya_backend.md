# EPR spectra: EPRAYA simulations

`eprdating.spectra.epraya_backend`

Optional EPRAYA backend: component line shapes from spin-Hamiltonian
simulations of powder samples.

Install with `pip install eprdating[spectra]`. EPRAYA is imported lazily, so
the rest of the library works without it (and without JAX).

Only two EPRAYA entry points are used (`Start` and `Powder`), which keeps
the coupling between both projects small.

**Contents:** [`Species`](#species), [`basis_from_species`](#basis_from_species), [`simulate_species`](#simulate_species)

### `Species`

*dataclass* `Species(g: Sequence[float], Hpp: Sequence[float] = (0.0, 0.3), S: float = 0.5, eta: float = 0.5, extra: dict[str, object] = <factory>)`

Spin-Hamiltonian description of one paramagnetic species (S=1/2 default).

- `g`: isotropic float or [gx, gy, gz] principal values.
- `Hpp`: [Gaussian, Lorentzian] peak-to-peak linewidths in mT, as in EPRAYA.
- `eta`: Gaussian weight of the Voigt profile (EPRAYA convention).
- `extra`: any other `Ham` attribute understood by EPRAYA (A, I, Nucl, ...).

### `basis_from_species`

`basis_from_species(B, species: dict[str, Species], freq_GHz: float, baseline_order: int = 1, points: int | None = None, grid: int = 70) -> ComponentBasis`

Build a `ComponentBasis` on the field grid `B` (mT) by
simulating every species with EPRAYA.

### `simulate_species`

`simulate_species(species: Species, freq_GHz: float, field_range_mT: tuple[float, float], points: int = 2048, grid: int = 70) -> tuple[np.ndarray, np.ndarray]`

Simulate the first-derivative powder spectrum of one species.

Returns `(B_mT, intensity)`. `grid` is EPRAYA's Delaunay `M`.
