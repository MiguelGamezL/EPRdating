# One-group beta attenuation

`eprdating.onegroup`

One-group (double-P0) beta transport in planar layers.

- `Theory`: O'Brien, Samson, Sanna & McLaughlin (1964) Nucl. Sci. Eng. 18, 90;
coefficients as in Prestwich & Chan (2000) Radiat. Phys. Chem. 59, 221.
Applied to tooth enamel by Brennan et al. (1997) Radiat. Meas. 27, 307 (ROSY).

For each beta emitter, the spectrum is replaced by a single group of
electrons at the mean energy E. The forward/backward fluences obey

    ± dΦ±/dz + (2 μa + μs) Φ± = Y(z) + μs Φ∓

with z the mass depth (g/cm²), Y the decays per gram, μa = S(E)/E (collision
stopping power over energy) and μs the Lewis transport cross section of a
screened Coulomb potential. The dose is D = μa E (Φ+ + Φ−); in an infinite
homogeneous medium D = E Y. Every layer is homogeneous, so the solution in a
layer is a constant plus a rising and a decaying exponential, with
attenuation ν = 2 sqrt(μa (μa + μs)). Fluences are continuous at interfaces.

**Contents:** [`Emitter`](#emitter), [`Layer`](#layer), [`Material`](#material), [`ToothLayers`](#toothlayers), [`bethe_stopping`](#bethe_stopping), [`lewis_mus`](#lewis_mus), [`mean_dose`](#mean_dose), [`mixture`](#mixture), [`solve_fluence`](#solve_fluence), [`weighted_fraction`](#weighted_fraction)

### `Emitter`

*dataclass* `Emitter(name: str, energy: float, weight: float)`

Emitter(name: 'str', energy: 'float', weight: 'float')

### `Layer`

*dataclass* `Layer(material: Material, thickness: float, source: float = 0.0)`

Layer(material: 'Material', thickness: 'float', source: 'float' = 0.0)

### `Material`

*dataclass* `Material(name: str, fractions: dict)`

Mass fractions of elements (normalised on construction).

**Members**

- `coefficients(self, E: float, scatter_factor: float = 1.0) -> tuple[float, float]` — (μa, μs) in cm²/g at energy E (MeV).
- `with_water(self, water: float) -> Material` — Add `water` grams of water per gram of dry material.

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
`source` ('enamel', 'dentine', 'cementum', 'sediment') for energy E.
- `nominal_values(self) -> dict`

### `bethe_stopping`

`bethe_stopping(z_over_a, I_eV, E)`

Collision stopping power for electrons, MeV cm²/g (no density effect).

### `lewis_mus`

`lewis_mus(Z, A, E)`

Lewis transport (scattering) cross section, cm²/g (Prestwich & Chan eq. 12-13).

### `mean_dose`

`mean_dose(sol, i: int, z0: float, z1: float, E: float) -> float`

Mean dose (E × decays/g units) over local depths [z0, z1] of layer i.

### `mixture`

`mixture(name: str, parts: list[tuple[Material, float]]) -> Material`

Mix materials by mass fraction: `[(material, mass_fraction), ...]`.

### `solve_fluence`

`solve_fluence(layers: list[Layer], E: float, scatter_factor: float = 1.0)`

Solve the double-P0 equations for one emitter of energy E.

`layers` runs from left to right; the first and last must be half-spaces
(`thickness = inf`). Returns a list of per-layer solutions
`(mua, mus, nu, P, A, B, t)` with, in local depth ζ ∈ [0, t]:

```
Φ+ = P + A e^{ν(ζ - t)} + B e^{-ν ζ}
Φ- = P + A g_A e^{ν(ζ - t)} + B g_B e^{-ν ζ}
```

### `weighted_fraction`

`weighted_fraction(emitters, frac_fn) -> float`

Energy-weighted mean of `frac_fn(E)` over a list of emitters.
