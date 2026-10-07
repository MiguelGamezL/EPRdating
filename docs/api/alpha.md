# Alpha efficiency

`eprdating.alpha`

Energy-dependent alpha efficiency.

The luminescence/ESR signal induced by an alpha particle is, to a good
approximation, proportional to the length of its track rather than to the
energy it deposits (Zimmerman 1971; Aitken & Bowman 1975; see Adamiec & Aitken
1998, p. 38). An alpha of initial energy E deposits E over its range R(E), so
its efficiency per unit dose scales as R(E)/E. Relative to a reference energy
(ROSY uses 5.3 MeV):

```
k(E) = k_ref · [R(E)/E] / [R(E_ref)/E_ref]
```

R(E) is the CSDA range in the medium from the Bethe stopping power for alpha
particles (below 1 MeV the range is extrapolated with R ∝ E^1.5). Alpha
energies per decay are those of Adamiec & Aitken (1998), Tables 2 and 3.

**Contents:** [`alpha_range`](#alpha_range), [`escape_fractions`](#escape_fractions), [`k_ratio`](#k_ratio), [`natural_u_k_ratio`](#natural_u_k_ratio), [`segment_k_ratios`](#segment_k_ratios), [`slab_fractions`](#slab_fractions), [`stopping_power`](#stopping_power), [`surface_fraction`](#surface_fraction), [`th232_k_ratio`](#th232_k_ratio)

### `alpha_range`

`alpha_range(E: float, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> float`

CSDA range of an alpha particle of energy E (MeV), in g/cm².

### `escape_fractions`

`escape_fractions(thickness_um: float, strip_outer_um: float = 0.0, strip_inner_um: float = 0.0, density: float = 3.0, energy: bool = False, e_ref: float = 5.3, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397}), outer_material: Material | None = None, inner_material: Material | None = None) -> dict`

Alpha escape for a layer, per U-series segment and for the 232Th chain.

Returns `{"own": {segment: f}, "outer": {...}, "inner": {...}}` with the
U-series segments of `ALPHA_EMITTERS` and `"Th232"`; each value is
averaged over the emitters of the segment weighted by their alpha dose
(times the relative efficiency k(E)/k_ref when `energy`). Ranges are
those in `material` at `density` (g/cm³); the incoming fractions are
scaled by the mass-range ratio of `outer_material` / `inner_material`
(the neighbouring media) to `material` when they are given.

### `k_ratio`

`k_ratio(E: float, e_ref: float = 5.3, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> float`

Alpha efficiency at energy E relative to that at `e_ref`.

### `natural_u_k_ratio`

`natural_u_k_ratio(e_ref: float = 5.3, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> float`

k/k_ref for natural U in secular equilibrium.

### `segment_k_ratios`

`segment_k_ratios(e_ref: float = 5.3, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> dict[str, float]`

Dose-weighted alpha-efficiency ratio k/k_ref for each U-series segment.

### `slab_fractions`

`slab_fractions(thickness_um: float, strip_outer_um: float, strip_inner_um: float, range_um: float) -> tuple[float, float, float]`

`(own, from_outer, from_inner)` for the measured part of a layer.

- `own`: fraction of the layer's own infinite-matrix alpha dose kept;
- `from_outer`: fraction of the outer medium's infinite-matrix alpha dose received (likewise `from_inner`).
The measured part is what remains after stripping `strip_outer_um` and
`strip_inner_um` from the two faces.

### `stopping_power`

`stopping_power(E: float, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> float`

Bethe electronic stopping power for alpha particles, MeV cm²/g.

### `surface_fraction`

`surface_fraction(a_um: float, b_um: float, range_um: float) -> float`

Mean of e(x/R) over distances x from a surface between `a_um` and
`b_um`: the fraction of the alpha track length lost through that
surface (or received from beyond it), averaged over the slice.

### `th232_k_ratio`

`th232_k_ratio(e_ref: float = 5.3, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> float`

k/k_ref for the 232Th chain in secular equilibrium.
