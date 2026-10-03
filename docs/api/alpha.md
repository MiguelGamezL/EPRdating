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

**Contents:** [`alpha_range`](#alpha_range), [`k_ratio`](#k_ratio), [`natural_u_k_ratio`](#natural_u_k_ratio), [`segment_k_ratios`](#segment_k_ratios), [`stopping_power`](#stopping_power)

### `alpha_range`

`alpha_range(E: float, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> float`

CSDA range of an alpha particle of energy E (MeV), in g/cm².

### `k_ratio`

`k_ratio(E: float, e_ref: float = 5.3, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> float`

Alpha efficiency at energy E relative to that at `e_ref`.

### `natural_u_k_ratio`

`natural_u_k_ratio(e_ref: float = 5.3, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> float`

k/k_ref for natural U in secular equilibrium.

### `segment_k_ratios`

`segment_k_ratios(e_ref: float = 5.3, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> dict[str, float]`

Dose-weighted alpha-efficiency ratio k/k_ref for each U-series segment.

### `stopping_power`

`stopping_power(E: float, material: Material = Material(name='hydroxyapatite', fractions={'Ca': 0.39893929409703627, 'P': 0.1849904540450362, 'O': 0.41406351096042854, 'H': 0.0020067408974989397})) -> float`

Bethe electronic stopping power for alpha particles, MeV cm²/g.
