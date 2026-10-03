# Fixed beta factors

`eprdating.beta`

Beta-dose geometry for tooth enamel: fixed factors.

Two ways to describe how much of the beta dose of each source reaches the
dated enamel:

* `BetaGeometry` (this module): fixed fractions of the infinite-matrix
  beta dose rate, e.g. from published tables or Monte Carlo (DosiVox).
* `ToothLayers`: a layered geometry whose
  attenuation is computed with one-group transport theory for every emitter,
  as in ROSY (Brennan et al. 1997). This is the recommended option.

**Contents:** [`BetaGeometry`](#betageometry)

### `BetaGeometry`

*dataclass* `BetaGeometry(internal: ValueLike, dentine: ValueLike, external: ValueLike)`

Fractions of infinite-matrix beta dose rate reaching the dated enamel.

- `internal`: enamel self-dose (U inside the enamel layer).
- `dentine`: from U in the adjacent dentine.
- `external`: from the sediment (or cement) on the outer side.

**Members**

- `values(self) -> dict`
