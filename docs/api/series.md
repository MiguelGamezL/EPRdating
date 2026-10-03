# U-series ingrowth

`eprdating.series`

U-series disequilibrium after uranium uptake.

Uranium taken up by a tooth arrives without its daughters: `230Th` (and
everything below it) grows in over ~10^5 years, and `234U` may be in excess
of `238U`. The dose rate of freshly incorporated U is therefore much lower
than the secular-equilibrium value given by the conversion factors.

We split the U-238 chain into segments whose activity relative to `238U`
evolves as (`tau` = time since that U was incorporated, `r0` = initial
234U/238U activity ratio, 230Th initially absent):

* `U238`  : 238U → 234Th → 234Pa                 ratio 1
* `U234`  : 234U                                  1 + (r0-1) e^{-λ4 τ}
* `Th230` : 230Th → 226Ra                         Bateman ingrowth
* `Rn222` : 222Rn and its daughters to 210Po      as Th230, times (1 - f_Rn)
* `U235`  : 235U → 231Th                          ratio 1
* `Pa231` : 231Pa and its daughters               1 - e^{-λPa τ}

226Ra is assumed in equilibrium with 230Th (half-life 1.6 ka) and `f_Rn`
is the fraction of 222Rn that escapes the tissue (radon loss).

The 234U/238U ratio can be given as the **present-day** (measured) value
(default) or as the initial value. For a present-day ratio, every U parcel is
assumed to show today the measured value, so a parcel incorporated a time
`tau` ago started with `r0 = 1 + (r_now - 1) e^{λ4 τ}`. This is exact for
early uptake and a consistent approximation for continuous uptake. For old
samples with a ratio far from 1 the implied `r0` grows quickly; check it
with `initial_ratio`.

The fraction of the equilibrium dose rate of natural U carried by each
segment, per radiation type, is bundled in `data/u_series_partition.json`.
It is derived from the energies per disintegration of Adamiec & Aitken
(1998, Tables 2, 3 and 5) by `tools/derive_u_series_partition.py`. A
different table can be passed to `USeries` (fractions must sum to 1).

For freshly incorporated U only ~20 % (alpha), ~39 % (beta) and ~2 %
(gamma) of the equilibrium dose rate is delivered.

**Contents:** [`USeries`](#useries), [`activity_ratio_Th230_U238`](#activity_ratio_th230_u238), [`load_partition`](#load_partition)

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

### `activity_ratio_Th230_U238`

`activity_ratio_Th230_U238(tau: float, r0: float = 1.0) -> float`

(230Th/238U) activity ratio after a time `tau` (ka) of a closed system
that started with 234U/238U = `r0` and no 230Th.

### `load_partition`

`load_partition(path: str | None = None) -> dict[str, dict[str, float]]`

Load segment fractions `{radiation: {segment: fraction}}`.
