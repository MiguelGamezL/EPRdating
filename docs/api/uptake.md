# Uranium uptake

`eprdating.uptake`

Uranium-uptake models for dental tissues.

The US model (Grün et al. 1988) describes the U content of a tissue at time
`t` after burial (`0 <= t <= T`, with `T` the age) as

    U(t) = U_m * (t / T) ** (p + 1),       p >= -1

where `U_m` is the U concentration measured today. Special cases:

* `p = -1`  early uptake (EU): all U present since burial.
* `p =  0`  linear uptake (LU).
* `p >  0`  increasingly recent uptake.

Uptake only affects the *U-derived* dose-rate contributions of the tissue
(alpha and beta of enamel, beta of dentine); sediment gamma and cosmic
contributions are taken as constant.

**Contents:** [`EarlyUptake`](#earlyuptake), [`LinearUptake`](#linearuptake), [`USModel`](#usmodel)

### `EarlyUptake`

`EarlyUptake() -> USModel`


### `LinearUptake`

`LinearUptake() -> USModel`


### `USModel`

*dataclass* `USModel(p: float = 0.0)`

Uptake following `(t/T)**(p+1)`; `p=-1` is EU, `p=0` is LU.

**Members**

- `accumulated(self, T: float, G = None) -> float` — Dose accumulated over `T` per unit present-day dose rate.
- `fraction(self, t, T: float)` — U(t)/U_m at time `t` after burial for a sample of age `T`.
- `is_early` *(property)*
