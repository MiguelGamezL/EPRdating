# Environmental histories

`eprdating.history`

Piecewise-constant histories of environmental parameters.

The burial environment of a tooth need not stay as it is today: a wetter
period raises the sediment water and lowers the external dose rate, and
erosion or deposition changes the overburden and so the cosmic dose rate.
A `History` gives one parameter as a sequence of constant values,
from the present backwards in time:

```
# 10 % water today and back to 12 ka, 25 % from 12 to 60 ka, 15 % before
water = History([(0.10, 0.03), (0.25, 0.05), (0.15, 0.05)], breaks=[12, 60])
```

Times are in ka before present and the last value holds back to any age.
In the Monte Carlo every segment value is sampled independently; the break
times are fixed.

Histories are accepted by `ToothSample` for the
sediment water (`Sediment(water=History(...))`), the cosmic dose rate and
the external gamma dose rate; `cosmic_history`
turns a burial-depth history into a cosmic dose-rate history. The internal
components (enamel, dentine and cementum) keep the present-day values.

**Contents:** [`History`](#history), [`integrate`](#integrate)

### `History`

*dataclass* `History(values: Sequence[ValueLike], breaks: Sequence[float] = ())`

Piecewise-constant parameter history, in ka before present.

`values[0]` holds from today back to `breaks[0]`, `values[i]`
between `breaks[i-1]` and `breaks[i]`, and the last value before
the last break.

**Members**

- `at(self, t: float) -> float` — Nominal value at `t` ka before present.
- `mean(self, T: float) -> float` — Time-weighted nominal mean over the last `T` ka.
- `nominal(self) -> list[float]`
- `sample(self, rng: np.random.Generator, n: int) -> list[np.ndarray]` — `n` draws of every segment value, sampled independently.

### `integrate`

`integrate(breaks: Sequence[float], rates: Sequence[float], T: float) -> float`

∫_0^T of a piecewise-constant function (`rates` between `breaks`).
