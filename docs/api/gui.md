# Interface

`eprdating.gui`

Interactive interface for Jupyter and Google Colab.

    from eprdating.gui import app
    app()

Three tabs: the equivalent dose from EPR spectra, the sediment U, Th and K
from HPGe spectra, and the age. A fitted De and an analysed sediment flow
into the Age tab. Needs `ipywidgets` and `matplotlib`
(`pip install "eprdating[gui]"`; both come with Google Colab).

Each tab is also available alone (`DePanel`, `GammaPanel`,
`AgePanel`) and can be driven from code, which is how it is tested.

**Contents:** [`AgePanel`](#agepanel), [`App`](#app), [`DePanel`](#depanel), [`GammaPanel`](#gammapanel), [`app`](#app)

### `AgePanel`

*class* `AgePanel()`

Dose rate and age of a tooth (one tab of `app`).
`tooth` builds the `ToothSample` from the form,
`compute` gives the nominal and Monte Carlo ages.

**Members**

- `compute(self)`
- `set_De(self, value: float, sigma: float)`
- `set_sediment(self, sed: Sediment)`
- `settings(self) -> dict`
- `tooth(self) -> ToothSample`

### `App`

*class* `App()`

The three panels in tabs. `app.de`, `app.gamma` and `app.age`
are the panels; `app.widget` is what is displayed.

### `DePanel`

*class* `DePanel(on_de: Callable[[float, float], None] | None = None)`

Interactive equivalent-dose analysis (one tab of `app`).

Every step is also a method, so the panel can be driven from code:
`load`, `compute`, `fit`, `settings`.

**Members**

- `compute(self) -> list[dict]` — Intensity of every aliquot (repeats averaged), then a fit.
- `fit(self)` — Dose-response fit of the ticked points.
- `group_by_dose(self)`
- `load(self, paths) -> list[str]` — Read the spectra among `paths` (companion files are skipped) and
add a row for each. Returns the names read.
- `results(self) -> list[dict]`
- `results_csv(self) -> str`
- `set_file(self, name: str, dose: float | None = None, aliquot: str | None = None, mass_mg: float | None = None, use: bool | None = None)` — Fill a file's row from code.
- `set_point(self, aliquot: str, in_fit: bool)`
- `settings(self) -> dict` — Everything needed to repeat the analysis.
- `window(self) -> IntensityWindow`

### `GammaPanel`

*class* `GammaPanel(on_sediment: Callable[[Sediment], None] | None = None)`

Sediment U, Th and K by the comparative method (one tab of
`app`). Drive it from code with `load`,
`set_roles` and `analyse`.

**Members**

- `analyse(self)`
- `load(self, paths) -> list[str]` — Add files and guess their roles from their names.
- `set_roles(self, sample: str, U: str, Th: str, K: str, background: str | None = None)`
- `settings(self) -> dict`

### `app`

`app() -> App`

Open the interface (in a notebook cell: `app()`).
