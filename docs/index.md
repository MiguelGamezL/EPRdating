# EPRdating

![EPRdating](img/logo/eprdating_logo.svg#only-light){ width="560" }
![EPRdating](img/logo/eprdating_logo_dark.svg#only-dark){ width="560" }

Open ESR (EPR) dating of tooth enamel in Python, from raw measurements to an
age with its uncertainty:

```text
EPR spectra ──► intensities ──► dose-response ──► De ─┐
                                                      ├──► age (Monte Carlo)
HPGe spectra ─► U, Th, K ─► dose rates ───────────────┘
```

Each step is a function you can call on its own, and every result comes with
its uncertainty.

| Step | What EPRdating does |
|---|---|
| **EPR spectra** | Reads Bruker (Xepr BES3T, ESP/WinEPR), ASCII and column files; corrects baselines, normalises by power, gain and mass; aligns spectra; fits a line-shape template (empirical, or simulated with EPRAYA and broadened by modulation and time constant) with errors by noise injection. |
| **Dose-response** | Linear, single and double saturating exponential, exponential + linear; weighting, Dmax test, Birge scaling, bootstrap. |
| **Sediment** | HPGe gamma spectra (ORTEC, Canberra, N42, ASCII…) analysed by the comparative method against IAEA reference materials, with automatic energy calibration and a 226Ra/238U equilibrium check. |
| **Dose rate** | Conversion factors (Adamiec & Aitken 1998, Guérin et al. 2011, Liritzis et al. 2013), water, cosmic (Prescott & Hutton 1994), U uptake (EU, LU, US), U-series ingrowth with radon loss, one-group beta attenuation in enamel/dentine/cementum layers, energy-dependent alpha efficiency, sediment disequilibrium. |
| **Age** | Early, linear and US uptake, combined U-series/ESR (US-ESR) and CSUS-ESR; Monte Carlo over every input, geometry included. |
| **Figures** | Spectra, dose-response, dose-rate budget, age distribution, gamma spectra and per-line contents. |

The physics is validated against ROSY 2.0 (ages within −0.6 to +0.7 %),
DATA (−2.6 to +1.6 % with DATA's beta convention), 58 published ages from
nine studies (all within 10 %, 56 within the published 1σ),
and third-party EPR and gamma files — see
[Validation](validation.md).

## Install

```bash
pip install eprdating                    # core: numpy + scipy
pip install "eprdating[plot]"            # + matplotlib, for the figures
pip install "eprdating[spectra]"         # + EPRAYA, simulated line shapes
pip install "eprdating[gamma-formats]"   # + becquerel: Canberra .cnf, ORTEC .spc, IEC 61455
```

Until the package is on PyPI, install it from GitHub:
`pip install "eprdating[plot] @ git+https://github.com/MiguelGamezL/EPRdating"`

Python 3.10–3.13. Works in scripts, Jupyter and Google Colab (figures show
inline); [`examples/quickstart.ipynb`](https://github.com/MiguelGamezL/EPRdating/blob/main/examples/quickstart.ipynb)
runs the whole chain on simulated data. In Colab, install with

```python
%pip install "eprdating[plot] @ git+https://github.com/MiguelGamezL/EPRdating"
```

and read your files after uploading them or mounting Google Drive
(`from google.colab import drive; drive.mount("/content/drive")`).

## Where to go next

- **[User guide](guide.md)** — the three steps with working code.
- **[Figures](plotting.md)** — what each plotting function draws.
- **[API reference](api/index.md)** — every public function and class.
- **[Validation](validation.md)** — what has been checked, against what.

## Units

Field in mT, frequency in GHz, dose in Gy, dose rate in Gy/ka, age in ka,
U and Th in µg/g (ppm), K in %. Inputs accept a number, a `(value, sigma)`
tuple or a `Value`.

## Authors

Miguel Enrique Gámez López ([ORCID](https://orcid.org/0000-0001-7831-3291)), Carol Jiseth Ospina Umaña,
Juan Sebastián Castro Millán and Ovidio Amado Almanza Montero —
Grupo de Física Aplicada, Departamento de Física, Universidad Nacional de
Colombia, Sede Bogotá. To cite the software: version 0.1.0,
[doi:10.5281/zenodo.23192195](https://doi.org/10.5281/zenodo.23192195); all
versions, [doi:10.5281/zenodo.23192194](https://doi.org/10.5281/zenodo.23192194).
