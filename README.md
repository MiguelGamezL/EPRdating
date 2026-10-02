# EPRdating

ESR/EPR dating in Python, from spectra to ages:

```
spectra ─► intensity ─► dose-response curve ─► De ─► dose rate ─► age
          (p-p, T1–B2,      (SSE, EXP+LIN,          (U, Th, K,     (EU / LU / US,
           deconvolution     DSE, LIN)               water, cosmic,  U-series ingrowth,
           with EPRAYA)                              beta geometry)  Monte Carlo)
```

Inspired by ROSY (Brennan et al. 1997, 1999) and DATA (Grün 2009), written as
an independent, open implementation from the published equations.

> **Status: 0.1.0.dev0 — pre-release.** Do not use it for published ages until
> the validation cases in `tests/validation/` pass. See *Roadmap* and
> *Known limitations*.

## Install

```bash
pip install -e .                 # core: numpy + scipy only
pip install -e '.[spectra]'      # + EPRAYA backend for simulated component shapes
pip install -e '.[dev]'          # + pytest, ruff
```

EPRAYA (Física Aplicada group, Universidad Nacional de Colombia) is optional
because it pulls JAX, numba, ipywidgets and tkinter.

## Quick start

```python
from eprdating import (fit_dose_response, ToothSample, Sediment, BetaGeometry,
                       LinearUptake, cosmic_dose_rate)

drc = fit_dose_response(dose, intensity, model="SSE", weighting="1/I^2")

sample = ToothSample(
    De=(drc.De, drc.De_sigma),
    enamel_U=(0.8, 0.08), dentine_U=(15.0, 1.5),
    sediment=Sediment(U=(2.1, 0.1), Th=(7.5, 0.4), K=(1.1, 0.05), water=(0.15, 0.05)),
    beta=BetaGeometry(internal=(0.6, 0.03), dentine=(0.3, 0.02), external=(0.25, 0.02)),
    cosmic=cosmic_dose_rate(depth_m=1.5, density=1.9, lat_deg=4.6, lon_deg=-74.1, altitude_m=2600),
    uptake_enamel=LinearUptake(), uptake_dentine=LinearUptake(),
)
print(sample.age().summary())
print(sample.age_mc(n=2000, seed=42).summary())
```

Every input accepts a number, a `(value, sigma)` tuple or a `Value`.
Units: Gy, Gy/ka, ka, ppm (U, Th), % (K), mT, GHz.

See `examples/quickstart.py` and `examples/deconvolution_epraya.py`.

## Modules

| module | content |
|---|---|
| `doseresponse` | additive-dose fits (SSE, EXPLIN, DSE, LIN), weighting, Dmax cut, covariance and bootstrap De |
| `dose_rate` | conversion factors (Adamiec & Aitken 1998, Guérin et al. 2011, Liritzis et al. 2013), water correction, alpha efficiency, cosmic dose rate (Prescott & Hutton 1994) |
| `uptake` | US model `U(t) = U_m (t/T)^(p+1)`; EU (p = −1), LU (p = 0) |
| `series` | 230Th / 234U ingrowth after uptake (Bateman, analytic time integrals) |
| `beta` | beta geometry factors; one-group solver planned |
| `age` | generic solver `∫₀ᵀ Ḋ(t) dt = De` and the `ToothSample` model with Monte Carlo |
| `spectra` | peak-to-peak, T1–B2, double integral, non-negative deconvolution, EPRAYA backend |

## Known limitations (v0.1)

1. **U-series ingrowth without radon loss.** Daughter ingrowth after uptake
   (234U, 230Th, 231Pa) is modelled with segment fractions derived from
   Adamiec & Aitken (1998) (`tools/derive_u_series_partition.py`); radon
   escape and a measured initial 234U/238U per tissue are still to come
   (`USeries(r0=...)` accepts a single ratio).
2. **No one-group beta attenuation yet.** Beta geometry factors are user
   inputs (published tables or Monte Carlo, e.g. DosiVox).
3. **No time-varying** water content or burial depth.
4. **Data provenance.** Conversion factors and the Prescott & Stefan F/J/H
   table were transcribed from the DRAC lookup tables. The Adamiec & Aitken
   (1998) set is checked against the paper; the others still need checking.

## Roadmap

- **v0.1** core: De, dose rate, EU/LU/US, age, Monte Carlo; validation against published ROSY/DATA ages.
- **v0.2** spectra: EPRAYA-based deconvolution → De; comparison of intensity methods.
- **v0.3** one-group beta attenuation in planar layers (Brennan et al. 1997).
- **later** US-ESR coupling with U-series data, radon loss, file readers (Bruker, JEOL), JOSS paper.

## Testing

```bash
pytest                      # unit tests
pytest -m validation -rs    # published cases; skips list what data is missing
pytest -m epraya            # EPRAYA integration (needs the extra)
```

## References

- Brennan B.J. et al. (1997) Beta doses in tooth enamel by "one-group" theory and the ROSY ESR dating software. *Radiation Measurements* 27, 307–314.
- Grün R. (2009) The DATA program for the calculation of ESR age estimates on tooth enamel. *Quaternary Geochronology* 4, 231–232.
- Grün R., Schwarcz H.P., Chadam J. (1988) ESR dating of tooth enamel: coupled correction for U-uptake and U-series disequilibrium. *Nuclear Tracks and Radiation Measurements* 14, 237–241.
- Durcan J.A., King G.E., Duller G.A.T. (2015) DRAC: Dose Rate and Age Calculator for trapped charge dating. *Quaternary Geochronology* 28, 54–61.
- Prescott J.R., Hutton J.T. (1994) Cosmic ray contributions to dose rates for luminescence and ESR dating. *Radiation Measurements* 23, 497–500.
- Adamiec G., Aitken M. (1998) Dose-rate conversion factors: update. *Ancient TL* 16, 37–50. doi:10.26034/la.atl.1998.292
- Guérin G., Mercier N., Adamiec G. (2011) Dose-rate conversion factors: update. *Ancient TL* 29, 5–8.

## License

MIT (see `LICENSE`).
