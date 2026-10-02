# Validation against ROSY 2.0 — findings

Reference outputs were produced with the original ROSY 2.0 program (Brennan et
al. 1999), run under Wine and driven by `tools/rosy_harness`. ROSY itself is
not part of this repository; only its inputs (`cases.json`), parsed outputs
(`results.json`) and raw output files (`raw/`) are stored here.

All comparisons are black-box: ROSY output vs. EPRdating's model evaluated at
the age ROSY found. "Time-averaged dose rate" means De/T, which is what ROSY
prints in its breakdown.

## Agreement (locked in `test_rosy_reference.py`)

| Piece | Result |
|---|---|
| Conversion factors | ROSY uses **Adamiec & Aitken (1998)**: external gamma agrees to < 0.1 %. |
| Water correction | Same Zimmerman coefficients (gamma 1.14): < 0.1 %. |
| Cosmic dose rate | Prescott & Hutton at **55°N and sea level** (no latitude or altitude correction): within 1 %. |
| U-series ingrowth (EU, LU, initial 234U/238U 1.0–3.0) | Enamel alpha implies k_eff = 0.144–0.154 for 50 ka–2 Ma with input k = 0.15; the drift comes from ROSY's energy-dependent alpha efficiency, not from ingrowth. |
| Radon | ROSY's "Fraction of Radon" is the fraction **retained** (1.0 = no loss). |
| End-to-end ages (EU and LU, 35–330 ka) | EPRdating within **+1.5 to +2.4 %** of ROSY using beta factors measured from ROSY. |

## Beta geometry factors measured from ROSY

Fraction of the infinite-matrix beta dose rate received by the enamel
(defaults: enamel density 3.0, dentine 2.82, dentine 2000 µm, no stripping,
EU, ages 0.6–2 Ma unless noted). These are the targets for the one-group
implementation.

| Source → enamel | 300 µm | 600 µm | 1000 µm | 1500 µm | 3000 µm |
|---|---|---|---|---|---|
| Enamel U (self-dose) | 0.391 | 0.563 | 0.692 | 0.781 | 0.888 |
| Dentine U | 0.280 | 0.200 | 0.140 | 0.099 | 0.051 |
| Sediment U | 0.303 | — | 0.152 | — | 0.056 |

Sediment, 1000 µm enamel: U 0.152, Th 0.143, K 0.148. Stripping 50 / 100 µm
from the outside lowers the sediment-U factor to 0.136 / 0.123.

The effective factors **depend on the chain segment**: enamel self-dose goes
from 0.645 (51 ka) to 0.696 (2 Ma) as 230Th daughters grow in, and drops by
~10 % with total radon loss. A single factor per source is therefore an
approximation; the one-group solver should attenuate each nuclide group with
its own coefficient.

## ROSY behaviour to be aware of

- **"Present Ratio" option looks unreliable.** With De = 300 Gy and enamel U
  only, present ratios 1.5 and 2.0 and initial ratio 2.0 all give exactly
  464.56 ka; present 1.4 gives the same age as initial 1.4 (226.15 ka at
  De = 100 Gy); present 1.2 does not converge. Present 3.0 gives an *older*
  age than initial 3.0, the opposite of what a back-correction implies.
  Reference cases therefore use initial ratios only.
- **Non-convergence** above roughly 4 Ma (and for present ratio 1.2): ROSY
  shows a dialog and writes no output. Those cases are listed under `failed`
  in `cases.json`.
- Imported files reset the alpha-efficiency option to "Varies with energy";
  "Constant" can only be set by hand, so all reference cases use the
  energy-dependent option.
- Ages ROSY reports for "combination uptake" (CU) use the per-tissue uptake
  models; EU and LU apply early or linear uptake to every tissue.

## Brennan et al. (1997) samples rerun in ROSY 2.0

`rosy_reference/brennan1997_tables.json` transcribes Tables 1 and 2 of
Brennan et al. (1997, *Radiation Measurements* 27, 307–314): six real teeth
with inputs and the outputs of the 1997 ROSY. The `BR97_*` cases rerun them in
ROSY 2.0 with the missing inputs left at defaults (dentine 2000 µm, 234U/238U
= 1.0, k = 0.15).

| Sample | EU age 1997 → 2.0 (ka) | LU | CU |
|---|---|---|---|
| 1 | 112 → 117.7 | 145 → 149.8 | 139 → 144.7 |
| 2 | 32 → 32.5 | 44 → 43.9 | 37 → 38.5 |
| 3 | 19 → 20.1 | 30 → 30.8 | 25 → 27.1 |
| 4 | 82 → 94.8 | 144 → 165.7 | 84 → 96.9 |
| 5 | 35 → 38.2 | 59 → 62.9 | 45 → 49.5 |
| 6 | 39 → 39.0 | 39 → 39.1 | 39 → 39.1 |

Ages agree within ~0–9 % except sample 4 (high enamel U, 80 m burial,
+16 %). The component breakdowns differ more, which points to differences
between the 1997 and 2.0 versions rather than to the inputs:

- Gamma is ~3–5 % higher in 1997, consistent with the older (Nambi & Aitken
  1986) conversion factors; ROSY 2.0 uses Adamiec & Aitken (1998).
- Cosmic dose rates differ (e.g. 89 vs 115 µGy/a at 4.3 m): 1997 used
  Prescott & Stephan (1982), 2.0 behaves like Prescott & Hutton (1994).
- Dentine beta in 1997 is ~2–3 times the ROSY 2.0 value. ROSY 2.0 prints
  time-averaged dose rates (De/T); the 1997 table may list present-day rates,
  which for linear uptake are about twice the average. Unconfirmed.

The paper describes the one-group method only qualitatively ("double-P0"
approximation after O'Brien et al. 1964; absorption cross section = stopping
power / energy at the mean beta energy; Lewis transport cross section for
scattering; flux continuity at interfaces; dose in each layer = sum of a rising
and a decaying exponential). The equations are in the follow-up report it
announces for *Ancient TL* and in O'Brien et al. (1964).

## One-group beta attenuation (`eprdating.onegroup`)

Implemented from O'Brien et al. (1964), with the absorption coefficient
μa = S(E)/E (Bethe collision stopping power) and the Lewis transport cross
section μs of Prestwich & Chan (2000, eqs. 11–13). The coefficients reproduce
the values Prestwich & Chan quote for H and O at 0.7 MeV (μs exactly, μa within
0.5 %). Each beta emitter of the U, Th and K chains (energies from Adamiec &
Aitken 1998) is transported separately at its mean energy; layers are
sediment | (cementum) | enamel | dentine | sediment, with hydroxyapatite enamel,
silica sediment and an indicative dentine (70 % mineral, 20 % collagen, 10 %
water). No parameter was fitted to ROSY: the theoretical scattering cross
section is used (Prestwich & Chan's factor of 0.42, which matches Monte Carlo,
reproduces ROSY much worse).

Against the factors measured from ROSY 2.0:

| Source → 1000 µm enamel | EPRdating | ROSY 2.0 |
|---|---|---|
| Enamel U (self) | 0.702 | 0.692 |
| Dentine U | 0.140 | 0.140 |
| Sediment U / Th / K | 0.147 / 0.136 / 0.145 | 0.152 / 0.143 / 0.148 |

Dentine factors agree to three decimals at all thicknesses (300–3000 µm),
self-dose within 1–2.5 % and sediment within 2–5 %, the remaining differences
probably coming from the enamel and sediment compositions ROSY assumes.

Full ages with this geometry (`test_ages_with_onegroup_within_2_percent_of_rosy`)
agree with ROSY 2.0 within **−1.2 % to +0.9 %** for EU, LU and CU, including
the six real teeth of Brennan et al. (1997) with their thicknesses, stripping
and water contents, using a constant alpha efficiency (ROSY's is energy
dependent).

## Monte Carlo uncertainties

`ToothSample.age_mc` samples the layer thicknesses, stripping, densities and
water contents and recomputes the one-group factors for every draw (~4–5 ms
per draw). For the six teeth of Brennan et al. (1997), with the input errors
of their Table 1, the Monte Carlo spread matches the age errors ROSY 2.0
reports (EU, 1000 draws):

| Sample | ROSY 2.0 (ka) | EPRdating MC (ka) | MC, geometry fixed (ka) |
|---|---|---|---|
| 1 | 117.7 ± 12.2 | 116.8 ± 12.4 | ± 7.9 |
| 2 | 32.5 ± 3.1 | 32.0 ± 3.0 | ± 2.9 |
| 3 | 20.1 ± 2.9 | 20.1 ± 2.8 | ± 2.6 |
| 4 | 94.8 ± 7.6 | 93.8 ± 7.6 | ± 7.6 |
| 5 | 38.2 ± 1.4 | 37.6 ± 1.3 | ± 0.8 |
| 6 | 39.0 ± 1.3 | 39.1 ± 1.3 | ± 1.2 |

Without sampling the geometry the error is underestimated by up to ~40 %
(sample 1, enamel 1828 ± 360 µm). `test_mc_uncertainty_matches_rosy_error`
checks the agreement within 25 % with 400 draws.
