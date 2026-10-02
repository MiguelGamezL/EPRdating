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
