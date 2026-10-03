# Validation against DATA (Grün 2009)

Two campaigns: EU/LU dose rates and ages (below), and the U-series/ESR
part (US-ESR and CS-US, [second half](#u-seriesesr-us-esr-and-cs-us)).

DATA is the DOS program most ESR tooth ages of the last 15 years were
calculated with. The original `DATA.EXE` was run unchanged as a black box in
DOSBox-X (`tools/data_harness`) on 41 cases built around a base tooth. Each
case gives EU and LU dose rates and ages, so there are 82 runs in total.

## Base case and what is varied

The base tooth (B0) has these inputs:

- De 100 Gy.
- Enamel U 1 ppm, 234U/238U 1.2, k 0.13.
- Dentine U 10 ppm on both sides, dentine water 5 %.
- Enamel 1000 µm, with 50 µm removed from each side, density 2.95 g/cm³.
- Sediment: U 2 ppm, Th 6 ppm, K 1 %, water 10 %.
- Depth 1 m.

The cases vary one input at a time:

| Group | Cases |
|---|---|
| Enamel thickness | 300, 600, 1500, 2500 µm |
| Removed thickness | 0, 200 µm |
| Neighbour of the enamel | dentine–sediment; sediment on both sides (the geometry cases repeated with sediment on both sides) |
| U | enamel 0.1 and 5 ppm; dentine 1 and 50 ppm; sediment only (no U in the tooth) |
| 234U/238U | 1.0, 1.6, 2.5 |
| Rn loss | 50 %, 100 % |
| De | 5, 30, 300, 1000, 3000 Gy (ages 4 ka – 2.4 Ma) |
| Water | sediment 0 and 25 %; dentine 0 and 20 % |
| Depth | 0, 5, 20 m |
| Other | k = 0.10; enamel density 2.8; "beta only", with the external gamma given directly |

## DATA's conventions

These were found from the outputs. The test translates the inputs this way.

1. **Water** is entered as % of the wet mass. EPRdating uses the dry mass,
   so `w_dry = w_wet / (1 − w_wet)`.
2. **234U/238U** is the ratio of the incoming U. It is not back-corrected
   from a present-day value: `u234_u238_is="initial"`.
3. **Radon loss** applies to the dentine only. The enamel keeps its radon.
4. **Beta attenuation**: there is one factor per tissue for the whole U
   chain. It is applied to the dose rate with ingrowth. The factors are
   ROSY's one-group values corrected by Marsh (1999).
5. **Cosmic** dose rate is Prescott & Hutton (1994) for an overburden of
   2 g/cm³ at sea level. DATA has no site coordinates.
6. **Conversion factors** are Adamiec & Aitken (1998).

## Results

All figures are deviations of EPRdating from DATA, over the 82 runs, for
components where DATA prints at least 20 µGy/a.

| Component | `beta_by_segment=False` (DATA's convention) | Default (per segment) |
|---|---|---|
| Gamma + cosmic | −0.3 to +0.1 % | same |
| Internal (enamel α + β) | −2.7 to +2.6 % | −4.5 to +1.3 % |
| Beta from dentine | −1.4 to +5.7 % (mean +2.6 %) | +1.0 to +39.8 % (mean +18.7 %) |
| Beta from sediment | +7.1 to +13.0 % | same |
| **Age** | **−2.6 to +1.6 % (mean −0.1 %)** | −7.5 to +1.1 % (mean −1.3 %) |

DATA's printed beta factors compare as follows:

- **Dentine U**: EPRdating's one-group factors are DATA's to the third
  decimal, for example 0.134 vs 0.134 for 1000 µm with 50 µm removed. The
  only exceptions are 300 µm (+2.7 %) and 2500 µm or 200 µm removed (−2.8 %).
- **Sediment U and Th**: EPRdating is 0–8 % higher.
- **Sediment 40K**: EPRdating is 5–16 % higher. The gap grows with enamel
  thickness and causes most of the sediment beta difference.

Some example ages, in ka:

| Case | DATA EU | EPRdating EU, single factor | DATA LU | EPRdating LU, single factor | EPRdating LU, default |
|---|---|---|---|---|---|
| B0 | 76 | 76.2 | 93 | 93.7 | 92.5 |
| De 5 Gy | 4.4 | 4.4 | 5.0 | 5.0 | 4.9 |
| De 3000 Gy | 1787 | 1789 | 2380 | 2381 | 2380 |
| Dentine 50 ppm | 49 | 49.0 | 70 | 69.6 | 65.7 |
| Sediment on both sides | 73 | 72.1 | 81 | 79.6 | 79.7 |
| "Beta only" | 83 | 81.4 | 94 | 91.6 | 91.8 |

## Findings

### Beta dose with ingrowth (`beta_by_segment`)

Before 226Ra grows in, the U chain's beta dose comes mostly from 234mPa,
whose betas are hard (2.27 MeV). These betas cross the enamel better than
the average beta of the chain in equilibrium.

- **DATA** applies the equilibrium chain factor to the ingrowth-corrected
  dose rate. This underestimates the dentine beta dose of young teeth.
- **EPRdating** by default attenuates each U-series segment with its own
  factor, as ROSY does. This gives the following, compared with DATA:
  - up to +40 % at 5 ka (LU);
  - +12–20 % around 100 ka;
  - +1 % at 2 Ma, as the chain approaches equilibrium.
- The age effect is up to −7.5 %, for a tooth with 50 ppm U in the dentine.
- With `ToothSample(beta_by_segment=False)` EPRdating reproduces DATA. It
  should be used only to compare with published DATA ages.

### Sediment beta

The one-group sediment factors were checked against ROSY 2.0 within 2–5 %
(`ROSY_FINDINGS.md`). EPRdating's sediment beta is 7–13 % above DATA's,
mostly from 40K. DATA's sediment factors are therefore lower than the
one-group ones. Since the sediment beta is usually 5–15 % of the total dose
rate, the effect on ages is about 1–2 %. This is visible in the cases with
sediment on both sides and in the "beta only" case.

### A bug found and fixed

`ToothLayers.chain_fraction("sediment", …)` counted only the sediment on the
outer side of the enamel. A geometry without dentine, that is with sediment
on both sides, therefore received half its sediment beta dose. It now covers
both sides; `"sediment_outer"` and `"sediment_inner"` give each side alone.

- Teeth with thick dentine are unaffected, because the inner sediment
  contributes < 10⁻⁴ of the outer one. This includes all the ROSY and US-ESR
  validation cases.
- The bug was found because the cases with sediment on both sides came out
  5–22 % too old: +10 % for 1000 µm and +22 % for 300 µm.

## Tests

`test_data_reference.py` checks the following:

| Check | Tolerance |
|---|---|
| Gamma + cosmic | 0.5 % |
| Internal dose rate | 3.5 % |
| Dentine beta, single factor | −2.5 to +7 % |
| Dentine beta, per segment | > DATA, converging at 2 Ma |
| Sediment beta | +5 to +15 % |
| Beta factors | dentine within 4 %; sediment U/Th 0–8 %, K 3–17 % |
| Every age, single factor | 3 % (DATA prints 2–3 digits) |
| Ages, default | −8.5 to +2 % |

The inputs are in `data_reference/cases.json`, the parsed outputs in
`results.json`, and DATA's raw printouts in `raw/`.

# U-series/ESR: US-ESR and CS-US

DATA's <F6> ("read enamel into U-series") loads an `.EPR` file into a screen
that adds 234U/238U and 230Th/234U for the enamel and for each dentine side.
It prints the following:

- EU and LU results.
- **US-ESR**: age, the uptake parameter p of each tissue, and dose rates.
- **CS-US**: the closed-system U-series model (Grün 2000). DATA prints an age
  and, in a second column, the dose delivered by the tooth's own U since its
  uptake (Gy).

The U-series data cannot be stored in the file. `tools/data_harness/run_useries.py`
types them in and checks them against the echo on the printout.

There are 40 cases around the base tooth, with 230Th/234U from 0.05 to
0.95. They vary:

- 234U/238U (1.0–2.0, and different in each tissue);
- the D_E (30–1500 Gy);
- U contents, thickness and geometry;
- late uptake into the enamel, the dentine or both;
- the distance to the closed-system bound (p of the dentine −0.95 to −0.7).

## DATA's conventions

These are in addition to those of the EU/LU part:

- **Radon loss is ignored.** The case with 50 % Rn loss prints exactly the
  same as the base case.
- **CS-US uptake ages:** each tissue took up all its U at its closed-system
  age, computed from its own ratios.
- **EU and LU in the U-series screen:** the EU/LU columns change with the
  entered 230Th/234U. DATA takes each tissue's measured 234U/238U back to
  the initial ratio over that tissue's closed-system U-series age:
  1 + (r − 1)·exp(λ234·t_cs). This reproduces DATA's EU internal dose within
  2 % for every case. Two alternatives fail: reading the ratio as initial
  gives up to −8 %, and back-correcting over the ESR age gives up to +74 %.
  Published DATA EU ages computed in this screen (e.g. Lovedale, see
  `PUBLISHED_FINDINGS.md`) follow this convention.
- **No dentine U:** with sediment on both sides, DATA's US-ESR crashes with
  "Illegal function call in line 114 of module DATA-EL". CS-US is printed
  before the crash: 83 ka, against 82.0 ka in EPRdating.

## Results

These use EPRdating with `beta_by_segment=False` and no radon loss.

| | Agreement |
|---|---|
| US-ESR age (29 cases) | −1.4 to +0.7 %, mean 0.0 % |
| p of each tissue | within 0.05 (2 % for p > 2), from p = −0.8 to p = 16.7 |
| CS-US age (37 cases) | within the rounding to 1 ka plus 1 % |
| Dose from the tooth's U in CS-US | DATA 1–5 % lower (the dentine beta difference of the EU/LU part) |
| Uncertainties (4 cases, Monte Carlo n = 300) | standard deviation within 5 % of DATA's errors: for example 99 +12 −12 vs ± 12.2 ka, and 471 +62 −59 vs ± 58.8 ka |

**Same U-series model.** In the 29 regular cases, DATA's own age and p
values, put into EPRdating's U-series equations, give back the 230Th/234U
that was entered, within ±0.004. This holds from p = −0.8 to p = 16.7 and for
ages from 30 ka to 1 Ma.

With EPRdating's default (`beta_by_segment=True`), US-ESR ages are −2.1
to +1.0 % from DATA's (mean −1.7 %). For the U-rich dentine case (50 ppm)
they are 13 % younger. This is the same per-segment beta effect described above.

## Findings

### DATA's dentine p does not always converge

When the dentine took up its U much later than the enamel, DATA's dentine
p does not reproduce the measured dentine ratio:

| Case (enamel / dentine 230Th/234U) | DATA p dentine → predicted ratio | EPRdating p dentine | Age DATA / EPRdating |
|---|---|---|---|
| 0.30 / 0.10 | 1.25 → 0.229 | 6.69 | 99 / 104.1 ka |
| 0.25 / 0.05 | 0.00 → 0.348 | 16.9 | 100 / 108.0 ka |
| 0.25 / 0.15 | 2.46 → 0.182 | 3.60 | 104 / 105.1 ka |
| 0.25 / 0.15 and 0.35 | 2.49 → 0.175 | 3.37 | 100 / 100.7 ka |

DATA converges normally in the other configurations:

- dentine 0.20 against enamel 0.25 (p = 1.98);
- late enamel (p up to 16.7);
- both tissues late (p = 7.4 in each);
- large p in both tissues at a large D_E (p = 10.7 and 7.9).

The failure therefore depends on the dentine being much later than the
enamel, not on the size of p. In these cases DATA's US-ESR ages are up to
8 % too young, because its dentine dose is too high for the measured
230Th/234U. The CS-US ages of the same cases agree, which confirms that only
the p solution fails.

### No result near the closed-system bound

DATA prints "U-SERIES TOO HIGH: NO RESULT" whenever the solution needs a p
below a cutoff between −0.9 and −0.8.

- **DATA solves:** cases with dentine p = −0.8 and −0.7, in agreement with
  EPRdating.
- **DATA refuses:** cases with dentine p = −0.9, −0.95, −0.97 and −0.99.
  Solutions exist for all four: 42.7, 40.6, 76.6 and 9.2 ka.

EPRdating solves down to p = −1. With uncertain inputs, `age_mc` reports the
fraction of draws that solve and flags the result as marginal (see
`USESR_FINDINGS.md`).

### CS-US ages younger than the uptake

When the external dose alone accounts for the D_E before the U arrives, the
CS-US model has no solution. DATA prints an age anyway. With D_E = 30 Gy, for
example, it gives 20 ka, younger than the uptake ages it assumes itself
(31 and 38 ka). EPRdating returns `no_solution` and gives the uptake age.

## Tests

`test_data_useries.py` checks the following:

| Check | Tolerance |
|---|---|
| US-ESR ages | 1.5 % or 0.5 ka |
| p | 0.05 or 2 % |
| DATA's (T, p) reproduce the measured ratios | ±0.004 |
| Unconverged dentine p in DATA | identified, with EPRdating reproducing the measured ratio |
| No result near the bound | DATA has none; EPRdating solves with p < −0.85 |
| CS-US ages | 0.5 ka plus 1.5 % |
| Dose from the tooth's U | DATA/EPRdating ratio between 0.94 and 1.00 |
| Radon loss | ignored by DATA |

The inputs are in `data_reference/useries_cases.json`, the parsed outputs in
`useries_results.json`, and the raw printouts in `raw_useries/`.
