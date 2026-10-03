# Validation against DATA (Grün 2009)

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
