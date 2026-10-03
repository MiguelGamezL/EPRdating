# Recomputing published ESR tooth ages

DRAC (Durcan et al. 2015) was validated by recomputing dose rates from
published studies. The same is done here for ESR dating of tooth enamel.
Published inputs are transcribed and fed to EPRdating, and the result is
compared with the published age. The inputs are translated with the
conventions of the program each study used (DATA, USESR or ROSY), since the
aim is to show that EPRdating reproduces the published calculation.

## Data

`published/*.json` holds one file per study, in the format of
`published/SCHEMA.md`. Each file records:

- the inputs, with their errors;
- the published dose rates, p values and ages;
- the conventions the paper states;
- the provenance (table and page) of every sample.

The tables were transcribed from open-access papers, supplements and
reports, and each table was read twice. A separate check re-read the key
values from the original sources without seeing the transcription. It found
no copy errors in the numbers, but it did correct some of the context:

- the side each removed thickness belongs to;
- that the thickness is the initial one;
- the cementum geometry of one Khok Sung tooth;
- a cosmic dose rate given only in the text.

**Complete studies, 58 published ages:**

| Study | Program | Ages | Site |
|---|---|---|---|
| Richard et al. 2023 | DATA | 4 EU | Lovedale, South Africa |
| Duval et al. 2019 | DATA (cited; not named) | 5 US-ESR (converged dentine p) + 11 CS-US | Khok Sung, Thailand |
| Bahain et al. 2020, CENIEH-RSES set | DATA (as given) | 2 US + 2 CS-US | Tourville, France |
| Falguères et al. 2025 | DATA | 7 US | Fumane, Italy |
| Bahain et al. 2020, MNHN set | USESR (US), DATA (CS-US) | 5 US + 5 CS-US | Tourville, France |
| Dirks et al. 2017 | USESR | 12 US (4 teeth, 2 laboratories, 2 radon scenarios) | Rising Star, South Africa |
| Rizal et al. 2020 | USESR | 3 US | Ngandong, Indonesia |
| Duval & Zhao 2021 | USESR | 2 US | Olieboomspoort, South Africa |
| Carvajal, Montes & Almanza 2011 | ROSY | EU = LU (external dose only) | Aguazuque, Colombia |

**Set apart, and why:**

- **Khok Sung US-ESR, 9 ages:** the published dentine p does not reproduce
  the measured dentine 230Th/234U (see below).
- **Tourville MNHN, TVL 929(a):** the paper prints 3.046 ppm U in the
  dentine, where the other teeth have 20–31 ppm. Its own Dβ needs about
  15–20 ppm, so the value is an error in the source.
- **Tourville CENIEH, T2:** no solution at the nominal inputs. The published
  p is −0.93 in both tissues, at the closed-system bound; with uncertain
  inputs, `age_mc` would flag it as marginal.
- **Yu et al. 2024:** internally inconsistent. The D_E divided by the
  printed total dose rate does not give the printed age (2103 Gy /
  698 µGy/a = 3013 ka, against 3264 ka printed), and the methods and the
  tables give different water contents.
- **Seven ROSY studies** (Kinoshita et al. 2008, 2011, 2014; Ribeiro et al.
  2013, 2021; Oliveira et al. 2010; Kerber et al. 2011; Qi et al. 2018):
  they give no enamel thickness. Several also lack water, k or 234U/238U,
  or have enamel U below detection. ROSY's one-group beta attenuation
  cannot be recomputed without the thickness. The reason for each is
  recorded in its file (`benchmark_reason`).
- **Models EPRdating does not compute** are not compared: AU-ESR (Fumane),
  ages combined over laboratories (Rising Star), and minimum ages.
- **Not reachable here:** Jebel Irhoud, Lida Ajer, Misliya, Mala Balanica,
  Callabonna, the Laos molar and El Mnasra. Their supplements were images,
  .doc files, or blocked.

## Conventions applied

- **DATA and USESR:** one beta attenuation factor for the U chain
  (`beta_by_segment=False`). Water is per wet mass, which is DATA's input.
- **ROSY:** beta attenuated per segment, water per dry mass.
- **Conversion factors** as stated in each paper; Adamiec & Aitken (1998)
  when not stated.
- **Dentine of unstated thickness** is thick (5 mm), as DATA assumes.
- **Removed thicknesses** go on the side the paper names. When the paper
  numbers the sides without naming them, side 1 is taken as the outer one.
- **Thickness:** the printed enamel thickness is the initial one, as all the
  checked papers state.
- **Radon loss:**
  - in the sediment, as a 226Ra deficit;
  - per tissue from 222Rn/230Th (Tourville MNHN);
  - none in DATA's US-ESR, which ignores it.
- **234U/238U:**
  - US-ESR and CS-US use the measured ratios.
  - DATA's EU/LU ages computed in its U-series screen use the measured ratio
    taken back to the initial one over each tissue's closed-system U-series
    age. This was found with the DATA harness and is tested in
    `test_data_useries.py`.
- **Gamma and cosmic:** where a paper gives no separate values, its
  published gamma+cosmic dose rate (or sediment beta+gamma+cosmic for
  Fumane) is used as an input. That component is then not recomputed.

## Results

![EPRdating vs published ages](../../docs/img/published_ages.png)

| | |
|---|---|
| Ages compared | 58, from 9 studies (DATA 29, USESR 27, ROSY 2), 3 ka – 720 ka |
| EPRdating − published | −9.5 to +3.4 %, mean −3.2 % |
| Within 5 % | 44 of 58 |
| Within the published 1σ | 56 of 58. The two exceptions are just outside: −6.7 % against ±5 %, and −9.5 % against ±9 %. |

| Study | Program | EPRdating − published |
|---|---|---|
| Lovedale | DATA, EU | −3.6 to −0.9 % |
| Khok Sung, CS-US and converged US-ESR | DATA | −5.8 to +0.6 % (mean −2.2 %) |
| Tourville CENIEH | DATA | −3.8 to −3.2 % |
| Fumane | DATA | −9.5 to −0.1 % (mean −6.5 %) |
| Tourville MNHN | USESR / DATA | −6.7 to +3.0 % (mean −3.0 %) |
| Rising Star | USESR | −7.5 to −1.3 % (mean −4.3 %) |
| Ngandong | USESR | −6.4 to −3.1 % |
| Olieboomspoort | USESR | −0.8 % |
| Aguazuque | ROSY, external only | +3.4 % |

## Findings

### 1. A small systematic offset from the beta doses

EPRdating's ages are on average 3 % younger. The published dose-rate
components show where the difference comes from:

- **Internal dose and gamma:**
  - Lovedale (DATA): internal −1 to −2 %.
  - Aguazuque (ROSY): gamma −2 %.
  - USESR studies: EPRdating's internal dose is higher, +1 to +21 %
    (+9 to +15 % at Ngandong).
- **Beta doses are higher in EPRdating:**
  - **Dentine beta** is 3–5 % higher at Lovedale (DATA), 7–15 % at
    Ngandong, 13–43 % at Fumane and 19–49 % at Rising Star.
  - **Sediment beta** is about 10 % higher. This matches the synthetic DATA
    comparison (`DATA_FINDINGS.md`), where the difference comes mostly from
    40K, and the USESR comparison (`USESR_FINDINGS.md`).

The one-group factors were themselves checked against ROSY 2.0 within 2–5 %.
These differences therefore reflect how DATA and USESR attenuate beta rays,
not an error in EPRdating's attenuation.

### 2. DATA's dentine p fails in a published study

Duval et al. (2019) give p values for every tissue. Putting their published
age and p back into the U-series equations shows two things:

- **Every enamel p** gives back the measured 230Th/234U within ±0.003.
- **The dentine p fails** for all nine US-ESR ages where the dentine is much
  younger than the enamel (dentine 230Th/234U 0.12–0.13 against enamel
  0.35–0.48). The predicted ratio is 0.22–0.45 against the 0.12–0.13
  measured. The tooth whose dentine is not late (3548) converges.

This is the behaviour found with the DATA harness. Where the dentine carries
U (3545, 3547, 3549), the published US-ESR ages are 5–13 % younger than the
ages that fit the measured ratios. For 3546A/B the dentine holds 2 ppm U, so
the wrong p hardly changes the age (−2.5 to −0.1 %). The CS-US ages of the
same teeth, which do not depend on p, agree within 2 %. EPRdating's US-ESR
ages for 3547 and 3549 lie between the published US-ESR and CS-US ages, as
expected.

### 3. Radon loss in the sediment

Rising Star gives two scenarios: 80 % radon loss in the sediment, and none.
The two scenarios bracket the comparison, with the sediment beta behaving as
follows:

| Scenario | EPRdating vs published sediment beta |
|---|---|
| No radon loss | −2 to +10 % |
| 80 % radon loss | +9 to +28 % |

Going from no loss to 80 % loss, EPRdating's sediment beta falls by about
10 % and the published USESR value by about 23 %. EPRdating applies the loss
to 226Ra and its daughters, about 59 % of the U-chain beta. The published
drop corresponds instead to removing 80 % of all the U-chain beta.

### 4. Aguazuque (ROSY)

The tooth's U was not measured and is taken as zero, so only the external
dose is compared:

- **Gamma:** 340 µGy/a against 347 µGy/a published.
- **Cosmic:** as given.
- **Sediment beta:** 33 µGy/a against 47 µGy/a. The published value
  corresponds to a sediment attenuation factor of 0.144. One-group theory
  gives that factor for enamel about 0.95 mm thick, not 1.4 mm.
- **Age:** +3.4 %.

The soil table prints "40K 1.66 ± 0.09 ppm" (by gamma spectrometry), and the
dose table has an explicit potassium row of zeros. If 1.66 ppm refers to
40K, the soil holds about 1.42 % K. With that K, the dose rate would rise
from 645 to about 1020 µGy/a and the age would fall from 3.3 ka to about
2.1 ka. Whether potassium was left out on purpose is worth checking with the
original data.

### 5. What a paper must give for its ages to be recomputable

Only 9 of the 18 studies collected could be recomputed. To recompute a
tooth age, a paper must give:

- the D_E;
- U, and for US-ESR the 234U/238U and 230Th/234U, of each tissue;
- the initial enamel thickness and the thickness removed from each side,
  naming which side;
- the geometry (which tissue or sediment is on each side);
- the sediment U/Th/K and water content, stating whether the water is per
  wet or dry mass;
- the gamma and cosmic dose rates separately;
- the radon loss;
- the program and its conversion factors.

Most of the older ROSY literature omits the thickness, and none of the
studies states the water basis. An open program that records every input
with the age, such as EPRdating's per-sample JSON, removes this problem.

## Files

| File | Content |
|---|---|
| `published/` | study files and the format |
| `published_cases.py` | conventions and recomputation |
| `test_published_ages.py` | tests |
| `tools/plot_published_ages.py` | figure |
