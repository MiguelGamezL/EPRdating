# Validation of US-ESR (`eprdating.usesr`) against published results

Model: Grün, Schwarcz & Chadam (1988), as formalised by Shao et al. (2015):
U(t) = U_m (t/T)^(p+1), uranium arriving with a constant 234U/238U (r0) and
no 230Th, no Th uptake and no isotope loss. The p–T relation of each tissue
follows from its measured 234U/238U and 230Th/234U; the age is where the ESR
dose equation also holds.

## U-series part (exact check)

| Case | Tissue | p at published age (EPRdating) | Published p |
|---|---|---|---|
| Shao et al. 2015, T = 726 ka | enamel | 0.49 | 0.50 ± 0.33 |
| | dentine | 0.32 | 0.32 ± 0.28 |
| | cementum | 0.32 | 0.33 ± 0.27 |
| De Nadale CN-857, T = 79 ka | enamel / dentine | −0.49 / +1.35 | −0.49 / +1.36 |
| De Nadale CN-1322, T = 72 ka | enamel / dentine | −0.25 / +0.09 | −0.27 / +0.07 |
| De Nadale CN-55, T = 55 ka | enamel / dentine | (see test) | −0.90 / −0.86 |

The residual differences (≤ 0.02) are what the rounding of the published
ages (1 ka) produces. The derivation in Shao et al. (eqs. 15–19) is
algebraically identical to the parcel average used here.

## Dose rates and ages (convention-limited)

At the published age and p values, EPRdating's internal (enamel) dose rate
matches USESR within 2–3 % (CN-857: 142 vs 145 µGy/a; CN-55: 188 vs 192).
Beta doses from dentine and sediment are 8–23 % higher in EPRdating:

| De Nadale | Sediment β (EPRdating / USESR) | Dentine β |
|---|---|---|
| CN-857 | 92 / 79 | 57 / 54 |
| CN-1322 | 84 / 78 | 165 / 134 |
| CN-55 | 77 / 71 | 145 / 118 |

EPRdating's sediment β for these geometries agrees with the factors measured
from ROSY 2.0 (one-group, Brennan et al. 1997), so the difference lies between
the USESR and ROSY implementations of beta attenuation (and the unreported
dentine thickness, 2000 µm assumed), not in the US-ESR algebra.

Resulting ages:

| Sample | EPRdating (nominal; MC) | Published |
|---|---|---|
| CN-857 | 77.1 ka; 77.3 ± 7.5 ka | 79 ± 9 ka |
| CN-1322 | 65.8 ka; 67.6 ± 11 ka | 72 ± 11 ka |
| CN-55 | no solution; 59.3 ± 7.9 ka, **marginal** (50 % of draws solve) | 55 ± 11 ka |
| Shao et al. 2015 example (cementum 1000 µm assumed) | 689 ka | 726 +87/−84 ka |

CN-55 sits at the closed-system U-series bound (p ≈ −0.9 in both tissues) and
its authors flag it as unreliable (younger than the U/Th minimum age of the
unit). Its closed-system ages are 49.0 ka (enamel) and 47.5 ka (dentine). At
49 ka the one-group doses give 41.9 Gy against De = 41.7 Gy: the dose is
reached just before the bound (it would need ≤ 851 µGy/a; EPRdating gives 856,
USESR 798), so the nominal inputs have no solution.

`USESRSample.age_mc` always runs and reports the fraction of draws that admit
a solution; below 80 % (`marginal_below`) the result is flagged as marginal
and the age is conditional on the solving draws. For CN-55 (1000 draws):
"59.26 ± 7.9 ka (68 %: 51.65–67.12; solution in 50 % of 1000 draws) —
WARNING marginal", consistent with the published 55 ± 11 ka and with the
authors' caveat.

The Aves and Milo's samples (Yu et al. 2024, PeerJ 12:e17478) have 230Th/234U
at or above secular equilibrium, which leaves p essentially unconstrained;
they are not used.

Tests: `test_usesr_published.py`.

## CSUS-ESR (Grün 2000)

`USESRSample.age(model="CSUS")` takes all the U of each tissue as taken up
at once at its closed-system U-series age, the alternative model DATA prints
next to US-ESR. Comparing both shows how much an age depends on the uptake
model. For the De Nadale cave teeth (Yu et al. 2026, J. Hum. Evol. 215–216, 103842):

| Sample | US-ESR | CSUS-ESR | Uptake ages (enamel / dentine) | Published (US) |
|---|---|---|---|---|
| CN-857 | 77.1 ka | 80.5 ka | 49.6 / 22.3 ka | 79 ± 9 ka |
| CN-1322 | 65.8 ka | 68.2 ka | 39.1 / 32.6 ka | 72 ± 11 ka |
| CN-55 | no solution | no solution | 49.0 / 47.5 ka | 55 ± 11 ka |

CSUS-ESR ages are 3–4 % older here: uptake as a late step delivers less
internal dose than the gradual uptake US-ESR finds. The model is checked by
a round trip in `tests/test_usesr.py`.

## Against DATA

Both models were also run against DATA's U-series/ESR part (Grün 2009) on 40
cases. US-ESR ages agree within −1.4 to +0.7 %, p within 0.05, and CS-US
ages within DATA's rounding; see the second half of `DATA_FINDINGS.md`.
