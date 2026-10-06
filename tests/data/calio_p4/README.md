# Calio P4 enamel powder: MS5000 spectra

Source: Hakim, B. et al. (2025) Hominins on Sulawesi during the Early
Pleistocene. *Nature*, doi:10.1038/s41586-025-09348-6. Data: Zenodo record
15515771, `Calio ESR spectra.zip`, licence CC-BY 4.0.

Fossil *Celebochoerus* P4 tooth enamel, powder (90–180 µm, ~95–100 mg
aliquots): natural and ⁶⁰Co gamma doses of 50, 100, 250, 600, 1200, 2400,
4000, 8000 and 15 000 Gy; each tube measured in up to three rotations (120°).
Freiberg MS5000, 2 mW, 0.1 mT modulation, 12 mT sweep, 10 accumulations per
run. Published De (T1–B2 peak-to-peak, SSE): 2267 ± 99 Gy.

What is here: for each rotation (`Calio_<dose>Gy_<rotation>`, `Calio_N`,
`Calio_N_2`) the `*_result.xml` file (average of the ten runs), or the ten
runs when the archive has no result file (`Calio_100Gy_1`). Folder and file
names are those of the archive. To keep the repository small, the raw
quadrature curves (`MWAbsorption sinus` / `cosinus`) were removed from each
`.xml`; the `BField` and `MWAbsorption` curves, the measurement attributes and
the recipe are unchanged, so the spectra read exactly as from the original
files. Files not used (individual runs when a result exists, duplicated
nested folders, `.csv` exports) are left out.

Note: the 8000 Gy aliquot is 24 % more intense than the 15 000 Gy one; with
the two labels exchanged the published De is reproduced (docs/validation.md).
The labels are kept as in the archive.
