# Published ESR tooth-enamel case format (one JSON file per study)

File: <study_id>.json, e.g. "dirks2017_rising_star.json". Top level:

{
  "study_id": "dirks2017_rising_star",
  "reference": "full citation with DOI",
  "url": "URL(s) actually read",
  "software": "DATA" | "ROSY" | "USESR" | "other: <name>",   # program the authors used for ages
  "conventions": {                     # as stated in the paper; null if not stated
    "conversion_factors": "e.g. Adamiec & Aitken 1998 / Guerin 2011",
    "k_alpha": [0.13, 0.02],
    "beta_attenuation": "e.g. Brennan 2003 / Marsh 1999 / ROSY one-group",
    "water_basis": "wet" | "dry" | null,       # % of wet or dry mass
    "u234_u238_meaning": "measured present-day" | "initial" | null,
    "radon_loss": "...",
    "cosmic": "Prescott & Hutton 1994 ... (lat/long/alt/depth if given)",
    "gamma_source": "in situ NaI / HPGe on sediment / ..."
  },
  "samples": [
    {
      "sample": "label exactly as in the paper",
      "De_Gy": [value, error],
      "geometry": "e.g. SED/EN/DE, DE/EN/DE, SED/EN/SED (outer side first)",
      "enamel":  {"U_ppm": [v,e], "u234_u238": [v,e] or null, "th230_u234": [v,e] or null,
                  "thickness_um": [v,e], "removed_outer_um": [v,e], "removed_inner_um": [v,e], "density": null},
      "dentine": {"U_ppm": [v,e], "u234_u238": ..., "th230_u234": ..., "water_pct": [v,e], "radon_loss_pct": ...} or null,
      "dentine2_or_cementum": {... same keys ...} or null,
      "sediment": {"U_ppm": [v,e], "Th_ppm": [v,e], "K_pct": [v,e], "water_pct": [v,e]} or null,
      "gamma_uGy_a": [v,e] or null,      # if gamma (or gamma+cosmic) was measured or given directly
      "cosmic_uGy_a": [v,e] or null,
      "depth_m": ..., "lat_deg": ..., "lon_deg": ..., "altitude_m": ...,   # if given
      "published": {
        "dose_rates_uGy_a": {"internal": [v,e], "beta_dentine": [v,e], "beta_sediment": [v,e],
                             "gamma_cosmic": [v,e], "total": [v,e]},    # any subset, as printed
        "p_enamel": [v,e], "p_dentine": [v,e],
        "ages_ka": {"US": [v, +err, -err], "EU": [v,e], "LU": [v,e], "CSUS": [v,e]}   # any subset
      },
      "provenance": "Table 2 row 3; Table S4; page 7",
      "notes": "anything unusual"
    }
  ],
  "completeness": "complete" | "partial: <what is missing>",
  "extraction_notes": "how values were read, any ambiguity"
}

Rules:
- Copy numbers exactly as printed (units as printed converted only where obvious; say so in notes).
- Never guess or fill a value that is not in the source; use null.
- Every sample must have provenance (table / page).
- U in ppm (µg/g). Dose rates in µGy/a. Ages in ka. Thicknesses in µm.
