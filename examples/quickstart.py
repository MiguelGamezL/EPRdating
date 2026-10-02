"""Quick start: from a dose-response series to a Monte Carlo age.

All numbers below are illustrative, not from a real sample.
"""

import numpy as np

from eprdating import (
    LinearUptake,
    Sediment,
    ToothLayers,
    ToothSample,
    USModel,
    cosmic_dose_rate,
    fit_dose_response,
)

# 1. Dose-response curve (added dose in Gy, T1-B2 intensity in a.u.)
dose = np.array([0, 50, 100, 200, 400, 800, 1600, 3200, 6400])
intensity = np.array([102.0, 126.0, 151.0, 196.0, 279.0, 410.0, 590.0, 780.0, 905.0])
drc = fit_dose_response(dose, intensity, model="SSE", weighting="1/I^2")
print(drc.summary(), "\n")

# 2. Cosmic dose rate for a site near Bogotá (2600 m a.s.l., 1.5 m deep)
cosmic = cosmic_dose_rate(depth_m=1.5, density=1.9, lat_deg=4.6, lon_deg=-74.1, altitude_m=2600)

# 3. Sample description
sample = ToothSample(
    De=(drc.De, drc.De_sigma),
    enamel_U=(0.8, 0.08),
    dentine_U=(15.0, 1.5),
    sediment=Sediment(U=(2.1, 0.1), Th=(7.5, 0.4), K=(1.1, 0.05), water=(0.15, 0.05)),
    # beta attenuation by one-group theory (as in ROSY)
    beta=ToothLayers(enamel_um=1100, dentine_um=2000, strip_outer_um=50, strip_inner_um=50),
    cosmic=cosmic,
    uptake_enamel=LinearUptake(),
    uptake_dentine=USModel(p=0.5),
    u234_u238_dentine=(1.25, 0.02),  # measured today
    radon_loss_dentine=(0.3, 0.1),
)

print(sample.age().summary(), "\n")
print(sample.age_mc(n=2000, seed=42).summary())
