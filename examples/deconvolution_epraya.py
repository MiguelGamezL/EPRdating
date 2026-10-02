"""Deconvolution of enamel spectra with component shapes simulated in EPRAYA.

Requires: pip install 'eprdating[spectra]'

Replace `load_series()` with your own reader. The g-tensors below are
placeholders: set them from the literature for each CO2- species you
want to separate.
"""


from eprdating import fit_dose_response
from eprdating.spectra import component_vs_dose
from eprdating.spectra.epraya_backend import Species, basis_from_species

FREQ_GHZ = 9.5


def load_series():
    """Return (doses, B_mT, list_of_spectra). Placeholder."""
    raise NotImplementedError("read your spectra here")


doses, B, spectra = load_series()

species = {
    "CO2- orthorhombic": Species(g=[2.0031, 1.9973, 2.0019], Hpp=[0, 0.25]),  # placeholder values
    "CO2- axial": Species(g=[2.0032, 2.0032, 1.9970], Hpp=[0, 0.25]),  # placeholder values
    "native": Species(g=2.0045, Hpp=[0, 0.8]),  # placeholder values
}
basis = basis_from_species(B, species, FREQ_GHZ)
amps, errs, fits = component_vs_dose(basis, spectra, "CO2- orthorhombic")
drc = fit_dose_response(doses, amps, model="SSE", weighting="1/I^2")
print(drc.summary())
print("min R² of spectral fits:", min(f.r2 for f in fits))
