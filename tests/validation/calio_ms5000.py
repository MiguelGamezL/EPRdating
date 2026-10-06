"""De of the Calio P4 enamel powder (Hakim et al. 2025, Nature) from its public spectra.

Data: Zenodo 15515771, ``Calio ESR spectra.zip`` (CC-BY 4.0). The spectra used
here are in ``tests/data/calio_p4`` (see its README); the full archive works too.

    python tests/validation/calio_ms5000.py ["path/to/Calio ESR spectra"]

Prints the peak-to-peak intensities (40 G window around g = 2.000, the three
rotations of each aliquot averaged) and SSE fits with the doses as labelled
in the archive and with the 8000 and 15 000 Gy labels exchanged. Published:
SSE De = 2267 ± 99 Gy. See docs/validation.md.
"""

import re
import sys
import warnings
from pathlib import Path

import numpy as np

from eprdating import fit_dose_response
from eprdating.spectra import IntensityWindow, combined_intensity, empirical_template, read_epr

PUBLISHED = (2267.0, 99.0)
DATA = Path(__file__).resolve().parents[1] / "data" / "calio_p4"
WINDOW = IntensityWindow(40, center_g=2.000)


def load(folder: Path) -> dict[float, list]:
    """Spectra by dose: the ``_result.xml`` of each rotation, or its runs."""
    groups: dict[float, list] = {}
    for d in sorted(folder.iterdir()):
        m = re.fullmatch(r"Calio_(?:(\d+)Gy_\d|N(?:_\d)?)", d.name)
        if not (d.is_dir() and m):
            continue
        r = d / f"{d.name}_result.xml"
        groups.setdefault(float(m.group(1) or 0), []).append(read_epr(r if r.exists() else d))
    return groups


def main(folder: Path) -> None:
    warnings.simplefilter("ignore")
    groups = load(folder)
    doses = np.array(sorted(groups))
    swapped = doses.copy()
    swapped[doses == 8000], swapped[doses == 15000] = 15000, 8000
    template = empirical_template([s for v in groups.values() for s in v], WINDOW)
    for method in ("peak_to_peak", "t1_b2", "template", "double_integral"):
        rs = [combined_intensity(groups[d], method, template if method == "template" else None, WINDOW,
                                 max_shift=0.3) for d in doses]
        I = np.array([r.value for r in rs])
        S = np.array([r.sigma for r in rs])
        print(f"\n{method}: " + "  ".join(f"{d:g}:{i:.0f}" for d, i in zip(doses, I, strict=True)))
        for label, D in (("as labelled", doses), ("8000 <-> 15000", swapped)):
            for wl, sig in (("1/I^2", 0.05 * I), ("data", S)):
                f = fit_dose_response(D, I, "SSE", sigma=sig)
                print(f"  SSE {label:15s} weights {wl:6s} De = {f.De:6.0f} ± {f.De_sigma:4.0f} Gy  "
                      f"(chi2_red {f.chi2_red:.2f})")
    print(f"\npublished (SSE): {PUBLISHED[0]:.0f} ± {PUBLISHED[1]:.0f} Gy")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else DATA)
