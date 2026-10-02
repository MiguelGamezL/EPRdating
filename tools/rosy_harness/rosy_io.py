"""Write ROSY 2.0 import files and parse ROSY output files (black-box I/O only)."""
import re

# order of the 70 numbers in the "DATA FOR INPUT" block (value, error pairs)
LAYERS = ("sed", "cem", "en", "den")
FIELDS = ["De", "ratio", "alpha_eff"]
for q in ("U", "Th", "K", "density", "radon"):
    FIELDS += [f"{q}_{l}" for l in LAYERS]
FIELDS += [f"water_{l}" for l in ("sed", "cem", "den")]
FIELDS += [f"thick_{l}" for l in ("cem", "en", "den")]
FIELDS += ["strip_out", "strip_in", "gamma_cosmic", "gamma", "depth", "overburden_density"]
assert len(FIELDS) == 35

DEFAULTS = dict(  # noqa: C408
    De=100.0, ratio=1.0, alpha_eff=0.15,
    U_sed=0, U_cem=0, U_en=0, U_den=0, Th_sed=0, Th_cem=0, Th_en=0, Th_den=0,
    K_sed=0, K_cem=0, K_en=0, K_den=0,
    density_sed=2.0, density_cem=2.54, density_en=3.0, density_den=2.82,
    radon_sed=1.0, radon_cem=1.0, radon_en=1.0, radon_den=1.0,
    water_sed=0, water_cem=0, water_den=0,
    thick_cem=0, thick_en=1000, thick_den=2000,
    strip_out=0, strip_in=0, gamma_cosmic=0, gamma=0, depth=0, overburden_density=2.13,
    uptake_cem=1, uptake_en=1, uptake_den=1, env_option=0,
    alpha_type="Varies with energy, Eref=5.3MeV", ratio_type="Initial Ratio",
)


def write_case(template_path, out_path, name, params):
    p = dict(DEFAULTS, **params)
    with open(template_path, newline="") as fh:
        t = fh.read()
    lines = t.split("\r\n")
    nums = []
    for f in FIELDS:
        v = p[f]
        val, err = (v if isinstance(v, (tuple, list)) else (v, 0.0))
        nums += [val, err]
    rows = ["  " + "   ".join(f"{x:10.4f}" for x in nums[i:i + 5]) for i in range(0, 70, 5)]
    i = lines.index("DATA FOR INPUT")
    lines[i + 1:i + 15] = rows
    lines[i + 15] = f"{p['uptake_cem']}  {p['uptake_en']}  {p['uptake_den']}  {p['env_option']}  0  "
    lines = [re.sub(r"(NAME OF SAMPLE:\s+).*", lambda m: m.group(1) + name, l) for l in lines]
    lines = [re.sub(r"(RATIO OF U234 TO U238:.*\d\s+)(Initial Ratio|Present Ratio)", lambda m: m.group(1) + p["ratio_type"], l) for l in lines]
    lines = [re.sub(r"(ALPHA EFFICIENCY:.*\d\s+)(Constant|Varies with energy, Eref=5.3MeV)", lambda m: m.group(1) + p["alpha_type"], l) for l in lines]
    j = next(k for k, l in enumerate(lines) if l.startswith(("(Dry)", "(Wet)")))
    lines[j] = "(Dry)   " + p["alpha_type"]
    if "trailer" in p:
        lines[j + 1] = "  " + "   ".join(f"{x:10.4f}" for x in p["trailer"])
    if "flag5" in p:
        lines[i + 15] = f"{p['uptake_cem']}  {p['uptake_en']}  {p['uptake_den']}  {p['env_option']}  {p['flag5']}  "
    with open(out_path, "w", newline="") as fh:
        fh.write("\r\n".join(lines))


def parse_output(path):
    with open(path, newline="") as fh:
        t = fh.read().replace("\r", "")
    out = {}
    blocks = re.split(r"\*\*\*\s+(EARLY|LINEAR|COMBINATION) UPTAKE\s+\*\*\*", t)
    for kind, body in zip(blocks[1::2], blocks[2::2]):
        key = {"EARLY": "EU", "LINEAR": "LU", "COMBINATION": "CU"}[kind]
        m = re.search(r"Age Estimate \(years\) =\s*([-\d.Ee+]+)\s*\+/-\s*([-\d.Ee+]+)", body)
        d = {"age": float(m.group(1)), "age_err": float(m.group(2))} if m else {"age": None}
        for row in ("Dentine", "Enamel", "Cementum", "Sediment", "Uranium", "Thorium", "Potassium", "Total"):
            m = re.search(rf"^\s*{row}\s+(.*)$", body, re.MULTILINE)
            if m:
                cells = m.group(1).split()
                d[row] = [None if c == "-" else float(c) for c in cells]
        m = re.search(r"Total dose rate \(microGy/a\) =\s*([-\d.]+)", body)
        d["total_rate"] = float(m.group(1)) if m else None
        out[key] = d
    return out
