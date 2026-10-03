"""Build EPRdating samples from published ESR tooth data (``published/*.json``).

Each study file follows ``published/SCHEMA.md``. Inputs are translated with
the conventions of the program the authors used:

- DATA (Grün 2009) and USESR (Shao et al. 2014): one beta attenuation factor
  for the whole U chain (``beta_by_segment=False``); water as % of the wet
  mass (DATA's input convention; USESR is assumed to follow it);
- ROSY (Brennan et al. 1999): beta attenuated per U-series segment (default);
- conversion factors as stated by the authors (Adamiec & Aitken 1998 when
  not stated, the default of both DATA and USESR before 2011);
- dentine of unstated thickness is taken as thick (5 mm), as DATA does;
- radon loss in the sediment is applied as a deficit of 226Ra and its
  daughters (``Sediment.U_ra226``);
- when a paper gives no separate gamma and cosmic dose rates, the published
  gamma+cosmic (or external) dose rate is used as an input and that
  component is not recomputed (``Recalc.external_from_paper``).
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field, replace
from pathlib import Path

from eprdating import Sediment, ToothLayers, ToothSample, USModel
from eprdating._types import as_value
from eprdating.dose_rate import cosmic_dose_rate
from eprdating.series import LAMBDA
from eprdating.usesr import UseriesData, USESRSample, closed_system_age

DIR = Path(__file__).parent / "published"
UPTAKE = {"EU": -1.0, "LU": 0.0}
#: models EPRdating recomputes (others, e.g. AU-ESR or ages combined over
#: laboratories, are listed but not compared)
MODELS = {"EU", "LU", "US", "CSUS"}


def load_studies() -> dict[str, dict]:
    return {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted(DIR.glob("*.json"))}


def _v(x, default=0.0):
    """Central value of [v, e] / number / None."""
    if x is None:
        return default
    if isinstance(x, (list, tuple)):
        return default if x[0] is None else float(x[0])
    return float(x)


def _ve(x):
    """(value, error) tuple, error 0 when missing."""
    if x is None:
        return None
    if isinstance(x, (list, tuple)):
        return (float(x[0]), float(x[1]) if len(x) > 1 and x[1] is not None else 0.0)
    return (float(x), 0.0)


def factors_of(study: dict) -> str:
    text = (study.get("conventions") or {}).get("conversion_factors") or ""
    if re.search(r"Gu[eé]rin", text) and not re.search(r"Adamiec.*DATA", text):
        return "guerin_2011"
    if re.search(r"Liritzis", text):
        return "liritzis_2013"
    return "adamiec_aitken_1998"


def program_of(study: dict) -> str:
    sw = study.get("software", "")
    for name in ("ROSY", "USESR", "DATA"):
        if name in sw:
            return name
    return "DATA"  # Khok Sung: DATA cited, program not named


@dataclass
class Recalc:
    study: str
    sample: str
    program: str
    tooth: ToothSample
    us: USESRSample | None
    published: dict
    external_from_paper: bool = False
    notes: list[str] = field(default_factory=list)


def _th(t: dict):
    """230Th/234U, also from a printed 230Th/238U."""
    if t.get("th230_u234") is not None:
        return t["th230_u234"]
    if t.get("th230_u238") is not None and t.get("u234_u238") is not None:
        return [_v(t["th230_u238"]) / _v(t["u234_u238"]), None]
    return None


def _removed(en: dict):
    """(outer, inner) removed thickness; papers that number the sides without
    naming them give side 1 as the outer one."""
    out, inn = en.get("removed_outer_um"), en.get("removed_inner_um")
    if out is None and inn is None:
        out, inn = en.get("removed_side1_um"), en.get("removed_side2_um")
    return _v(out), _v(inn)


def _r0(t: dict, prog: str, conv: dict) -> float:
    """234U/238U for EU/LU ages. DATA's U-series screen turns the measured
    ratio into the initial one with the tissue's closed-system U-series age
    (found with the DATA harness, see DATA_FINDINGS.md); elsewhere the ratio
    is used as entered."""
    r = _v(t.get("u234_u238"), 1.0)
    th = _th(t)
    if prog == "DATA" and "U-series screen" in (conv.get("u234_u238_meaning") or "") and th is not None:
        tc = closed_system_age(_v(th), r)
        if math.isfinite(tc):
            return 1.0 + (r - 1.0) * math.exp(LAMBDA["U234"] * tc)
    return r


def _water(pct, basis: str) -> float:
    w = _v(pct) / 100.0
    if w <= 0:
        return 0.0
    return w / (1.0 - w) if basis == "wet" else w


def build(study: dict, s: dict, beta_by_segment: bool | None = None, water_basis: str | None = None) -> Recalc:
    prog = program_of(study)
    conv = study.get("conventions") or {}
    # ROSY takes water per dry mass (checked against ROSY 2.0), DATA per wet mass
    basis = water_basis or conv.get("water_basis") or ("dry" if prog == "ROSY" else "wet")
    seg = (prog == "ROSY") if beta_by_segment is None else beta_by_segment
    en, de = s["enamel"], s.get("dentine") or {}
    c2 = s.get("dentine2_or_cementum") or {}
    sed = s.get("sediment") or {}
    geom = (s.get("geometry") or "").upper()
    notes = []

    # a tissue counts as present when the paper describes it, with or without U data
    inner_dentine = bool(de) and not geom.endswith("/SED")
    outer_tissue = bool(c2)
    if geom.startswith(("DE", "CEM", "T2")) and not outer_tissue:
        notes.append("outer tissue without U data: treated as sediment")
    layers = ToothLayers(
        enamel_um=_v(en.get("thickness_um")),
        dentine_um=_v(de.get("thickness_um"), 5000.0) if inner_dentine else 0.0,
        cementum_um=_v(c2.get("thickness_um"), 5000.0) if outer_tissue else 0.0,
        strip_outer_um=_removed(en)[0],
        strip_inner_um=_removed(en)[1],
        enamel_density=_v(en.get("density"), 2.95),
        dentine_density=2.82, cementum_density=2.82,
    )
    rn_sed = _v(sed.get("radon_loss_pct")) / 100.0
    U_s = _v(sed.get("U_ppm"))
    sediment = Sediment(U=U_s, Th=_v(sed.get("Th_ppm")), K=_v(sed.get("K_pct")),
                        water=_water(sed.get("water_pct"), basis),
                        U_ra226=U_s * (1 - rn_sed) if rn_sed > 0 else None)

    pub = s.get("published") or {}
    dr = pub.get("dose_rates_uGy_a") or {}
    gamma = _v(s.get("gamma_uGy_a"), None)
    cosmic = _v(s.get("cosmic_uGy_a"), None)
    external_from_paper = False
    if gamma is not None and cosmic is None:
        if s.get("depth_m") is not None and s.get("lat_deg") is not None:
            cosmic = 1000 * as_value(cosmic_dose_rate(depth_m=_v(s["depth_m"]), density=2.0, lat_deg=_v(s["lat_deg"]),
                                                      lon_deg=_v(s.get("lon_deg")),
                                                      altitude_m=_v(s.get("altitude_m")))).value
            notes.append(f"cosmic from depth and site: {cosmic:.0f} µGy/a")
        elif dr.get("gamma_cosmic") is not None:
            cosmic = _v(dr["gamma_cosmic"]) - gamma
        else:
            cosmic = 0.0
    ext_with_beta = None
    gamma_from_sediment = re.search(r"computed|calculated", conv.get("gamma_source") or "", re.IGNORECASE) and sed \
        and s.get("gamma_uGy_a") is None and cosmic is not None
    if gamma is None and gamma_from_sediment:
        notes.append("gamma computed from the sediment, as in the paper")
    elif gamma is None:
        gc = dr.get("gamma_cosmic")
        if gc is not None:
            gamma, cosmic = _v(gc), 0.0
            external_from_paper = True
            notes.append("gamma+cosmic taken from the paper")
        elif dr.get("beta_gamma_sediment_plus_cosmic") is not None:
            ext_with_beta = _v(dr["beta_gamma_sediment_plus_cosmic"])
            gamma, cosmic = 0.0, 0.0
            external_from_paper = True
            notes.append("sediment beta + gamma + cosmic taken from the paper")
    k = _v(conv.get("k_alpha"), 0.13) if conv.get("k_alpha") else 0.13

    rn_d = _v(de.get("radon_loss_pct")) / 100.0
    if not rn_d and de.get("rn222_th230") is not None:
        rn_d = 1.0 - _v(de["rn222_th230"])
    rn_e = 1.0 - _v(en["rn222_th230"]) if en.get("rn222_th230") is not None else 0.0
    tooth = ToothSample(
        De=_v(s["De_Gy"]),
        enamel_U=_v(en.get("U_ppm")),
        dentine_U=_v(de.get("U_ppm")) if inner_dentine else 0.0,
        cementum_U=_v(c2.get("U_ppm")) if outer_tissue else 0.0,
        sediment=sediment,
        beta=layers,
        gamma=None if gamma is None else gamma / 1000.0,
        cosmic=0.0 if cosmic is None else cosmic / 1000.0,
        k_alpha=k,
        dentine_water=_water(de.get("water_pct"), basis),
        cementum_water=_water((c2 or de).get("water_pct"), basis),
        u234_u238_enamel=_r0(en, prog, conv),
        u234_u238_dentine=_r0(de, prog, conv),
        u234_u238_cementum=_r0(c2, prog, conv),
        u234_u238_is="initial",
        radon_loss_enamel=rn_e, radon_loss_dentine=rn_d, radon_loss_cementum=rn_d,
        factors=factors_of(study),
        beta_by_segment=seg,
    )
    if ext_with_beta is not None:
        # the paper gives sediment beta + gamma + cosmic as one number: keep its
        # total, with EPRdating's sediment beta inside it
        sed_beta = {c.name: c.rate for c in tooth.components()}.get("sediment beta", 0.0)
        tooth = replace(tooth, gamma=max(ext_with_beta / 1000.0 - sed_beta, 0.0))
    us = None
    if _th(en) is not None and en.get("u234_u238") is not None:
        def d(t):
            return UseriesData(_v(_th(t)), _v(t["u234_u238"]))
        us = USESRSample(
            tooth, enamel=d(en),
            dentine=d(de) if inner_dentine and _th(de) is not None else None,
            cementum=d(c2) if outer_tissue and _th(c2) is not None else None,
        )
    if gamma is None and not gamma_from_sediment:
        notes.append("no gamma information")
    return Recalc(study["study_id"], s["sample"], prog, tooth, us, pub, external_from_paper, notes)


def components(r, T: float) -> dict[str, float]:
    """Mean dose rates (µGy/a) over T ka, grouped as papers print them."""
    acc = r.detail.accumulated if hasattr(r, "detail") else r.accumulated
    a = {k: v / T * 1000 for k, v in acc.items()}
    return {
        "internal": a.get("enamel alpha", 0) + a.get("enamel beta", 0),
        "beta_dentine": a.get("dentine beta", 0) + a.get("cementum beta", 0),
        "beta_sediment": a.get("sediment beta", 0),
        "gamma_cosmic": a.get("gamma", 0) + a.get("cosmic", 0),
        "beta_total": a.get("dentine beta", 0) + a.get("cementum beta", 0) + a.get("sediment beta", 0),
        "alpha_internal": a.get("enamel alpha", 0),
        "beta_all": a.get("enamel beta", 0) + a.get("dentine beta", 0) + a.get("cementum beta", 0)
        + a.get("sediment beta", 0),
        "beta_gamma_sediment_plus_cosmic": a.get("sediment beta", 0) + a.get("gamma", 0) + a.get("cosmic", 0),
        "total": sum(a.values()),
    }


def recompute(rc: Recalc) -> dict:
    """Ages (ka) and dose rates for every model the paper reports."""
    out = {}
    ages = rc.published.get("ages_ka") or {}
    for model, val in ages.items():
        if val is None or not isinstance(val, (list, tuple)) or val[0] is None:
            continue
        m = model.upper()
        if m in ("US", "CSUS") and rc.us is not None:
            res = rc.us.age(model=m)
            out[model] = {"age": res.age, "status": res.status,
                          "p_enamel": res.p_enamel, "p_dentine": res.p_dentine,
                          "dose_rates": components(res, res.age) if res.age else None}
        elif m in UPTAKE:
            t = rc.tooth
            res = replace(t, uptake_enamel=USModel(UPTAKE[m]), uptake_dentine=USModel(UPTAKE[m]),
                          uptake_cementum=USModel(UPTAKE[m])).age()
            out[model] = {"age": res.age, "status": "ok", "dose_rates": components(res, res.age)}
    return out


def isfinite(x) -> bool:
    return x is not None and math.isfinite(x)


def unconverged_dentine_p(s: dict, model: str = "US") -> bool:
    """True when the paper's own (age, dentine p) do not reproduce the measured
    dentine 230Th/234U (DATA's dentine p does not always converge, see
    DATA_FINDINGS.md)."""
    from eprdating.usesr import th230_u234

    pub = s.get("published") or {}
    age = (pub.get("ages_ka") or {}).get(model)
    pd, de = pub.get("p_dentine"), s.get("dentine") or {}
    if not age or not pd or pd[0] is None or _th(de) is None or de.get("u234_u238") is None:
        return False
    return abs(th230_u234(age[0], pd[0], _v(de["u234_u238"])) - _v(_th(de))) > 0.02


def benchmark() -> list[dict]:
    """One row per published age: study, sample, model, published and
    recomputed age, and why a row is set apart (``excluded``)."""
    rows = []
    for sid, st in load_studies().items():
        if not st.get("benchmark"):
            continue
        for s in st["samples"]:
            ages = (s.get("published") or {}).get("ages_ka") or {}
            if not any(isinstance(v, list) and v and v[0] is not None for v in ages.values()):
                continue
            base = {"study": sid, "sample": s["sample"], "program": program_of(st)}
            if s.get("exclude"):
                rows += [{**base, "model": m, "published": v[0], "ours": None, "excluded": s["exclude"]}
                         for m, v in ages.items() if isinstance(v, list) and v and v[0] is not None]
                continue
            try:
                rc = build(st, s)
                res = recompute(rc)
            except Exception as e:  # noqa: BLE001
                res, rc = {}, None
                err = repr(e)
            for m, v in ages.items():
                if not (isinstance(v, list) and v and v[0] is not None) or m.upper() not in MODELS:
                    continue
                r = res.get(m)
                row = {**base, "model": m, "published": v[0],
                       "ours": r["age"] if r else None, "status": r["status"] if r else "not computed",
                       "notes": rc.notes if rc else [err]}
                if r is None:
                    row["excluded"] = "model not reproducible from the published inputs"
                elif m.upper() == "US" and unconverged_dentine_p(s, m):
                    row["excluded"] = "published dentine p does not reproduce the measured dentine ratio"
                rows.append(row)
    return rows
