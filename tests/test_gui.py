"""The interface, driven from code with synthetic files."""

import json

import numpy as np
import pytest

pytest.importorskip("ipywidgets")
import matplotlib

matplotlib.use("Agg")

from eprdating.gui import app
from eprdating.gui.de import dose_from_name
from eprdating.gui.gamma import guess_role
from eprdating.spectra import field_for_g

FREQ = 9.43
DE_TRUE = 40.0


def write_dat(path, amp, rng, nscans=4, shift_mT=0.0, n=512, sweep_G=500.0):
    """A .dat/.par pair in the UNAL format (field in G, by default 512 points over 50 mT)."""
    BG = np.linspace(3355.0 - sweep_G / 2, 3355.0 + sweep_G / 2, n)
    B = BG / 10
    c = field_for_g(2.0006, FREQ) + shift_mT
    x = (B - c) / (0.5 * np.sqrt(3) / 2)
    shape = -2 * x / (1 + x**2) ** 2
    shape /= np.ptp(shape)
    lines, k = [], 0
    for _ in range(nscans):
        y = 3.6 + amp * shape + 0.05 * rng.standard_normal(B.size)
        for i, (b, v) in enumerate(zip(BG, y, strict=True)):
            lines.append(f"{k} {i} {b:.6f} {v * 1e-6:.6E} 1")
            k += 1
    path.write_text("\n".join(lines))
    path.with_suffix(".par").write_text(
        f"N : {n}\nCF : 3355.0\nCF_ : 3355.0\nSW : {sweep_G}\nNscans : {nscans}\nFreq : {FREQ}\nRG : 30\n")
    return path


@pytest.fixture
def series(tmp_path):
    rng = np.random.default_rng(3)
    paths = []
    for d in (0, 50, 100, 150, 200, 250):
        amp = 0.004 * (d + DE_TRUE)
        paths.append(write_dat(tmp_path / f"S_{d}Gy_19mW.dat", amp, rng, shift_mT=rng.normal(0, 0.05)))
    paths.append(write_dat(tmp_path / "S_100Gy_19mW_repeat.dat", 0.004 * (100 + DE_TRUE), rng, nscans=1))
    return tmp_path, paths


def test_names():
    assert dose_from_name("M18_200Gy_19mW.dat") == "200"
    assert dose_from_name("tooth_12,5 Gy.DSC") == "12.5"
    assert dose_from_name("M18_3_19mW_4SCAN.dat") == ""
    assert guess_role("fondo_24h.txt") == "background"
    assert guess_role("RGU1_24h.Spe") == "U"
    assert guess_role("RGTh1.Spe") == "Th"
    assert guess_role("K_24h.txt") == "K"
    assert guess_role("corte0_24h.txt") is None


def test_de_tab_from_files_to_age_tab(series):
    folder, _paths = series
    a = app()
    de = a.de
    de.folder.value = str(folder)
    de.load_folder.click()  # .par companions are skipped
    assert len(de.spectra) == 7 and "S_0Gy_19mW.dat" in de.rows
    de.group_by_dose()  # the 100 Gy repeat joins its aliquot
    res = de.compute()
    assert len(res) == 6
    rep = next(r for r in res if r["dose_Gy"] == 100)
    assert rep["n_repeats"] == 2 and rep["chi2_red"] is not None
    assert all(r["p_noise"] is not None for r in res)
    d = de.drc
    assert d.model == "LIN" and abs(d.De - DE_TRUE) < 3 * d.De_sigma
    assert a.age.De.value_widget.value == pytest.approx(d.De, rel=1e-3)  # flowed into the Age tab
    assert "De =" in de.fit_text.value

    de.set_point("150 Gy", False)
    de.fit()
    assert de.drc.dose.size == 5
    de.model.value = "SSE"
    de.fit_btn.click()
    assert de.drc.model == "SSE"

    s = de.settings()
    json.dumps(s)
    assert s["window"] == {"width": 100.0, "unit": "G", "center": 2.0023, "center_as": "g"}
    # the centre can be given as a field, in G or mT, and keeps its place when the unit changes
    de.center_mode.value = "G"
    assert de.center_g.value == pytest.approx(3364.9, abs=0.1)
    assert de.window().center(FREQ) == pytest.approx(336.49, abs=0.01)
    de.center_mode.value = "mT"
    assert de.center_g.value == pytest.approx(336.49, abs=0.01)
    de.center_g.value = 336.0
    de.center_mode.value = "g"
    assert de.center_g.value == pytest.approx(2.0052, abs=1e-4)
    de.center_mode.value = "G"
    de.center_g.value = 3360.0
    assert de.window().bounds(None) == pytest.approx((331.0, 341.0))
    de.compute_btn.click()
    assert "331.00-341.00 mT" in de.status.value
    assert s["fit"]["points"]["150 Gy"] is False
    assert de.results_csv().startswith("aliquot,dose_Gy")
    de.export_btn.click()
    assert "download=" in de.export_html.value


def test_compare_methods(series):
    _folder, paths = series
    de = app().de
    de.load(paths)
    de.group_by_dose()
    de.compute()
    de.set_point("250 Gy", False)
    de.compare_btn.click()
    rows = {r["method"]: r for r in de.comparison}
    assert set(rows) == {"template", "peak_to_peak", "t1_b2", "double_integral"}
    for r in rows.values():
        assert r["points"] == 5 and not r["problem"]
        assert abs(r["De_Gy"] - DE_TRUE) < 3 * r["sigma_Gy"] + 10, r
    assert "within the errors" in de.compare_html.value
    de.export_btn.click()
    assert "eprdating_De_methods.csv" in de.export_html.value
    assert de.comparison_csv().startswith("method,De_Gy")
    de.compute()  # new intensities: the comparison is cleared
    assert de.comparison == [] and de.compare_html.value == ""


def test_de_tab_reports_problems_in_the_panel(series):
    _folder, paths = series
    de = app().de
    de.load(paths)
    de.set_file("S_0Gy_19mW.dat", aliquot="x")
    de.rows["S_0Gy_19mW.dat"]["dose"].value = ""
    de.compute_btn.click()  # no traceback: the message is shown
    assert "no dose" in de.status.value
    de.rows["S_0Gy_19mW.dat"]["dose"].value = "0"
    for name in de.rows:
        de.set_file(name, use=name == "S_0Gy_19mW.dat")
    de.compute_btn.click()
    assert "at least two aliquots" in de.status.value
    de.fit_btn.click()
    assert "compute the intensities first" in de.status.value


def test_gamma_tab_to_age_tab(tmp_path):
    from test_gamma import _write, synth

    from eprdating.gamma import Calibration

    cal = Calibration((2.2, 0.35175), (0.33, 0.0004))
    specs = {
        "RGU1_24h.txt": synth({"Ra226": 4.0, "U238": 4.0}, seed=1, cal=cal),
        "RGTh1_24h.txt": synth({"Th232": 8.0}, seed=2, cal=cal),
        "RGK1_24h.txt": synth({"K": 0.448}, seed=3, cal=cal),
        "fondo_24h.txt": synth({}, seed=4, cal=cal, continuum=0.5),
        "soil_24h.txt": synth({"Ra226": 0.024, "U238": 0.024, "Th232": 0.096, "K": 0.036}, seed=5, cal=cal,
                              live=86400 * 4),
    }
    paths = []
    for name, sp in specs.items():
        p = _write(tmp_path, sp, cal)
        paths.append(p.rename(tmp_path / name))
    a = app()
    g = a.gamma
    g.load(paths)
    assert (g.sample.value, g.background.value) == ("soil_24h.txt", "fondo_24h.txt")
    assert (g.ref["U"].value, g.ref["Th"].value, g.ref["K"].value) == ("RGU1_24h.txt", "RGTh1_24h.txt",
                                                                      "RGK1_24h.txt")
    g.sample_mass.value = 600.0
    g.water.value_widget.value = 0.12
    g.analyse_btn.click()
    assert g.result is not None, g.status.value
    for el, true in (("U", 2.0), ("Th", 8.0), ("K", 3.0)):
        got = getattr(g.result, el)
        assert abs(got.value - true) < 3 * got.sigma + 0.05 * true
    assert a.age.sed_U.value_widget.value == pytest.approx(g.result.U.value, rel=1e-3)
    assert a.age.sed_water.value_widget.value == pytest.approx(0.12)
    json.dumps(g.settings())


def test_age_tab():
    a = app()
    age = a.age
    age.set_De(60.0, 6.0)
    age.n_mc.value = 50
    age.run.click()
    assert age.result is not None, age.status.value
    assert 10 < age.result.age < 200 and age.mc.samples.size > 40
    age.up_e.value = "US (give p)"
    age.p_e.value = 1.5
    age.cosmic_mode.value = "given"
    age.gamma_mode.value = "measured in situ"
    age.gamma.value_widget.value = 0.6
    age.alpha_escape.value = True
    age.compute()
    assert "dentine alpha" in age.result.components
    assert age.settings()["alpha_escape"] is True
    json.dumps(age.settings())
    age.export_btn.click()
    assert "download=" in age.export_html.value
    age.set_De(0.0, 0.0)
    age.run.click()
    assert "positive De" in age.status.value


def test_upload_values_of_ipywidgets_7_and_8(series):
    from types import SimpleNamespace

    from eprdating.gui._common import uploaded

    v8 = SimpleNamespace(value=({"name": "a.dat", "content": memoryview(b"xy")},))
    v7 = SimpleNamespace(value={"a.dat": {"metadata": {}, "content": b"xy"}})
    assert uploaded(v8) == uploaded(v7) == [("a.dat", b"xy")]
    # an upload goes through the pool, so a .dat finds its .par
    _folder, paths = series
    de = app().de
    for p in paths:
        de.pool.add_bytes(p.name, p.read_bytes())
        de.pool.add_bytes(p.with_suffix(".par").name, p.with_suffix(".par").read_bytes())
    de.load(list(de.pool.paths.values()))
    assert len(de.spectra) == 7 and all(s.freq_GHz == 9.43 for s in de.spectra.values())


def test_wide_sweep_natural_is_bridged(tmp_path):
    rng = np.random.default_rng(8)
    paths = [write_dat(tmp_path / f"S_{d}Gy.dat", 0.004 * (d + DE_TRUE), rng) for d in (50, 100, 150, 200, 250)]
    wide = {"n": 4096, "sweep_G": 5000.0}
    paths.append(write_dat(tmp_path / "S_0Gy_wide.dat", 0.004 * DE_TRUE * 1.3, rng, nscans=8, **wide))
    paths.append(write_dat(tmp_path / "S_250Gy_wide.dat", 0.004 * (250 + DE_TRUE) * 1.3, rng, nscans=8, **wide))
    de = app().de
    de.load(paths)
    de.group_by_dose()
    res = {r["aliquot"]: r for r in de.compute()}
    assert res["0 Gy"]["scale"].startswith("×")  # the wide-sweep natural is rescaled
    assert "factor" in de.status.value
    assert res["250 Gy"]["n_repeats"] == 1  # the bridge aliquot keeps its main-sweep measurement
    assert abs(de.drc.De - DE_TRUE) < 3 * de.drc.De_sigma
