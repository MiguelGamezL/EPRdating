"""Freiberg MS5000 (ESRStudio) files, written here in the instrument's layout."""

import base64

import numpy as np
import pytest

from eprdating.spectra import (
    DEFAULT_WINDOW,
    IntensityWindow,
    combine_spectra,
    field_for_g,
    read_epr,
    read_ms5000,
)
from eprdating.spectra.freiberg import is_ms5000

F = 9.38625


def shape(B, amp):
    x = (B - field_for_g(2.0006, F)) / 0.3
    return amp * -2 * x / (1 + x**2) ** 2


def b64(v):
    # ESRStudio writes each double as its own padded Base64 chunk
    return "".join(base64.b64encode(np.float64(x).tobytes()).decode() for x in v)


def curves(amp, sweep_s=31.0, rng=None, offset=0.0):
    """Field and signal sampled in time, as the instrument does."""
    rng = rng or np.random.default_rng(0)
    tB = 0.7 + 0.02 * np.arange(int((sweep_s - 0.7) / 0.02))
    B = 331.0 + 12.0 * (tB - tB[0]) / (tB[-1] - tB[0])
    tY = 0.5855 + 0.026 * np.arange(int((sweep_s + offset - 0.5855) / 0.026))
    Bi = np.interp(tY, tB, B)
    y = shape(Bi, amp) + 5.0 * rng.standard_normal(tY.size)
    return tB, B, tY, y


def write_xml(path, amp, name=None, rng=None, offset=0.0):
    tB, B, tY, y = curves(amp, rng=rng, offset=offset)
    name = name or path.stem
    path.write_text(f"""<?xml version="1.0" encoding="utf-8"?>
<ESRXmlFile Version="1" Timestamp="2023-11-02T00:08:50Z">
  <Data>
    <Measurement Name="{name}" Device="MS-5000" MwFreq="{F}" GonAngle="180" Temperature="32.4" QFactor="1561.9">
      <Recipe Name="powder measurments 10accu" Type="single">
        <Parameters>
          <Param Name="Accumulations">10</Param>
          <Param Name="Bfrom">331</Param>
          <Param Name="Bto">343</Param>
          <Param Name="Modulation">0.1</Param>
          <Param Name="MicrowavePower">2</Param>
        </Parameters>
      </Recipe>
      <DataCurves>
        <Curve Name="BField" YType="BField" XOffset="0.7" XSlope="0.02" Compression="Base64">{b64(B)}</Curve>
        <Curve Name="BField_Temp" YType="BField_TempVoltage" XOffset="0" XSlope="0" Compression="Base64" />
        <Curve Name="MWAbsorption" YType="MW_Absorption" XOffset="0.5855" XSlope="0.026" Compression="Base64">{b64(y)}</Curve>
      </DataCurves>
    </Measurement>
  </Data>
</ESRXmlFile>""", encoding="utf-8")
    inside = (tY >= tB[0]) & (tY <= tB[-1])  # what the CSV export keeps
    return path, np.interp(tY, tB, B)[inside], y[inside]


def write_csv(path, B, y, name):
    keep = (B >= 331.0) & (B <= 343.0)
    rows = "\n".join(f"{b};{v + 235.8}" for b, v in zip(B[keep], y[keep], strict=True))
    path.write_text(f"﻿Name;{name}\n\nRecipe\nAccumulations;10;Number of accumulations\nBfrom;331;B from\n"
                    f"Bto;343;B to\nModulation;0.1;Modulation\nMicrowavePower;2;Microwave power\n\nMeas\n"
                    f"BField [mT];MW_Absorption []\n{rows}\n", encoding="utf-8")
    return path


def test_xml_and_csv(tmp_path):
    p, B, y = write_xml(tmp_path / "S_50Gy_1_result.xml", 4000.0)
    s = read_epr(p)
    assert is_ms5000(p) and s.name == "S_50Gy_1_result"
    assert (s.freq_GHz, s.power_mW, s.mod_amp_mT, s.gain) == (F, 2.0, 0.1, None)
    assert s.params["Accumulations"] == 10 and s.params["QFactor"] == pytest.approx(1561.9)
    assert s.B[0] >= 331.0 and s.B[-1] <= 343.0 and np.all(np.diff(s.B) > 0)
    assert np.ptp(s.y) == pytest.approx(np.ptp(shape(s.B, 4000.0)), rel=0.02)
    # the CSV export: same spectrum plus an offset; frequency from the .xml next to it
    c = read_epr(write_csv(tmp_path / "S_50Gy_1_result.csv", B, y, "S_50Gy_1_result"))
    assert is_ms5000(tmp_path / "S_50Gy_1_result.csv")
    assert c.freq_GHz == F and c.power_mW == 2.0
    assert np.allclose(c.y - c.y.mean(), np.interp(c.B, s.B, s.y) - np.interp(c.B, s.B, s.y).mean(), atol=1e-6)
    alone = tmp_path / "other"
    alone.mkdir()
    lone = read_epr(write_csv(alone / "x.csv", B, y, "x"))
    assert lone.freq_GHz is None and read_ms5000(alone / "x.csv", freq_GHz=9.4).freq_GHz == 9.4
    with pytest.raises(ValueError):
        (tmp_path / "bad.xml").write_text('<?xml version="1.0"?><ESRXmlFile><Data/></ESRXmlFile>')
        read_ms5000(tmp_path / "bad.xml")


def test_folder_of_runs_and_different_lengths(tmp_path):
    rng = np.random.default_rng(1)
    d = tmp_path / "S_100Gy_2"
    d.mkdir()
    for k in range(1, 11):  # runs of slightly different length, as the instrument gives
        write_xml(d / f"S_100Gy_2_{k}.xml", 5000.0, rng=rng, offset=0.03 * (k % 3))
    write_xml(d / "S_100Gy_2_result.xml", 5000.0, rng=rng)
    s = read_epr(d)
    assert s.n_scans == 10 and s.params["files"][0] == "S_100Gy_2_1.xml" and s.params["files"][-1] == "S_100Gy_2_10.xml"
    assert read_ms5000(d, include_result=True).n_scans == 11
    # repeats with a point more or less still combine (on the grid of the first)
    a, _, _ = write_xml(tmp_path / "a.xml", 5000.0, rng=rng)
    b, _, _ = write_xml(tmp_path / "b.xml", 5000.0, rng=rng, offset=0.05)
    sa, sb = read_epr(a), read_epr(b)
    assert sa.B.size != sb.B.size
    win = IntensityWindow(40, center_g=2.000)
    c = combine_spectra([sa, sb], win)
    assert c.spectrum.B.size == sa.B.size and len(c.members) == 2
    with pytest.warns(UserWarning, match="beyond the sweep"):
        DEFAULT_WINDOW.mask(sa.B, F)  # 100 G around 2.0023 does not fit a 12 mT sweep at this frequency


def test_gui_reads_ms5000_files(tmp_path):
    pytest.importorskip("ipywidgets")
    import matplotlib

    matplotlib.use("Agg")
    from eprdating.gui import app

    rng = np.random.default_rng(2)
    for d in (0, 200, 400, 800):
        for rot in (1, 2):
            name = f"S_N_{rot}_result.xml" if d == 0 else f"S_{d}Gy_{rot}_result.xml"
            write_xml(tmp_path / name, 10.0 * (d + 300) * (1 + 0.01 * rng.standard_normal()), rng=rng)
    de = app().de
    de.folder.value = str(tmp_path)
    de.load_folder.click()
    assert len(de.spectra) == 8
    for n in de.rows:
        if "_N_" in n:
            de.set_file(n, dose=0)
    de.group_by_dose()
    de.width.value, de.center_g.value = 40.0, 2.000
    de.method.value = "peak_to_peak"
    res = de.compute()
    assert len(res) == 4 and all(r["n_repeats"] == 2 for r in res)
    assert abs(de.drc.De - 300) < 3 * de.drc.De_sigma + 30
