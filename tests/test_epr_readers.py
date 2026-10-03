import numpy as np
import pytest

from eprdating.spectra import read_epr


def _deriv(n=512):
    x = np.linspace(-5, 5, n)
    return -x / (1 + x**2) ** 2


def write_bes3t(tmp_path, y, order="BIG", fmt="D", cplx=False, unit="G"):
    dsc = f"""#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*	Dataset Type and Format:
DSRC	EXP
BSEQ	{order}
IKKF	{"CPLX" if cplx else "REAL"}
XTYP	IDX
YTYP	NODATA
IRFMT	{fmt}
XPTS	{y.size}
XMIN	3300.000000
XWID	100.000000
TITL	'test sample'
XUNI	'{unit}'
#SPL	1.2 * STANDARD PARAMETER LAYER
EXPT    CW
MWFQ    9.43e+09
MWPW    0.019
RCAG    40
B0MA    0.0004
RCTC    0.04096
AVGS    4
#DSL	1.0 * DEVICE SPECIFIC LAYER
.DVC     acqStart, 1.0
"""
    (tmp_path / "s.DSC").write_text(dsc)
    dt = (">" if order == "BIG" else "<") + {"D": "f8", "F": "f4", "I": "i4"}[fmt]
    data = np.c_[y, np.zeros_like(y)].ravel() if cplx else y
    data.astype(dt).tofile(tmp_path / "s.DTA")
    return tmp_path / "s.DSC"


def write_esp(tmp_path, y, winepr=True):
    lines = (["DOS  Format"] if winepr else []) + [
        "JSS 0", f"RES {y.size}", "HCF 3350.0", "HSW 100.0", "MF  9.43", "MP  19.0",
        "RRG 1.0e+04", "RMA 4.0", "RTC 40.96", "JSD 4", "JEX field-sweep"]
    (tmp_path / "e.par").write_text("\r\n".join(lines))
    (y.astype("<f4") if winepr else np.round(y * 1e6).astype(">i4")).tofile(tmp_path / "e.spc")
    return tmp_path / "e.par"


@pytest.mark.parametrize("order,fmt,cplx", [("BIG", "D", False), ("LIT", "F", False), ("BIG", "I", True)])
def test_bes3t(tmp_path, order, fmt, cplx):
    y = np.round(_deriv() * 1000)
    s = read_epr(write_bes3t(tmp_path, y, order, fmt, cplx))
    np.testing.assert_allclose(s.y, y)
    assert s.B[0] == pytest.approx(330.0) and s.B[-1] == pytest.approx(340.0)
    assert s.freq_GHz == pytest.approx(9.43) and s.power_mW == pytest.approx(19.0)
    assert s.gain == pytest.approx(100.0) and s.mod_amp_mT == pytest.approx(0.4)
    assert s.time_constant_ms == pytest.approx(40.96)
    assert s.params["TITL"] == "test sample"


def test_bes3t_rejects_non_field_axis(tmp_path):
    with pytest.raises(ValueError, match="magnetic field"):
        read_epr(write_bes3t(tmp_path, _deriv(), unit="ns"))


@pytest.mark.parametrize("winepr", [True, False])
def test_esp_winepr(tmp_path, winepr):
    y = _deriv()
    s = read_epr(write_esp(tmp_path, y, winepr))
    scale = 1.0 if winepr else 1e6
    np.testing.assert_allclose(s.y / scale, y, atol=1e-6)
    assert s.B[0] == pytest.approx(330.0) and s.B[-1] == pytest.approx(340.0)
    assert (s.freq_GHz, s.power_mW, s.gain, s.mod_amp_mT) == pytest.approx((9.43, 19.0, 1e4, 0.4))
    assert read_epr(tmp_path / "e.spc").B.size == y.size  # either file of the pair


def test_columns(tmp_path):
    y = _deriv(100)
    B = np.linspace(3300, 3400, 100)
    np.savetxt(tmp_path / "c.csv", np.c_[B[::-1], y[::-1], 2 * y[::-1]], delimiter=",", header="B,s1,s2")
    s = read_epr(tmp_path / "c.csv", freq_GHz=9.43)
    assert s.n_scans == 2 and s.B[0] == pytest.approx(330.0)  # sorted, G -> mT
    np.testing.assert_allclose(s.scans[0], y)
