import struct

import numpy as np
import pytest
from test_gamma import ROOM, TRUE, synth

from eprdating.gamma import (
    IAEA_RGK_1,
    IAEA_RGTH_1,
    IAEA_RGU_1,
    Calibration,
    analyse_files,
    auto_calibrate,
    read_chn,
    read_gamma,
    read_n42,
    read_spe,
)
from eprdating.gamma.spectrum import GammaSpectrum


def write_spe(path, spec, coef=(0.5, 0.35)):
    n = spec.counts.size
    txt = ["$SPEC_ID:", "test", "$DATE_MEA:", "08/01/2024 12:28:55", "$MEAS_TIM:",
           f"{spec.live_time:.0f} {spec.real_time or spec.live_time:.0f}", "$DATA:", f"0 {n - 1}"]
    txt += [f"{int(c):8d}" for c in spec.counts]
    txt += ["$ENER_FIT:", f"{coef[0]} {coef[1]}", "$MCA_CAL:", "3", f"{coef[0]} {coef[1]} 0.0 keV"]
    path.write_text("\n".join(txt) + "\n")
    return path


def write_chn(path, spec, coef=(0.5, 0.35)):
    n = spec.counts.size
    head = struct.pack("<hhh2sii8s4shh", -1, 1, 1, b"55", round((spec.real_time or spec.live_time) / 0.02),
                       round(spec.live_time / 0.02), b"01AUG241", b"1228", 0, n)
    body = spec.counts.astype("<u4").tobytes()
    trailer = struct.pack("<hh3f", -102, 0, coef[0], coef[1], 0.0) + b"\0" * 496
    path.write_bytes(head + body + trailer)
    return path


def write_n42(path, spec, coef=(0.5, 0.35), compress=True):
    vals = []
    zeros = 0
    for c in spec.counts.astype(int):
        if compress and c == 0:
            zeros += 1
            continue
        if zeros:
            vals += ["0", str(zeros)]
            zeros = 0
        vals.append(str(c))
    if zeros:
        vals += ["0", str(zeros)]
    code = ' compressionCode="CountedZeroes"' if compress else ""
    path.write_text(f"""<?xml version="1.0"?>
<RadInstrumentData xmlns="http://physics.nist.gov/N42/2011/N42">
 <EnergyCalibration id="ec"><CoefficientValues>{coef[0]} {coef[1]} 0</CoefficientValues></EnergyCalibration>
 <RadMeasurement id="m1">
  <StartDateTime>2024-08-01T12:28:55Z</StartDateTime>
  <RealTimeDuration>PT{spec.real_time or spec.live_time}S</RealTimeDuration>
  <Spectrum id="s1" energyCalibrationReference="ec">
   <LiveTimeDuration>PT{spec.live_time}S</LiveTimeDuration>
   <ChannelData{code}>{' '.join(vals)}</ChannelData>
  </Spectrum>
 </RadMeasurement>
</RadInstrumentData>
""")
    return path


@pytest.fixture
def spec():
    s = synth({"K": 1.0, "Ra226": 0.5, "Th232": 1.0}, seed=7)
    s.real_time = 86500.0
    return s


@pytest.mark.parametrize("writer,ext", [(write_spe, ".Spe"), (write_chn, ".Chn"), (write_n42, ".n42")])
def test_readers_round_trip(tmp_path, spec, writer, ext):
    p = writer(tmp_path / f"x{ext}", spec)
    r = read_gamma(p)
    np.testing.assert_array_equal(r.counts, spec.counts)
    assert r.live_time == pytest.approx(spec.live_time)
    assert r.real_time == pytest.approx(spec.real_time)
    assert r.calibration.coef[:2] == pytest.approx((0.5, 0.35))
    assert r.channels[0] == 0


def test_specific_readers_and_details(tmp_path, spec):
    assert read_spe(write_spe(tmp_path / "a.spe", spec)).start == "08/01/2024 12:28:55"
    assert read_chn(write_chn(tmp_path / "a.chn", spec)).start == "2024-08-01 12:28:55"
    r = read_n42(write_n42(tmp_path / "a.xml", spec, compress=False))
    np.testing.assert_array_equal(r.counts, spec.counts)


def test_column_files_need_live_time(tmp_path, spec):
    p = tmp_path / "c.csv"
    np.savetxt(p, np.c_[np.arange(spec.counts.size), spec.counts], delimiter=",", fmt="%d")
    with pytest.raises(ValueError, match="live_time"):
        read_gamma(p)
    r = read_gamma(p, live_time=1000.0)
    np.testing.assert_array_equal(r.counts, spec.counts)
    assert r.live_time == 1000.0


def test_auto_calibration_without_guess():
    cal = Calibration((-3.0, 0.4012), (0.36, 0.0004))
    s = synth({"K": 1.0, "Ra226": 0.5, "Th232": 1.0}, cal=cal, seed=8)
    s.calibration = None
    got = auto_calibrate(s)
    for e in (300.0, 1460.82, 2614.5):
        assert got.energy(cal.channel(e)) == pytest.approx(e, abs=0.15)
    assert got.sigma_keV(1332) == pytest.approx(cal.sigma_keV(1332), rel=0.15)


def test_scintillator_resolution_is_rejected():
    cal = Calibration((0.0, 1.0), (0.0, 1.5))  # FWHM ≈ 74 keV at 662 keV
    s = synth({"K": 3.0, "Ra226": 2.0, "Th232": 2.0}, cal=cal, seed=9)
    with pytest.raises(RuntimeError):
        auto_calibrate(s)


def test_analyse_files_end_to_end(tmp_path):
    def save(name, act, seed, cal):
        sp = synth(act, seed=seed, cal=cal, name=name)
        return write_spe(tmp_path / f"{name}.Spe", sp, coef=(0.0, 0.25))  # stored calibration is wrong

    U = save("U", {"Ra226": 4.0, "U238": 4.0}, 1, Calibration((2.3, 0.35173), (0.35, 0.0004)))
    Th = save("Th", {"Th232": 8.0}, 2, Calibration((2.4, 0.35177), (0.36, 0.0004)))
    K = save("K", {"K": 0.448}, 3, Calibration((3.0, 0.35156), (0.37, 0.0004)))
    bg = save("bg", {}, 4, Calibration((2.4, 0.35154), (0.31, 0.0004)))
    soil = save("soil", {"Ra226": 0.02, "U238": 0.02, "Th232": 0.08, "K": 0.03}, 5, TRUE)
    r = analyse_files(soil, 500.0, {"U": (U, IAEA_RGU_1), "Th": (Th, IAEA_RGTH_1), "K": (K, IAEA_RGK_1)},
                      background=bg)
    for got, true in ((r.K, 3.0), (r.U, 2.0), (r.Th, 8.0)):
        assert abs(got.value - true) < 3 * got.sigma
    assert ROOM["K"] > 0  # background lines exist in every synthetic spectrum


def test_gamma_spectrum_from_arrays():
    s = GammaSpectrum(np.ones(10), 100.0)
    assert s.channels[0] == 1 and s.dead_time_fraction is None
