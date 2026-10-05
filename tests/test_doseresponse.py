import numpy as np
import pytest

from eprdating import bootstrap_De, fit_dose_response


def synthetic(model_De=150.0, Imax=1000.0, D0=2000.0, noise=0.01, seed=1):
    rng = np.random.default_rng(seed)
    D = np.array([0, 50, 100, 200, 400, 800, 1600, 3200, 6400, 10000], float)
    I = Imax * (1 - np.exp(-(D + model_De) / D0))
    return D, I * (1 + noise * rng.standard_normal(D.size))


def test_sse_recovers_De():
    D, I = synthetic(noise=0.0)
    r = fit_dose_response(D, I, "SSE")
    assert r.De == pytest.approx(150.0, rel=1e-4)
    assert r.params["D0"] == pytest.approx(2000.0, rel=1e-4)


def test_noisy_fit_within_uncertainty():
    D, I = synthetic(noise=0.01, seed=7)
    r = fit_dose_response(D, I, "SSE", weighting="1/I^2")
    assert abs(r.De - 150.0) < 4 * r.De_sigma
    assert r.r2 > 0.999


@pytest.mark.parametrize("model", ["EXPLIN", "DSE", "LIN"])
def test_other_models_run(model):
    D, I = synthetic(noise=0.005)
    if model == "LIN":
        D, I = D[:5], I[:5]
    r = fit_dose_response(D, I, model)
    assert 100 < r.De < 200


def test_max_dose_cut_and_errors():
    D, I = synthetic()
    r = fit_dose_response(D, I, "SSE", max_dose=3200)
    assert r.dose.max() == 3200
    with pytest.raises(ValueError):
        fit_dose_response(D[:3], I[:3], "SSE")
    with pytest.raises(ValueError):
        fit_dose_response(D, I, "NOPE")


def test_samples_and_bootstrap():
    D, I = synthetic(noise=0.01, seed=3)
    r = fit_dose_response(D, I, "SSE")
    s = r.De_samples(2000, np.random.default_rng(0))
    assert np.std(s) == pytest.approx(r.De_sigma, rel=0.1)
    b = bootstrap_De(r, n=100, seed=0)
    assert b.size > 90 and abs(np.median(b) - r.De) < 3 * r.De_sigma


def test_linear_fit_with_a_noisy_low_dose_end():
    # the first four points alone fall with dose; the fit must still start inside its bounds
    D = np.arange(0, 181, 20.0)
    I = np.array([0.9, 0.3, 0.2, 0.5, 1.37, 1.2, 1.38, 2.25, 1.99, 1.80])
    f = fit_dose_response(D, I, "LIN", De_min=-np.inf)
    assert np.isfinite(f.De)
    assert list(f.params.values())[1] > 0
