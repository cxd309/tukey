"""Conditioning and reductions: detrending, interpolation, filling gaps, RMS"""

from typing import Final

import numpy as np
import scipy as sc

from dsp.cases import SIGNAL_LENGTH
from dsp.signals import synthetic_signals
from write import write_vector

NAN: Final = float("nan")


def generate_detrend() -> None:
    """scipy.signal.detrend, linear and constant, on signals with and without a trend"""
    signals = synthetic_signals(SIGNAL_LENGTH)
    rng = np.random.default_rng(seed=1)
    signals["ramp_noise"] = (
        0.02 * np.arange(SIGNAL_LENGTH) + 3 + rng.standard_normal(SIGNAL_LENGTH)
    )
    signals["short2"] = np.array(
        [1.0, 4.0]
    )  # two points: the line through them is exact
    for kind in ("linear", "constant"):
        for name, x in signals.items():
            write_vector(
                "detrend",
                f"{name}_{kind}",
                f"scipy.signal.detrend(x, type={kind!r})",
                {"x": x.tolist(), "type": kind},
                {"y": sc.signal.detrend(x, type=kind).tolist()},
            )


def generate_interp1() -> None:
    """linear interpolation, NaN outside the data range: scipy.interpolate.interp1d
    with bounds_error=False, which is also MATLAB's interp1"""
    rng = np.random.default_rng(seed=2)
    nonuniform = np.sort(rng.uniform(0, 10, 20))
    cases = {
        # query points below, inside, on and above the knots
        "uniform": (np.arange(10.0), np.sin(np.arange(10.0)), np.linspace(-1, 10, 45)),
        "nonuniform": (
            nonuniform,
            rng.standard_normal(20),
            np.sort(rng.uniform(-1, 11, 50)),
        ),
        "on_knots": (np.arange(10.0), np.arange(10.0) ** 2, np.arange(10.0)),
        # resampling a series onto a finer grid, as the app does
        "resample": (
            np.arange(float(SIGNAL_LENGTH)),
            synthetic_signals(SIGNAL_LENGTH)["chirp"],
            np.arange(0, SIGNAL_LENGTH - 1 + 1e-9, 0.25),
        ),
    }
    for name, (x, y, xi) in cases.items():
        yi = sc.interpolate.interp1d(x, y, bounds_error=False)(xi)
        write_vector(
            "interp1",
            name,
            "scipy.interpolate.interp1d(x, y, bounds_error=False)(xi)",
            {"x": x.tolist(), "y": y.tolist(), "xi": xi.tolist()},
            {"yi": yi.tolist()},
        )


def fill_linear(x: np.ndarray) -> np.ndarray:
    """MATLAB's fillmissing(x, 'linear'), from its documentation, as SciPy has none:
    gaps interpolated linearly and the ends extrapolated linearly from the known
    values, with fewer than two known values, x is returned unchanged.
    testdata/matlab/fillmissing checks this definition"""
    known = ~np.isnan(x)
    if known.sum() < 2:
        return x.copy()
    i = np.arange(len(x), dtype=float)
    return sc.interpolate.interp1d(i[known], x[known], fill_value="extrapolate")(i)


def generate_fillmissing() -> None:
    """MATLAB's fillmissing, linear and constant, from its documented definition"""
    cases = {
        "none_missing": [1.0, 2.0, 4.0, 8.0],
        "one_gap": [1.0, 2.0, NAN, 4.0, 5.0],
        "long_gap": [0.0, NAN, NAN, NAN, 4.0, 3.0],
        "several_gaps": [1.0, NAN, 3.0, NAN, NAN, 0.0, NAN, 2.0],
        "leading": [NAN, NAN, 2.0, 3.0, 5.0],
        "trailing": [1.0, 2.0, 4.0, NAN, NAN],
        "both_ends": [NAN, 1.0, NAN, 3.0, NAN],
        "one_known": [NAN, 2.0, NAN],
        "all_missing": [NAN, NAN, NAN],
        "empty": [],
    }
    for name, values in cases.items():
        x = np.array(values)
        write_vector(
            "fillmissing",
            f"{name}_linear",
            "fillmissing(x, 'linear') from MATLAB's documentation",
            {"x": x.tolist(), "method": "linear", "value": None},
            {"y": fill_linear(x).tolist()},
        )
        write_vector(
            "fillmissing",
            f"{name}_constant",
            "fillmissing(x, 'constant', 0) from MATLAB's documentation",
            {"x": x.tolist(), "method": "constant", "value": 0.0},
            {"y": np.where(np.isnan(x), 0.0, x).tolist()},
        )


def generate_rms() -> None:
    """root mean square, sqrt(mean(x^2)), on the synthetic signals; SciPy has none,
    so NumPy's definition, which MATLAB's rms is too"""
    for name, x in synthetic_signals(SIGNAL_LENGTH).items():
        write_vector(
            "rms",
            name,
            "numpy.sqrt(numpy.mean(x**2))",
            {"x": x.tolist()},
            {"rms": [float(np.sqrt(np.mean(x**2)))]},
        )
