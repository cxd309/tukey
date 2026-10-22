"""Finding a template in a signal: MATLAB's findsignal, with its default options"""

import math

import numpy as np

from dsp.cases import SIGNAL_LENGTH
from dsp.signals import synthetic_signals
from write import write_vector


def find_signal(data: np.ndarray, signal: np.ndarray) -> tuple[int, np.ndarray]:
    """
    MATLAB's findsignal(data, signal), as SciPy has none: the window of data
    closest to signal by the sum of squared differences, with no normalisation
    and fixed alignment, MATLAB's documented defaults. Ties go to the earliest
    window. Each window is summed exactly with math.fsum, so the distances are
    the true ones rather than one implementation's rounding.
    testdata/matlab/findsignal checks this definition
    """
    m = len(signal)
    profile = np.array(
        [
            math.fsum(((data[i : i + m] - signal) ** 2).tolist())
            for i in range(len(data) - m + 1)
        ]
    )
    return int(np.argmin(profile)), profile


def findsignal_cases() -> list[tuple[str, np.ndarray, np.ndarray]]:
    """(name, data, signal) for each vector"""
    signals = synthetic_signals(SIGNAL_LENGTH)
    chirp, noise, sines = signals["chirp"], signals["noise"], signals["sines"]
    rng = np.random.default_rng(seed=7)
    # as the app: a slowly wandering reference, and the same signal delayed and
    # noisier, searched for a 1000-sample feature cut from the reference
    walk = np.cumsum(rng.standard_normal(11000))
    delayed = walk[500:10500] + 0.1 * rng.standard_normal(10000)
    return [
        # the template cut from the data, as the app's reference calls: distance exactly 0
        ("self_chirp", chirp, chirp[60:100]),
        ("self_noise", noise, noise[120:150]),
        ("self_walk", walk[:10000], walk[4000:5000]),
        ("at_start", chirp, chirp[:30]),
        ("at_end", chirp, chirp[-30:]),
        ("whole", sines, sines),  # one window only
        ("one_sample", noise, np.array([0.5])),
        # a template that's no exact match
        ("noisy_copy", noise, noise[50:90] + 0.05 * rng.standard_normal(40)),
        (
            "offset",
            chirp,
            chirp[100:140] + 3.0,
        ),  # no normalisation, so the offset counts
        # exact ties: the template matches every period, so which start wins
        ("ties", np.tile([0.0, 1, 2, 1], 10), np.array([0.0, 1, 2, 1])),
        # as the app's second call: the feature found in a delayed, noisier copy
        ("delayed_walk", delayed, walk[4000:5000]),
    ]


def generate_findsignal() -> None:
    """findsignal with its default options, on each case"""
    for name, data, signal in findsignal_cases():
        start, profile = find_signal(data, signal)
        write_vector(
            "findsignal",
            name,
            "findsignal(data, signal), summed exactly with math.fsum",
            {"data": data.tolist(), "signal": signal.tolist()},
            {
                "start": start,
                "stop": start + len(signal) - 1,
                "dist": [profile[start]],
                "profile": profile.tolist(),
            },
        )
