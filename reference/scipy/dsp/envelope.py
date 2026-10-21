"""The analytic signal and envelopes: hilbert, envelope, and MATLAB's RMS envelope"""

from typing import Final

import numpy as np
import scipy as sc

from dsp.cases import SIGNAL_LENGTH
from dsp.signals import synthetic_signals
from write import complex_json, write_vector

# lengths for hilbert: odd and even differ at the Nyquist bin
HILBERT_LENGTHS: Final = (1, 2, 3, 8, 9, 64, 255)

# RMS envelope window lengths: odd, even, and longer than the 200-sample signals
RMS_WINDOWS: Final = (1, 5, 10, 50, 250)


def generate_hilbert() -> None:
    """scipy.signal.hilbert, the analytic signal, at odd and even lengths"""
    for n in HILBERT_LENGTHS:
        x = np.random.default_rng(seed=n + 100).standard_normal(n)
        write_vector(
            "hilbert",
            f"n{n}",
            "scipy.signal.hilbert(x)",
            {"x": x.tolist()},
            {"analytic": complex_json(sc.signal.hilbert(x))},
        )


def generate_envelope() -> None:
    """upper and lower envelopes, analytic and RMS, on the synthetic signals"""
    for name, x in synthetic_signals(SIGNAL_LENGTH).items():
        # MATLAB's envelope(x): mean +- |hilbert(x - mean)|, which keeps the
        # Nyquist bin, unlike scipy.signal.envelope's default band
        mean = x.mean()
        env = np.abs(sc.signal.hilbert(x - mean))
        write_vector(
            "envelope",
            f"{name}_analytic",
            "mean +- abs(scipy.signal.hilbert(x - mean))",
            {"x": x.tolist(), "mode": "analytic", "n": None},
            {"upper": (mean + env).tolist(), "lower": (mean - env).tolist()},
        )
        for n in RMS_WINDOWS:
            upper, lower = rms_envelope(x, n)
            write_vector(
                "envelope",
                f"{name}_rms{n}",
                f"envelope(x, {n}, 'rms') from MATLAB's documentation",
                {"x": x.tolist(), "mode": "rms", "n": n},
                {"upper": upper.tolist(), "lower": lower.tolist()},
            )


def rms_envelope(x: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    """
    MATLAB's envelope(x, n, 'rms'), as SciPy has none: the mean plus and minus
    the moving RMS of x less its mean, over a window of n samples. For even n,
    the extra sample is after the current one, so n = 10 covers k-4 to k+5,
    the opposite of movmean's convention. The window shrinks at the ends.
    testdata/matlab/envelope checks this definition
    """
    mean = x.mean()
    squared = (x - mean) ** 2
    before, after = (n - 1) // 2, n // 2
    rms = np.array(
        [
            np.sqrt(squared[max(0, k - before) : k + after + 1].mean())
            for k in range(len(x))
        ]
    )
    return mean + rms, mean - rms
