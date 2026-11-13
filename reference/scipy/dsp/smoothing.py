"""
Smoothing: the median filter

The median filter is computed from its definition rather than taken from
scipy.signal.medfilt, which only zero-pads and whose NaN handling is the sort
order's accident. Each sample is the median of the n samples centred on it:

- zeropad: the window runs past the ends onto zeros, as scipy.signal.medfilt
  and MATLAB's medfilt1 default
- truncate: the window keeps only the samples inside x, as MATLAB's
  medfilt1(x, n, 'truncate'), so near the ends it's shorter and uneven

A window with an even count, which only truncate gives, takes the midpoint of
its two middle values, exactly, rounded once. A window holding a NaN gives NaN,
as MATLAB's default 'includenan'. Without NaN, zeropad is asserted to match
scipy.signal.medfilt exactly, which checks the definition.
"""

import math
import warnings
from fractions import Fraction
from typing import Final

import numpy as np
import scipy as sc

from dsp.cases import SIGNAL_LENGTH
from dsp.signals import synthetic_signals
from write import write_vector

NAN: Final = float("nan")
TINY: Final = 5e-324  # the smallest subnormal


def median(window: list[float]) -> float:
    """the exact median, rounded once: NaN if any value is NaN"""
    if any(math.isnan(v) for v in window):
        return NAN
    ordered = sorted(window)
    middle = len(ordered) // 2
    if len(ordered) % 2 == 1:
        return ordered[middle]
    return float((Fraction(ordered[middle - 1]) + Fraction(ordered[middle])) / 2)


def medfilt(x: list[float], n: int, edge: str) -> list[float]:
    """the median filter of odd length n, with the ends handled by edge"""
    half = n // 2
    y = []
    for i in range(len(x)):
        window = x[max(i - half, 0) : i + half + 1]
        if edge == "zeropad":
            window = window + [0.0] * (n - len(window))
        y.append(median(window))
    return y


def generate_medfilt() -> None:
    """the median filter, truncated and zero-padded at the ends"""
    signals = {name: x.tolist() for name, x in synthetic_signals(SIGNAL_LENGTH).items()}
    rng = np.random.default_rng(seed=30)
    t = np.arange(SIGNAL_LENGTH) / SIGNAL_LENGTH
    spikes = np.sin(2 * np.pi * 3 * t) + 20 + 0.1 * rng.standard_normal(SIGNAL_LENGTH)
    hits = rng.choice(SIGNAL_LENGTH, size=12, replace=False)
    spikes[hits] += rng.choice([-1.0, 1.0], size=12) * rng.uniform(5, 15, size=12)
    # the use case: spikes on a slow signal far from 0, which zero padding pulls down
    signals["offset_spikes"] = spikes.tolist()
    special = {
        # even-count windows at the ends: midpoints of the two middle values
        "short7": ([5.0, 1.0, 4.0, 2.0, 8.0, 3.0, 9.0], (3, 5)),
        # longer than x: truncate takes the whole of x near the middle
        "longer_window": ([3.0, 1.0, 2.0], (5, 9)),
        "one": ([2.5], (1, 3)),
        "empty": ([], (3,)),
        "ties": ([2.0, 2.0, 1.0, 2.0, 1.0, 1.0, 2.0], (3, 5)),
        # midpoints that overflow if formed as (a + b)/2
        "huge": ([1.7e308, 1.6e308, 1.5e308, -1.7e308, -1.6e308], (3, 5)),
        # midpoints of subnormals that round wrongly if formed as a/2 + b/2
        "subnormal": ([TINY, 2 * TINY, 3 * TINY, TINY, 2 * TINY], (3, 5)),
        "nan_middle": ([1.0, 2.0, 3.0, NAN, 5.0, 6.0, 7.0, 8.0, 9.0], (3, 5)),
        "nan_end": ([NAN, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], (3, 5)),
    }
    cases = [(name, x, n) for name, x in signals.items() for n in (1, 3, 5, 11)]
    cases += [(name, x, n) for name, (x, ns) in special.items() for n in ns]
    for name, x, n in cases:
        for edge in ("truncate", "zeropad"):
            y = medfilt(x, n, edge)
            if edge == "zeropad" and x and not any(math.isnan(v) for v in x):
                with warnings.catch_warnings():  # it warns of n longer than x
                    warnings.simplefilter("ignore", UserWarning)
                    assert y == sc.signal.medfilt(x, n).tolist(), (name, n)
            write_vector(
                "medfilt",
                f"{name}_n{n}_{edge}",
                f"the median filter, n = {n}, {edge} at the ends",
                {"x": x, "n": n, "edge": edge},
                {"y": y},
            )
