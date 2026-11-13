"""
Smoothing: the median filter and the Savitzky–Golay filter

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

The Savitzky–Golay filter is computed exactly too. Each output is a polynomial
of degree polyorder, fitted by least squares to a window of n samples, or its
deriv-th derivative divided by delta^deriv, as scipy.signal.savgol_filter's
default mode='interp':

- inside, the window is centred on the sample
- within n//2 of the ends, the window is the first or last n samples, and
  the one fit is evaluated at each edge sample

The fit to integer positions has a rational solution, so each output's
weights are exact Fractions, and the output, a sum of weights times x, is
exact and rounded once. A window holding a NaN gives NaN. scipy's output is
asserted to agree, to within its own rounding, which checks the definition.

scipy fits the raw powers of the positions by lstsq, so its accuracy falls with
the window and degree: 1e-12 at n = 21, 3e-9 at n = 51, degree 6, and at
n = 101, degree 10, where the powers reach 1e17, lstsq drops the small singular
values and its coefficients are nonsense, wrong by 100%. Those vectors are
checked against mpmath instead: least squares at 80 digits, in the middle and
at an edge.
"""

import math
import warnings
from fractions import Fraction
from functools import cache
from typing import Final

import mpmath as mp
import numpy as np
import scipy as sc

from dsp.cases import SIGNAL_LENGTH
from dsp.signals import synthetic_signals
from write import write_vector

NAN: Final = float("nan")
TINY: Final = 5e-324  # the smallest subnormal
FLOAT_SHIFT: Final = 1074  # 2^1074·x is an integer for every finite float64
# how far scipy.signal.savgol_filter may be from the exact output, relative to
# max |x| / delta^deriv: its rounding, which grows with the window and degree.
# None where scipy can't check the case, which mpmath or MATLAB checks instead
SAVGOL_SCIPY_ACCURACY: Final = 1e-11
SAVGOL_SCIPY_ACCURACY_BY_CASE: Final = {
    "long51": 1e-8,
    "long101": None,  # scipy's fit fails
    "nan_end": None,  # scipy's edge fit rejects a NaN
}
SAVGOL_MPMATH_DPS: Final = 80  # the independent check's working precision


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


@cache
def fit_matrix(n: int, polyorder: int) -> tuple[tuple[Fraction, ...], ...]:
    """
    (AᵀA)⁻¹Aᵀ for the Vandermonde matrix A of positions 0 to n - 1 and powers
    0 to polyorder, exactly: row j maps a window's samples to the fitted
    polynomial's coefficient of t^j
    """
    powers = polyorder + 1
    a = [[Fraction(t) ** j for j in range(powers)] for t in range(n)]
    # Gauss–Jordan on [AᵀA | Aᵀ]: exact, so no pivoting for accuracy is needed
    rows = [
        [sum(a[t][i] * a[t][j] for t in range(n)) for j in range(powers)]
        + [a[t][i] for t in range(n)]
        for i in range(powers)
    ]
    for col in range(powers):
        pivot = next(r for r in range(col, powers) if rows[r][col] != 0)
        rows[col], rows[pivot] = rows[pivot], rows[col]
        lead = rows[col][col]
        rows[col] = [v / lead for v in rows[col]]
        for r in range(powers):
            if r != col and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [
                    v - factor * p for v, p in zip(rows[r], rows[col], strict=True)
                ]
    return tuple(tuple(row[powers:]) for row in rows)


@cache
def savgol_weights(
    n: int, polyorder: int, deriv: int, at: int
) -> tuple[tuple[int, ...], int]:
    """
    (weights, denominator): the fitted polynomial's deriv-th derivative at
    position at, as integer weights on the window's samples over one common
    denominator
    """
    fit = fit_matrix(n, polyorder)
    weights = [Fraction(0)] * n
    for j in range(deriv, polyorder + 1):
        # d^deriv/dt^deriv t^j at t = at
        factor = Fraction(math.perm(j, deriv)) * Fraction(at) ** (j - deriv)
        weights = [w + factor * f for w, f in zip(weights, fit[j], strict=True)]
    denominator = math.lcm(*(w.denominator for w in weights))
    return tuple(int(w * denominator) for w in weights), denominator


def savgol_filter(
    x: list[float], n: int, polyorder: int, deriv: int, delta: float
) -> list[float]:
    """the exact Savitzky–Golay filter, mode='interp', each output rounded once"""
    half = n // 2
    scaled = [None if math.isnan(v) else int(Fraction(v) * 2**FLOAT_SHIFT) for v in x]
    scale = Fraction(delta) ** deriv * 2**FLOAT_SHIFT

    def output(start: int, at: int) -> float:
        window = scaled[start : start + n]
        if None in window:
            return NAN
        weights, denominator = savgol_weights(n, polyorder, deriv, at)
        total = sum(w * v for w, v in zip(weights, window, strict=True))
        return float(Fraction(total, denominator) / scale)

    y = []
    for i in range(len(x)):
        if i < half:
            y.append(output(0, i))
        elif i >= len(x) - half:
            y.append(output(len(x) - n, i - (len(x) - n)))
        else:
            y.append(output(i - half, half))
    return y


def fit_by_mpmath(
    x: np.ndarray, n: int, polyorder: int, deriv: int, delta: float, i: int
) -> float:
    """
    output i of the Savitzky–Golay filter, independently: least squares by
    mpmath's LU on the normal equations, with the window's positions centred
    on output i, at SAVGOL_MPMATH_DPS digits
    """
    with mp.workdps(SAVGOL_MPMATH_DPS):
        start = min(max(i - n // 2, 0), len(x) - n)
        ts = [mp.mpf(t - i) for t in range(start, start + n)]
        a = mp.matrix([[t**j for j in range(polyorder + 1)] for t in ts])
        ys = mp.matrix([mp.mpf(float(v)) for v in x[start : start + n]])
        coeffs = mp.lu_solve(a.T * a, a.T * ys)
        # at t = 0 the deriv-th derivative is deriv!·coefficient of t^deriv
        return float(mp.factorial(deriv) * coeffs[deriv] / mp.mpf(delta) ** deriv)


def generate_savgol_filter() -> None:
    """the Savitzky–Golay filter, smoothing and differentiating, fitted at the ends"""
    signals = synthetic_signals(SIGNAL_LENGTH)
    rng = np.random.default_rng(seed=31)
    noise = rng.standard_normal(SIGNAL_LENGTH)
    k = np.arange(50.0)
    cubic = 0.5 * k**3 - 3 * k**2 + k - 2  # exact: a cubic fit reproduces it
    nan_inside, nan_end = noise.copy(), noise.copy()
    nan_inside[100] = NAN
    nan_end[3] = NAN  # in the first fit window: every edge sample is NaN

    # (name, x, n, polyorder, deriv, delta)
    cases = [
        (name, x, n, p, 0, 1.0)
        for name, x in signals.items()
        for n, p in ((5, 2), (11, 3), (21, 4))
    ]
    # derivatives, with the sample spacing of a signal of 200 samples over 1 s
    cases += [
        (name, signals[name], n, p, d, 1 / SIGNAL_LENGTH)
        for name in ("sines", "chirp", "noise")
        for n, p, d in ((11, 3, 1), (11, 4, 2))
    ]
    cases += [("cubic", cubic, 7, 3, d, 1.0) for d in range(5)]  # d = 4: zero
    cases += [
        ("moving_average", noise, 5, 0, 0, 1.0),
        ("identity", noise, 5, 4, 0, 1.0),  # polyorder n - 1: the fit is exact
        ("n1", noise, 1, 0, 0, 1.0),
        ("whole", noise[:9], 9, 2, 0, 1.0),  # one window: all of x
        ("short", noise[:3], 3, 1, 0, 1.0),
        ("offset", 1e6 + 0.01 * noise, 11, 3, 0, 1.0),  # far from 0
        ("long51", noise, 51, 6, 0, 1.0),
        ("long101", noise, 101, 10, 0, 1.0),  # high degree: the fit is ill-conditioned
        ("long101", noise, 101, 10, 1, 0.1),
        ("nan", nan_inside, 11, 3, 0, 1.0),
        ("nan_end", nan_end, 11, 3, 0, 1.0),
    ]
    for name, x, n, p, d, delta in cases:
        x = np.asarray(x, dtype=float)
        y = savgol_filter(x.tolist(), n, p, d, delta)
        accuracy = SAVGOL_SCIPY_ACCURACY_BY_CASE.get(name, SAVGOL_SCIPY_ACCURACY)
        if accuracy is not None:
            got = sc.signal.savgol_filter(x, n, p, deriv=d, delta=delta)
            scale = np.nanmax(np.abs(x)) / delta**d
            diff = np.nanmax(np.abs(got - np.array(y))) / scale
            assert diff < accuracy, (name, n, p, d, diff)
            assert np.array_equal(np.isnan(got), np.isnan(y)), name
        elif name == "long101":
            for i in (0, len(x) // 2):
                assert y[i] == fit_by_mpmath(x, n, p, d, delta, i), (name, d, i)
        write_vector(
            "savgol_filter",
            f"{name}_n{n}_p{p}" + (f"_d{d}" if d else ""),
            f"the Savitzky–Golay filter, n = {n}, polyorder = {p}, deriv = {d}, mode='interp'",
            {"x": x.tolist(), "n": n, "polyorder": p, "deriv": d, "delta": delta},
            {"y": y},
        )
