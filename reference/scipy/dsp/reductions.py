"""
Reductions: sum and mean, differences and gradients, running sums and
trapezoidal integrals, and percentiles

tukey sums with compensation, so the sums here are exact (math.fsum, and
Shewchuk's partials for running sums) rather than NumPy's: each generator
asserts NumPy agrees to well within its own rounding, which checks the
definition, and writes the exact value as the reference.
"""

import math
from collections.abc import Callable
from typing import Final

import numpy as np
import scipy as sc

from write import write_vector

NAN: Final = float("nan")
LONG: Final = 20000  # long enough for a plain sum's rounding to build up


def summed_signals() -> dict[str, np.ndarray]:
    """signals for the sums, from exact to hard for a plain sum"""
    rng = np.random.default_rng(seed=21)
    return {
        "ramp": np.arange(1.0, 11.0),  # small integers: every sum is exact
        "noise": rng.standard_normal(1000),
        "tenths_long": np.full(LONG, 0.1),  # 0.1 isn't exact, so a plain sum drifts
        # one large value among many small ones, then cancelled: a plain sum
        # loses the small values' low bits against the large one
        "mixed_magnitudes": np.concatenate(([1e8], np.full(LONG, 1e-3), [-1e8])),
        "one": np.array([2.5]),
        "empty": np.array([]),
    }


def running_exact(
    x: np.ndarray, extra: Callable[[int], list[float]] = lambda i: []
) -> list[float]:
    """
    each prefix's exact sum, x[0] + ... + x[i] + sum(extra(i)), rounded once,
    keeping Shewchuk's exact partials as it goes, as math.fsum does internally
    """
    partials: list[float] = []
    sums = []
    for index, v in enumerate(x.tolist()):
        i = 0
        for y in partials:
            if abs(v) < abs(y):
                v, y = y, v
            high = v + y
            low = y - (high - v)
            if low:
                partials[i] = low
                i += 1
            v = high
        partials[i:] = [v]
        sums.append(math.fsum(partials + extra(index)))
    return sums


def generate_sum() -> None:
    """the exact sum, which numpy.sum's pairwise summation is close to"""
    for name, x in summed_signals().items():
        exact = math.fsum(x.tolist())
        assert math.isclose(np.sum(x), exact, rel_tol=1e-12, abs_tol=1e-6), name
        write_vector(
            "sum",
            name,
            "math.fsum(x), the exact sum, rounded once",
            {"x": x.tolist()},
            {"sum": [exact]},
        )


def generate_mean() -> None:
    """the exact sum divided by the length. An empty signal has no mean"""
    for name, x in summed_signals().items():
        if len(x) == 0:
            continue
        mean = math.fsum(x.tolist()) / len(x)
        assert math.isclose(np.mean(x), mean, rel_tol=1e-12, abs_tol=1e-10), name
        write_vector(
            "mean",
            name,
            "math.fsum(x) / len(x)",
            {"x": x.tolist()},
            {"mean": [mean]},
        )


def generate_cumsum() -> None:
    """each prefix's exact sum, which numpy.cumsum sums plainly"""
    for name, x in summed_signals().items():
        exact = running_exact(x)
        if len(x) > 0:
            scale = max(1.0, np.abs(np.cumsum(np.abs(x))).max())
            assert np.allclose(np.cumsum(x), exact, rtol=0, atol=1e-9 * scale), name
        write_vector(
            "cumsum",
            name,
            "each prefix of x summed exactly",
            {"x": x.tolist()},
            {"s": exact},
        )


def generate_diff() -> None:
    """numpy.diff: one subtraction each, so exact"""
    rng = np.random.default_rng(seed=22)
    cases = {
        "squares": np.array([1.0, 4, 9, 16, 25]),
        "noise": rng.standard_normal(500),
        "two": np.array([3.0, 1.0]),
        "one": np.array([2.0]),  # no gaps: empty
        "empty": np.array([]),
    }
    for name, x in cases.items():
        write_vector(
            "diff", name, "numpy.diff(x)", {"x": x.tolist()}, {"d": np.diff(x).tolist()}
        )


def generate_gradient() -> None:
    """numpy.gradient: central differences inside, one-sided at the ends"""
    rng = np.random.default_rng(seed=23)
    t = np.arange(400) * 0.01
    cases = [
        ("squares", np.array([1.0, 4, 9, 16]), 1.0),
        ("noise", rng.standard_normal(500), 0.01),
        ("sine", np.sin(2 * np.pi * t), 0.01),
        ("two", np.array([3.0, 1.0]), 0.5),  # the fewest a slope needs
    ]
    for name, x, dx in cases:
        write_vector(
            "gradient",
            name,
            "numpy.gradient(x, dx)",
            {"x": x.tolist(), "dx": dx},
            {"g": np.gradient(x, dx).tolist()},
        )


def trapezoid_cases() -> list[tuple[str, np.ndarray, float]]:
    """signals and sample spacings for both trapezoidal integrals"""
    rng = np.random.default_rng(seed=24)
    t = np.arange(1000) * 0.001
    signals = summed_signals()
    return [
        ("noise", rng.standard_normal(1000), 0.01),
        ("sine_period", np.sin(2 * np.pi * t), 0.001),  # integrates to about 0
        ("ramp", signals["ramp"], 0.5),
        ("tenths_long", signals["tenths_long"], 0.25),
        ("mixed_magnitudes", signals["mixed_magnitudes"], 1.0),
        ("one", signals["one"], 1.0),  # no interval: 0
        ("empty", signals["empty"], 1.0),
    ]


def generate_trapezoid() -> None:
    """dx·(y[0]/2 + y[1] + ... + y[n-1]/2), the sum exact and multiplied by dx once"""
    for name, y, dx in trapezoid_cases():
        terms = y.tolist()
        if len(terms) >= 2:
            terms[0] /= 2
            terms[-1] /= 2
        exact = dx * math.fsum(terms) if len(terms) >= 2 else 0.0
        reference = sc.integrate.trapezoid(y, dx=dx) if len(y) >= 2 else 0.0
        assert math.isclose(reference, exact, rel_tol=1e-9, abs_tol=1e-9), name
        write_vector(
            "trapezoid",
            name,
            "dx * exact sum of y[0]/2, y[1], ..., y[n-1]/2",
            {"y": y.tolist(), "dx": dx},
            {"integral": [exact]},
        )


def generate_cumulative_trapezoid() -> None:
    """
    the running trapezoidal integral, starting at 0 and the same length as y,
    each value dx times the exact sum of y[0]/2, y[1], ..., y[i]/2
    """
    for name, y, dx in trapezoid_cases():
        exact = []
        if len(y) > 0:
            # y[0] + ... + y[i], less y[0]/2 and y[i]/2, summed exactly
            halves = running_exact(y, lambda i, y=y: [-y[0] / 2, -y[i] / 2])
            exact = [0.0] + [dx * total for total in halves[1:]]
            reference = sc.integrate.cumulative_trapezoid(y, dx=dx, initial=0)
            scale = max(1.0, dx * np.abs(y).sum())
            assert np.allclose(reference, exact, rtol=0, atol=1e-9 * scale), name
        write_vector(
            "cumulative_trapezoid",
            name,
            "dx * exact sum of y[0]/2, y[1], ..., y[i]/2, starting at 0",
            {"y": y.tolist(), "dx": dx},
            {"integral": exact},
        )


def generate_percentile() -> None:
    """numpy.percentile, linear (its default) and hazen"""
    rng = np.random.default_rng(seed=25)
    cases = [
        ("four", np.array([1.0, 2, 3, 4]), [0, 10, 25, 50, 75, 90, 100]),
        ("noise101", rng.standard_normal(101), [0, 1, 5, 25, 50, 95, 99, 100]),
        ("noise100", rng.standard_normal(100), [0, 1, 5, 25, 50, 95, 99, 100]),
        ("ties", np.array([3.0, 1, 2, 2, 1, 2]), [0, 20, 50, 80, 100]),
        ("one", np.array([7.0]), [0, 50, 100]),
        ("nan", np.array([1.0, NAN, 3, 2]), [50]),  # any NaN gives NaN
    ]
    for name, x, ps in cases:
        for method in ("linear", "hazen"):
            values = np.percentile(x, ps, method=method)
            write_vector(
                "percentile",
                f"{name}_{method}",
                f"numpy.percentile(x, p, method={method!r})",
                {"x": x.tolist(), "p": [float(p) for p in ps], "method": method},
                {"values": values.tolist()},
            )
