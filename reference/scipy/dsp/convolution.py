"""
Convolution and correlation, in full, same and valid modes, with each
correlation output's lag

The references are exact: each output is the exact sum of exact products
(Dekker's two-product, summed with math.fsum), rounded once. tukey offers
direct and FFT methods, so both are judged against the true answer rather than
one method's rounding. Each generator asserts scipy.signal's direct and FFT
results agree with it to within their own rounding, which checks the
definitions and the mode slicing.

Lags are the lags each output actually has, by the definition
c[k] = sum a[n + k]·v[n]. SciPy 1.18.1's correlation_lags(mode="same") is one
lag off correlate(mode="same") when len(a) is odd and len(v) is even, so it's
checked only where the two agree.
"""

import math
from typing import Final

import numpy as np
import scipy as sc

from write import write_vector

SPLIT: Final = 134217729.0  # 2^27 + 1, Dekker's splitter for float64


def split(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """x as hi + lo exactly, each with at most 26 significant bits"""
    c = SPLIT * x
    hi = c - (c - x)
    return hi, x - hi


def two_product(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """x·y as p + e exactly, p the rounded product (Dekker), elementwise"""
    p = x * y
    xh, xl = split(x)
    yh, yl = split(y)
    e = ((xh * yh - p) + xh * yl + xl * yh) + xl * yl
    return p, e


def exact_convolve(a: np.ndarray, v: np.ndarray) -> list[float]:
    """the full convolution, each output sum(a[i]·v[k - i]) exact, rounded once"""
    p, e = two_product(a[:, None], v[None, :])
    # output k sums the anti-diagonal i + j = k: flip v's axis to make it a diagonal
    p, e = p[:, ::-1], e[:, ::-1]
    nv = len(v)
    return [
        math.fsum(
            np.diagonal(p, offset=nv - 1 - k).tolist()
            + np.diagonal(e, offset=nv - 1 - k).tolist()
        )
        for k in range(len(a) + nv - 1)
    ]


def same_slice(na: int, nv: int) -> slice:
    """the central len(a) outputs of the full result, as SciPy takes them"""
    start = (nv - 1) // 2
    return slice(start, start + na)


def valid_slice(na: int, nv: int) -> slice:
    """
    the outputs where the shorter input lies wholly inside the longer: SciPy
    swaps the inputs when v is the longer, where MATLAB's 'valid' is empty
    """
    return slice(min(na, nv) - 1, max(na, nv))


def convolution_cases() -> list[tuple[str, np.ndarray, np.ndarray]]:
    """(name, a, v) for each vector, shared by convolve and correlate"""
    rng = np.random.default_rng(seed=31)
    spike = 1e-3 * rng.standard_normal(1000)
    spike[500] = 1e8  # an FFT's error is relative to the largest output
    return [
        ("integers", np.array([1.0, 2, 3, 4, 5]), np.array([1.0, 10, 100, 1000])),
        ("odd_kernel", rng.standard_normal(50), np.array([0.25, 0.5, 0.25])),
        ("even_kernel", rng.standard_normal(51), rng.standard_normal(6)),
        ("kernel_longer", rng.standard_normal(5), rng.standard_normal(12)),
        ("equal_lengths", rng.standard_normal(16), rng.standard_normal(16)),
        ("one_sample_kernel", rng.standard_normal(20), np.array([2.5])),
        ("one_sample_each", np.array([3.0]), np.array([-2.0])),
        # a realistic FIR filter on a long signal, where FFT starts to pay
        ("fir_lowpass", rng.standard_normal(2000), sc.signal.firwin(101, 0.2)),
        ("long_both", rng.standard_normal(1500), rng.standard_normal(1000)),
        ("mixed_magnitudes", spike, rng.standard_normal(64)),
    ]


def assert_close(name: str, mode: str, got: np.ndarray, exact: list[float]) -> None:
    """got agrees with exact to within rounding at the scale of the result"""
    assert len(got) == len(exact), (name, mode)
    scale = max(1.0, float(np.abs(exact).max())) if exact else 1.0
    assert np.allclose(got, exact, rtol=0, atol=1e-10 * scale), (name, mode)


def generate_convolve() -> None:
    """the exact convolution, in full, same and valid modes"""
    for name, a, v in convolution_cases():
        full = exact_convolve(a, v)
        modes = {
            "full": full,
            "same": full[same_slice(len(a), len(v))],
            "valid": full[valid_slice(len(a), len(v))],
        }
        for mode, exact in modes.items():
            for method in ("direct", "fft"):
                got = sc.signal.convolve(a, v, mode=mode, method=method)
                assert_close(name, f"{mode} {method}", got, exact)
        write_vector(
            "convolve",
            name,
            "the exact convolution of a and v; same and valid as scipy.signal.convolve",
            {"a": a.tolist(), "v": v.tolist()},
            modes,
        )


def generate_correlate() -> None:
    """the exact correlation, c[k] = sum a[n + k]·v[n], with each output's lag"""
    for name, a, v in convolution_cases():
        na, nv = len(a), len(v)
        full = exact_convolve(a, v[::-1])
        lags = list(range(-(nv - 1), na))
        output: dict[str, list] = {}
        for mode, part in (
            ("full", slice(None)),
            ("same", same_slice(na, nv)),
            ("valid", valid_slice(na, nv)),
        ):
            exact, mode_lags = full[part], lags[part]
            for method in ("direct", "fft"):
                got = sc.signal.correlate(a, v, mode=mode, method=method)
                assert_close(name, f"{mode} {method}", got, exact)
            scipy_lags = sc.signal.correlation_lags(na, nv, mode=mode).tolist()
            known_off_by_one = mode == "same" and na % 2 == 1 and nv % 2 == 0
            assert (scipy_lags == mode_lags) != known_off_by_one, (name, mode)
            output[mode] = exact
            output[f"{mode}_lags"] = mode_lags
        write_vector(
            "correlate",
            name,
            "the exact correlation sum a[n + k]·v[n]; same and valid as "
            "scipy.signal.correlate, and the lag k of each output",
            {"a": a.tolist(), "v": v.tolist()},
            output,
        )
