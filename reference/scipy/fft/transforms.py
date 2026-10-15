"""scipy.fft transforms, their norm and n options, and the frequency and shift helpers"""

import numpy as np
import scipy as sc

from write import complex_json, write_vector

# lengths the FFT must handle, chosen to exercise each algorithm path
FFT_LENGTHS = (
    1,
    2,
    4,
    8,
    16,
    128,
    1024,
    4096,  # powers of two: radix-2
    3,
    5,
    9,
    15,
    30,
    60,
    100,
    1000,  # smooth: mixed radix over 2, 3, 5
    6,
    12,
    210,  # more mixed factors; 210 = 2*3*5*7 brings in a 7
    7,
    17,
    97,
    1009,  # primes: the generic butterfly or Bluestein
)


def generate_fft() -> None:
    """scipy.fft transforms of seeded random signals, at every length in FFT_LENGTHS."""
    for n in FFT_LENGTHS:
        rng = np.random.default_rng(seed=n)  # per-length seed
        x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        real = rng.standard_normal(n)
        spectrum = sc.fft.rfft(rng.standard_normal(n))  # Hermitian as irfft expects

        write_vector(
            "fft",
            f"n{n}",
            f"scipy.fft.fft(x), n={n}",
            {"x": complex_json(x)},
            {"y": complex_json(sc.fft.fft(x))},
        )
        write_vector(
            "ifft",
            f"n{n}",
            f"scipy.fft.ifft(x), n={n}",
            {"x": complex_json(x)},
            {"y": complex_json(sc.fft.ifft(x))},
        )
        write_vector(
            "rfft",
            f"n{n}",
            f"scipy.fft.rfft(x), n={n}",
            {"x": real.tolist()},
            {"y": complex_json(sc.fft.rfft(real))},
        )
        write_vector(
            "irfft",
            f"n{n}",
            f"scipy.fft.irfft(x, n={n})",
            {"x": complex_json(spectrum), "n": n},
            {"y": sc.fft.irfft(spectrum, n=n).tolist()},
        )
