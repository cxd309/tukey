"""scipy.fft transforms, their norm and n options, and the frequency and shift helpers"""

import numpy as np
import scipy as sc

from .cases import FFT_LENGTHS
from .write import complex_json, write_vector


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


def generate_fft_options() -> None:
    """scipy.fft's norm and n options on each transform, and its frequency and shift helpers."""
    for n in (12, 1009):  # mixed radix and Bluestein
        rng = np.random.default_rng(seed=n + 1)
        x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        spectrum = sc.fft.rfft(x.real)
        for norm in ("backward", "ortho", "forward"):
            write_fft_option(
                f"fft_n{n}_{norm}", "fft", x, sc.fft.fft(x, norm=norm), norm, None
            )
            write_fft_option(
                f"ifft_n{n}_{norm}", "ifft", x, sc.fft.ifft(x, norm=norm), norm, None
            )
            write_fft_option(
                f"rfft_n{n}_{norm}",
                "rfft",
                x.real,
                sc.fft.rfft(x.real, norm=norm),
                norm,
                None,
            )
            write_fft_option(
                f"irfft_n{n}_{norm}",
                "irfft",
                spectrum,
                sc.fft.irfft(spectrum, n=n, norm=norm),
                norm,
                n,
            )

    # n pads (16) or truncates (8) a length-12 input before transforming
    x = np.random.default_rng(seed=12).standard_normal(12) + 0j
    for n in (8, 16):
        write_fft_option(
            f"fft_n{n}_from12", "fft", x, sc.fft.fft(x, n=n), "backward", n
        )
        write_fft_option(
            f"ifft_n{n}_from12", "ifft", x, sc.fft.ifft(x, n=n), "backward", n
        )
        write_fft_option(
            f"rfft_n{n}_from12",
            "rfft",
            x.real,
            sc.fft.rfft(x.real, n=n),
            "backward",
            n,
        )

    for n in (1, 8, 9):  # smallest, even, odd
        for d in (1.0, 0.001):
            write_vector(
                "fftfreq",
                f"n{n}_d{d}",
                f"scipy.fft.fftfreq({n}, {d}) and rfftfreq({n}, {d})",
                {"n": n, "d": d},
                {
                    "freq": sc.fft.fftfreq(n, d).tolist(),
                    "rfreq": sc.fft.rfftfreq(n, d).tolist(),
                },
            )
        x = np.arange(n, dtype=float)
        write_vector(
            "fftshift",
            f"n{n}",
            f"scipy.fft.fftshift and ifftshift of arange({n})",
            {"x": x.tolist()},
            {
                "shift": sc.fft.fftshift(x).tolist(),
                "ishift": sc.fft.ifftshift(x).tolist(),
            },
        )


def write_fft_option(name: str, transform: str, x, y, norm: str, n: int | None) -> None:
    """One option case; real inputs and outputs are stored as complex, so every case shares a shape."""
    write_vector(
        "fft_options",
        name,
        f"scipy.fft.{transform}(x, n={n!r}, norm={norm!r})",
        {"transform": transform, "x": complex_json(x), "norm": norm, "n": n},
        {"y": complex_json(y)},
    )
