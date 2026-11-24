"""
The exact chirp, at 50 significant digits with mpmath

scipy.signal.chirp is cos(2·pi·phase + phi), with phase in cycles:

- linear: f0·t + (f1 - f0)·t²/(2·t1)
- quadratic, vertex at 0: f0·t + (f1 - f0)·t³/(3·t1²)
- quadratic, vertex at t1: f1·t + (f1 - f0)·((t1 - t)³ - t1³)/(3·t1²)
- logarithmic: t1·f0·((f1/f0)^(t/t1) - 1)/ln(f1/f0), or f0·t when f0 = f1
- hyperbolic: -s·f0·ln|1 - t/s| with s = f1·t1/(f1 - f0), or f0·t when f0 = f1

scipy multiplies by 2·pi before taking the cosine, so at large t the phase
loses digits.
The exact phase is reduced to its fractional cycle first, then the cosine taken.

Each sweep is checked against its specification, independently of the
formulas: the phase starts at 0, and the frequency, its derivative, is f0 at
t = 0 and f1 at t = t1.
"""

from typing import Final

import mpmath as mp
import numpy as np
import scipy as sc

from write import write_vector

DPS: Final = 50  # working precision, in significant digits
SPEC_TOLERANCE: Final = mp.mpf("1e-30")  # how far a sweep may miss its spec

# the sample times, and the sweep over them
# unit: one second at 200 Hz, 1 to 80 Hz, as dsp's synthetic chirp
# beyond: before 0 and past t1, where each sweep carries on,
#   and the hyperbolic up sweep passes its singularity at t = 80/79
# audio: ten seconds sweeping 20 Hz to 20 kHz, at 1001 of its times,
#   where the phase reaches 10^5 cycles
# empty: no samples
GRIDS: Final = {
    "unit": (np.arange(201) / 200, 1.0, 1.0, 80.0),
    "beyond": (np.arange(-100, 301) / 100, 1.0, 1.0, 80.0),
    "audio": (np.linspace(0, 10, 1001), 20.0, 10.0, 20000.0),
    "empty": (np.array([]), 1.0, 1.0, 80.0),
}

# (method, vertex_zero)
SWEEPS: Final = [
    ("linear", True),
    ("quadratic", True),
    ("quadratic", False),
    ("logarithmic", True),
    ("hyperbolic", True),
]

# how far scipy.signal.chirp may be from the exact chirp: its rounding grows
# with the phase, so it's measured per grid: at most 1.2e-13, 8.5e-11 and 1.8e-10
SCIPY_ACCURACY: Final = {
    "unit": 2e-13,
    "beyond": 1e-10,
    "audio": 1e-9,
    "empty": 0.0,
}


def phase(t, f0, t1, f1, method: str, vertex_zero: bool):
    """the chirp's phase in cycles at t, every argument an mpf"""
    if method == "linear":
        return f0 * t + (f1 - f0) * t**2 / (2 * t1)
    if method == "quadratic":
        beta = (f1 - f0) / t1**2
        if vertex_zero:
            return f0 * t + beta * t**3 / 3
        return f1 * t + beta * ((t1 - t) ** 3 - t1**3) / 3
    if f0 == f1:
        return f0 * t
    if method == "logarithmic":
        return t1 * f0 * (mp.power(f1 / f0, t / t1) - 1) / mp.log(f1 / f0)
    s = f1 * t1 / (f1 - f0)
    return -s * f0 * mp.log(abs(1 - t / s))


def check_sweep(f0, t1, f1, method: str, vertex_zero: bool) -> None:
    """asserts the sweep meets its specification"""
    name = (f0, t1, f1, method, vertex_zero)

    def cycles(t):
        return phase(t, f0, t1, f1, method, vertex_zero)

    assert abs(cycles(mp.mpf(0))) < SPEC_TOLERANCE, name
    assert abs(mp.diff(cycles, mp.mpf(0)) - f0) < SPEC_TOLERANCE, name
    assert abs(mp.diff(cycles, t1) - f1) < SPEC_TOLERANCE, name


def exact_chirp(t, f0, t1, f1, method, vertex_zero, phi) -> list[float]:
    """the chirp at every t, each value exact rounded once to float64"""
    with mp.workdps(DPS):
        f0, t1, f1 = mp.mpf(f0), mp.mpf(t1), mp.mpf(f1)
        check_sweep(f0, t1, f1, method, vertex_zero)
        offset = mp.mpf(phi) / 360
        x = []
        for ti in t:
            cycles = phase(mp.mpf(float(ti)), f0, t1, f1, method, vertex_zero)
            x.append(float(mp.cospi(2 * mp.frac(cycles + offset))))
        return x


def cases():
    """(name, grid, f0, f1, method, vertex_zero, phi) of every chirp"""
    for method, vertex_zero in SWEEPS:
        sweep = method if vertex_zero else f"{method}_vertex_t1"
        for grid, (_, f0, _, f1) in GRIDS.items():
            yield f"{sweep}_{grid}", grid, f0, f1, method, vertex_zero, 0.0
        _, f0, _, f1 = GRIDS["unit"]
        # down: the same sweep from 80 to 1 Hz
        yield f"{sweep}_unit_down", "unit", f1, f0, method, vertex_zero, 0.0
        # phi: 90 degrees is a sine, -45 is between
        for phi in (90.0, -45.0):
            yield f"{sweep}_unit_phi{phi:g}", "unit", f0, f1, method, vertex_zero, phi
        # constant: f0 = f1, a special case for logarithmic and hyperbolic
        yield f"{sweep}_unit_constant", "unit", 5.0, 5.0, method, vertex_zero, 0.0


def generate_chirp() -> None:
    """scipy.signal.chirp's exact output, after checking scipy is within SCIPY_ACCURACY"""
    for name, grid, f0, f1, method, vertex_zero, phi in cases():
        t, _, t1, _ = GRIDS[grid]
        x = exact_chirp(t, f0, t1, f1, method, vertex_zero, phi)
        got = sc.signal.chirp(t, f0, t1, f1, method, phi, vertex_zero)
        assert np.all(np.abs(got - x) <= SCIPY_ACCURACY[grid]), (
            name,
            np.max(np.abs(got - x)),
        )
        write_vector(
            "waveform/chirp",
            name,
            f"scipy.signal.chirp(t, {f0:g}, {t1:g}, {f1:g}, {method!r}, {phi:g}, "
            f"vertex_zero={vertex_zero}), exact",
            {
                "t": t.tolist(),
                "f0": f0,
                "t1": t1,
                "f1": f1,
                "method": method,
                "phi": phi,
                "vertex_zero": vertex_zero,
            },
            {"x": x},
        )
