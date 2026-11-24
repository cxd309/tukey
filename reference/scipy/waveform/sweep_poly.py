"""
The exact polynomial sweep, at 50 significant digits with mpmath

scipy.signal.sweep_poly is cos(2·pi·phase + phi), where the frequency is the
polynomial p(t), its coefficients highest power first, and the phase in cycles
is its integral from 0:

- phase = sum of c·t^(k+1)/(k+1), for each coefficient c of t^k

As in chirp.py, the exact phase has its whole cycles removed before the cosine
is taken.

Each sweep is checked against its specification, independently of the formula:
the phase starts at 0, and its derivative is p(t).
The linear and quadratic sweeps are also chirp.py's exact chirps, so each
generator checks the other.
"""

from typing import Final

import mpmath as mp
import numpy as np
import scipy as sc

from waveform.chirp import DPS, GRIDS, SPEC_TOLERANCE, exact_chirp
from write import write_vector

# polynomial -> (chirp.py's grid of times, coefficients highest power first)
# constant: a 5 Hz tone
# linear and quadratic: chirp.py's unit sweeps from 1 to 80 Hz
# cubic: scipy.signal.sweep_poly's docs example, over ten seconds
# audio: chirp.py's linear sweep from 20 Hz to 20 kHz over ten seconds,
#   where the phase reaches 10^5 cycles
POLYS: Final = {
    "constant": ("unit", [5.0]),
    "linear": ("unit", [79.0, 1.0]),
    "quadratic": ("unit", [79.0, 0.0, 1.0]),
    "cubic": ("audio", [0.025, -0.36, 1.25, 2.0]),
    "audio": ("audio", [1998.0, 20.0]),
    "empty": ("empty", [79.0, 1.0]),
}

# how far scipy.signal.sweep_poly may be from the exact sweep: its rounding grows
# with the phase, at most 5.8e-14, 8.2e-14 and 1.2e-10
SCIPY_ACCURACY: Final = {
    "constant": 1e-13,
    "linear": 2e-13,
    "quadratic": 2e-13,
    "cubic": 2e-13,
    "audio": 1e-9,
    "empty": 0.0,
}

# the linear and quadratic polynomials as chirp.py's (f0, t1, f1, method)
AS_CHIRP: Final = {
    "linear": (1.0, 1.0, 80.0, "linear"),
    "quadratic": (1.0, 1.0, 80.0, "quadratic"),
}


def phase(t, coeffs):
    """the sweep's phase in cycles at t, the integral of the polynomial from 0, every argument an mpf"""
    degree = len(coeffs) - 1
    return mp.fsum(
        c * t ** (degree - i + 1) / (degree - i + 1) for i, c in enumerate(coeffs)
    )


def check_sweep(coeffs) -> None:
    """asserts the sweep meets its specification"""
    assert abs(phase(mp.mpf(0), coeffs)) < SPEC_TOLERANCE, coeffs
    for t in (mp.mpf("0.5"), mp.mpf(2)):
        frequency = mp.diff(lambda t: phase(t, coeffs), t)
        assert abs(frequency - mp.polyval(coeffs, t)) < SPEC_TOLERANCE, (coeffs, t)


def exact_sweep_poly(t, coeffs, phi) -> list[float]:
    """the sweep at every t, each value exact rounded once to float64"""
    with mp.workdps(DPS):
        coeffs = [mp.mpf(c) for c in coeffs]
        check_sweep(coeffs)
        offset = mp.mpf(phi) / 360
        return [
            float(mp.cospi(2 * mp.frac(phase(mp.mpf(float(ti)), coeffs) + offset)))
            for ti in t
        ]


def cases():
    """(name, poly, phi) of every sweep: each at phi 0, and linear at 90 and -45"""
    for name in POLYS:
        yield name, name, 0.0
    for phi in (90.0, -45.0):
        yield f"linear_phi{phi:g}", "linear", phi


def generate_sweep_poly() -> None:
    """scipy.signal.sweep_poly's exact output, after checking scipy is within SCIPY_ACCURACY"""
    for name, poly, phi in cases():
        grid, coeffs = POLYS[poly]
        t = GRIDS[grid][0]
        x = exact_sweep_poly(t, coeffs, phi)
        if poly in AS_CHIRP:
            f0, t1, f1, method = AS_CHIRP[poly]
            assert x == exact_chirp(t, f0, t1, f1, method, True, phi), name
        got = sc.signal.sweep_poly(t, coeffs, phi)
        assert np.all(np.abs(got - x) <= SCIPY_ACCURACY[poly]), (
            name,
            np.max(np.abs(got - x)),
        )
        write_vector(
            "waveform/sweep_poly",
            name,
            f"scipy.signal.sweep_poly(t, {coeffs}, {phi:g}), exact",
            {"t": t.tolist(), "poly": coeffs, "phi": phi},
            {"x": x},
        )
