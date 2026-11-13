"""
The exact group delay, at 50 significant digits with mpmath

For a filter B(z)/A(z), with B and A polynomials in z^-1 in the coefficients as
given, the group delay in samples at w rad/sample is

    Re(Σ k·b_k·e^(-jwk) / B) - Re(Σ k·a_k·e^(-jwk) / A)

evaluated here at freqz's grid, w = pi·i/n, from the coefficients exactly.

Where B or A is zero on the unit circle, at a grid point, the phase jumps by pi
and the delay is undefined. Either side of it, each such zero adds exactly half
a sample, so the delay has a limit from both sides, which the vectors give,
with the point's index listed as singular: MATLAB's grpdelay returns the limit,
SciPy's group_delay returns 0, and a caller may want NaN.
"""

from typing import Final

import mpmath as mp

from dsp.elliptic import DPS

SINGULAR: Final = mp.mpf("1e-40")  # |B|/Σ|b| below which B is zero at a point
# either side of a singular point, for its limit: far enough that a double
# zero leaves |B| near 1e-24, so 26 digits, and the limit is good to 1e-24
STEP: Final = mp.mpf("1e-12")


def polynomial_delay(
    coeffs: list, w: mp.mpf, check: bool = True
) -> tuple[mp.mpf | None, bool]:
    """
    (Re(Σ k·c_k·x^k / Σ c_k·x^k) at x = e^(-jw), whether it's singular
    there), or without check, the value even where it's nearly singular
    """
    x = mp.expj(-w)
    value = mp.fsum(c * x**k for k, c in enumerate(coeffs))
    if check and abs(value) < SINGULAR * mp.fsum(abs(c) for c in coeffs):
        return None, True
    slope = mp.fsum(k * c * x**k for k, c in enumerate(coeffs))
    return mp.re(slope / value), False


def delay_at(b: list, a: list, w: mp.mpf) -> tuple[mp.mpf, bool]:
    """(the group delay at w, or its two-sided limit, whether w is singular)"""
    num, num_singular = polynomial_delay(b, w)
    den, den_singular = polynomial_delay(a, w)
    if not (num_singular or den_singular):
        return num - den, False
    either_side = [
        polynomial_delay(b, w + s, check=False)[0]
        - polynomial_delay(a, w + s, check=False)[0]
        for s in (-STEP, STEP)
    ]
    return (either_side[0] + either_side[1]) / 2, True


def exact_group_delay(sections: list, n: int) -> tuple[list, list, list]:
    """
    (w, delay, singular indices) at freqz's n points for the cascade of
    sections, each (b, a): the sum of each section's delay
    """
    with mp.workdps(DPS):
        b_a = [([mp.mpf(c) for c in b], [mp.mpf(c) for c in a]) for b, a in sections]
        ws, delays, singular = [], [], []
        for i in range(n):
            w = mp.pi * i / n
            total, is_singular = mp.mpf(0), False
            for b, a in b_a:
                d, s = delay_at(b, a, w)
                total += d
                is_singular = is_singular or s
            ws.append(float(w))
            delays.append(float(total))
            if is_singular:
                singular.append(i)
        return ws, delays, singular
