"""
The exact second-order notch and peak designs, at 50 significant digits with
mpmath

scipy.signal.iirnotch and iirpeak in closed form, for a centre w0 normalised to
Nyquist and quality factor Q, whose -3 dB bandwidth is w0/Q:

- beta = tan(pi·w0/(2Q)) and g = 1/(1 + beta)
- denominator z² - 2g·cos(pi·w0)·z + (1 - beta)/(1 + beta), for both
- notch: g·(z² - 2·cos(pi·w0)·z + 1), zeros at e^(±j·pi·w0)
- peak: (1 - g)·(z² - 1), zeros at ±1

Each design is checked against its specification, independently of the
formula: a notch is 0 dB at DC and Nyquist and has a null at w0, a peak is
0 dB at w0 and nulls at DC and Nyquist, and the two -3 dB edges, found from
the response itself, are exactly w0/Q apart (Orfanidis, "Introduction to
Signal Processing", section 11.3: beta = tan(width/2) for the digital width).
"""

from typing import Final

import mpmath as mp

from dsp.elliptic import DPS, poly

SPEC_TOLERANCE: Final = mp.mpf("1e-35")  # how far a design may miss its spec


def response(b: list, a: list, w: mp.mpf) -> mp.mpc:
    """H at w, normalised to Nyquist, of real coefficients, highest power first"""
    z = mp.expjpi(w)
    return mp.polyval(b, z) / mp.polyval(a, z)


def design(w0: mp.mpf, q: mp.mpf, ftype: str) -> tuple[list, list, list, mp.mpf]:
    """(zeros, poles, gain, a) of the notch or peak, a highest power first"""
    beta = mp.tan(mp.pi * w0 / (2 * q))
    g = 1 / (1 + beta)
    c = mp.cospi(w0)
    a = [mp.mpf(1), -2 * g * c, (1 - beta) / (1 + beta)]
    disc = mp.sqrt(mp.mpc((g * c) ** 2 - a[2]))
    poles = [g * c + disc, g * c - disc]
    if ftype == "notch":
        zeros = [mp.expjpi(w0), mp.expjpi(-w0)]
        return zeros, poles, g, a
    return [mp.mpc(1), mp.mpc(-1)], poles, beta / (1 + beta), a


def half_power_edges(b: list, a: list, w0: mp.mpf, ftype: str) -> list:
    """the two frequencies either side of w0 where |H|² is 1/2"""

    def excess(w):
        return abs(response(b, a, w)) ** 2 - mp.mpf(1) / 2

    # |H| is monotonic on each side of w0, between w0 and an end
    return [
        mp.findroot(excess, (mp.mpf(0), w0), solver="anderson"),
        mp.findroot(excess, (w0, mp.mpf(1)), solver="anderson"),
    ]


def check_design(w0: mp.mpf, q: mp.mpf, ftype: str) -> None:
    """asserts the design meets its specification"""
    zeros, poles, gain, a = design(w0, q, ftype)
    b = [gain * c for c in poly(zeros)]
    name = (w0, q, ftype)
    assert all(abs(p) < 1 for p in poles), name
    for c in b + a:
        assert abs(mp.im(c)) < SPEC_TOLERANCE, name
    b, a = [mp.re(c) for c in b], [mp.re(c) for c in a]

    passes, stops = ([0, 1], [w0]) if ftype == "notch" else ([w0], [0, 1])
    for w in passes:
        assert abs(abs(response(b, a, mp.mpf(w))) - 1) < SPEC_TOLERANCE, (name, w)
    for w in stops:
        assert abs(response(b, a, mp.mpf(w))) < SPEC_TOLERANCE, (name, w)

    low, high = half_power_edges(b, a, w0, ftype)
    assert low < w0 < high, name
    assert abs((high - low) - w0 / q) < SPEC_TOLERANCE, name


def exact_notch(w0: float, q: float, ftype: str):
    """
    (z, p, k, b, a) of the notch or peak, every value exact rounded once to
    float64, after checking it meets its specification
    """
    with mp.workdps(DPS):
        w0_mp, q_mp = mp.mpf(w0), mp.mpf(q)
        check_design(w0_mp, q_mp, ftype)
        zeros, poles, gain, a = design(w0_mp, q_mp, ftype)
        b = [mp.re(gain * c) for c in poly(zeros)]
        return (
            [complex(z) for z in zeros],
            [complex(p) for p in poles],
            float(gain),
            [float(c) for c in b],
            [float(c) for c in a],
        )
