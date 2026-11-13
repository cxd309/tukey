"""
The exact Bessel filter design, at 50 significant digits with mpmath

The analog Bessel lowpass is H(s) = θ(0)/θ(s), for θ the reverse Bessel
polynomial of order N, whose coefficient of s^k is
(2N - k)! / (2^(N - k)·k!·(N - k)!). Its poles are θ's roots, found here by
mpmath's polyroots, and give a group delay at DC of exactly 1 s. That is one of
three normalisations, as scipy.signal.besselap's norm:

- "delay": as is, the group delay at DC 1 s
- "phase": the poles scaled to multiply to magnitude 1, so the response falls
  as 1/w^N at high frequency, as Butterworth's does from the same cutoff, and
  the phase is near -N·pi/4 at 1 rad/s. MATLAB's besselap.
- "mag": the poles scaled so the gain is -3 dB at 1 rad/s, as Butterworth's

Each prototype is checked against its normalisation, independently of how it
was found, then moved onto the band and through the bilinear transform with
the elliptic design's exact steps.
"""

from typing import Final

import mpmath as mp

from dsp.elliptic import DPS, bilinear, poly, transform

NORMS: Final = ("phase", "delay", "mag")
# how far a prototype may miss its normalisation
SPEC_TOLERANCE: Final = mp.mpf("1e-35")


def reverse_bessel(order: int) -> list:
    """θ's coefficients, highest power first"""
    n = order
    return [
        mp.factorial(2 * n - k) / (2 ** (n - k) * mp.factorial(k) * mp.factorial(n - k))
        for k in range(n, -1, -1)
    ]


def gain_sq(poles: list, gain: mp.mpf, w: mp.mpf) -> mp.mpf:
    """|H(jw)|² of an all-pole analog filter"""
    s = mp.mpc(0, w)
    return abs(gain / mp.fprod(s - p for p in poles)) ** 2


def prototype(order: int, norm: str) -> tuple[list, mp.mpf]:
    """(poles, gain) of the analog Bessel lowpass prototype, normalised by norm"""
    coeffs = reverse_bessel(order)
    poles = mp.polyroots(coeffs, maxsteps=500, extraprec=4 * DPS)
    poles = [mp.mpc(p) for p in poles]
    dc = coeffs[-1]  # θ(0)

    if norm == "delay":
        return poles, dc
    if norm == "phase":
        scale = dc ** (mp.mpf(-1) / order)
        return [scale * p for p in poles], mp.mpf(1)
    # "mag": the -3 dB frequency of the delay-normalised filter, where the gain
    # falls monotonically from 1, bracketed between 0 and 2N + 1
    cutoff = mp.findroot(
        lambda w: gain_sq(poles, dc, w) - mp.mpf(1) / 2,
        (mp.mpf("1e-3"), mp.mpf(2 * order + 1)),
        solver="anderson",
    )
    return [p / cutoff for p in poles], dc / cutoff**order


def check_prototype(order: int, norm: str) -> None:
    """asserts the prototype is a Bessel filter, normalised as norm says"""
    poles, gain = prototype(order, norm)
    name = (order, norm)
    assert len(poles) == order, name
    assert all(mp.re(p) < 0 for p in poles), name

    # θ, scaled by the normalisation, has the poles as its roots
    scale = mp.fprod(-p for p in poles) / reverse_bessel(order)[-1]
    rebuilt = poly(poles)
    for k, c in enumerate(reverse_bessel(order)):
        want = c * (scale ** (mp.mpf(k) / order))
        assert abs(rebuilt[k] - want) < SPEC_TOLERANCE * abs(want), (name, k)

    # 0 dB at DC
    assert abs(gain_sq(poles, gain, mp.mpf(0)) - 1) < SPEC_TOLERANCE, name
    if norm == "delay":  # the group delay at DC, Σ Re(-1/p)
        assert abs(mp.fsum(mp.re(-1 / p) for p in poles) - 1) < SPEC_TOLERANCE, name
    elif norm == "phase":
        assert abs(mp.fprod(abs(p) for p in poles) - 1) < SPEC_TOLERANCE, name
        assert abs(gain - 1) < SPEC_TOLERANCE, name
    else:
        half_power = gain_sq(poles, gain, mp.mpf(1)) - mp.mpf(1) / 2
        assert abs(half_power) < SPEC_TOLERANCE, name


def exact_bessel(order: int, wn, btype: str, norm: str):
    """
    (z, p, k, b, a) of the digital Bessel filter, every value its exact design
    rounded once to float64, after checking the prototype's normalisation
    """
    with mp.workdps(DPS):
        check_prototype(order, norm)
        poles, gain = prototype(order, norm)
        zeros, poles, gain = bilinear(*transform([], poles, gain, wn, btype))
        b = [mp.re(gain * c) for c in poly(zeros)]
        a = [mp.re(c) for c in poly(poles)]
        return (
            [complex(z) for z in zeros],
            [complex(p) for p in poles],
            float(gain),
            [float(c) for c in b],
            [float(c) for c in a],
        )
