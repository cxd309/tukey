"""
Exact filter order selection, at 50 significant digits with mpmath

For a specification (passband edge wp, stopband edge ws, rp dB of passband
ripple, at least rs dB of stopband attenuation), each design's order is the
smallest whole number meeting it. That order usually has margin to spare, and
the natural frequency Wn decides where it goes:

- "passband": the passband edge is met exactly, at -rp dB, and the stopband
  gets the margin. scipy.signal's buttord, cheb1ord, cheb2ord and ellipord.
- "stopband": the stopband edge is met exactly, at -rs dB, and the passband
  extends past wp. MATLAB's buttord and cheb2ord; neither tool has it for
  Chebyshev type I or elliptic designs.

Each is found as rho, the prototype's Wn in units where the passband edge
is 1 rad/s and the stopband edge nat, then moved onto the band.

A bandstop specification's passband edges are first moved inward until
P1·P2 = S1·S2, where both stopband edges bind equally, which is the lowest
order. SciPy and MATLAB find the same edges by a numerical minimiser, to about
1e-5; here they're exact.
"""

from typing import Final

import mpmath as mp

DPS: Final = 50  # working precision, in significant digits
INTEGER_MARGIN: Final = mp.mpf("1e-9")  # how close a raw order may be to a whole number


def pow10m1(x: mp.mpf) -> mp.mpf:
    """10^x - 1"""
    return mp.expm1(x * mp.log(10))


def band_type(wp, ws) -> str:
    """the band type wp and ws describe, as scipy.signal's *ord"""
    if isinstance(wp, (list, tuple)):
        return "bandstop" if wp[0] < ws[0] else "bandpass"
    return "lowpass" if wp < ws else "highpass"


def prewarp(w):
    """digital edges, normalised to Nyquist, as analog ones"""
    if isinstance(w, (list, tuple)):
        return [mp.tan(mp.pi * mp.mpf(x) / 2) for x in w]
    return mp.tan(mp.pi * mp.mpf(w) / 2)


def unwarp(w: mp.mpf) -> mp.mpf:
    """an analog edge as a digital one, normalised to Nyquist"""
    return 2 * mp.atan(w) / mp.pi


def stopband_ratio(btype: str, p, s) -> mp.mpf:
    """
    nat: the stopband edge of the equivalent lowpass, whose passband edge is 1,
    taking the tighter of two stopband edges
    """
    if btype == "lowpass":
        return s / p
    if btype == "highpass":
        return p / s
    if btype == "bandpass":
        return min(abs((x * x - p[0] * p[1]) / (x * (p[1] - p[0]))) for x in s)
    return min(abs(x * (p[1] - p[0]) / (p[0] * p[1] - x * x)) for x in s)


def symmetric_bandstop(p: list, s: list) -> list:
    """
    the bandstop passband edges moved inward until P1·P2 = S1·S2: the edge on
    the side with margin moves, the other stays
    """
    product = s[0] * s[1]
    if p[0] * p[1] < product:
        return [product / p[1], p[1]]
    return [p[0], product / p[0]]


def raw_order(kind: str, nat: mp.mpf, gpass: mp.mpf, gstop: mp.mpf) -> mp.mpf:
    """
    the real-valued order that exactly meets the spec, where gpass and gstop
    are 10^(dB/10) - 1
    """
    if kind == "butter":
        return mp.log(gstop / gpass) / (2 * mp.log(nat))
    if kind in ("cheby1", "cheby2"):
        return mp.acosh(mp.sqrt(gstop / gpass)) / mp.acosh(nat)
    m, m1 = 1 / nat**2, gpass / gstop
    return mp.ellipk(m) * mp.ellipk(1 - m1) / (mp.ellipk(1 - m) * mp.ellipk(m1))


def prototype_wn(kind, match, order, nat, gpass, gstop) -> mp.mpf:
    """
    rho: the order's prototype Wn, in units where the passband edge is 1 and
    the stopband edge nat, putting the margin in the other band
    """
    n = mp.mpf(order)
    selectivity = mp.cosh(mp.acosh(mp.sqrt(gstop / gpass)) / n)
    if kind == "butter":  # Wn is the -3 dB point
        if match == "passband":
            return gpass ** (-1 / (2 * n))
        return nat * gstop ** (-1 / (2 * n))
    if kind == "cheby1":  # Wn is the passband edge
        return mp.mpf(1) if match == "passband" else nat / selectivity
    if kind == "cheby2":  # Wn is the stopband edge
        return selectivity if match == "passband" else nat
    # ellip: Wn is the passband edge, and the stopband edge Wn/sqrt(m), with m
    # the order's own modulus squared, from the degree equation
    if match == "passband":
        return mp.mpf(1)
    m = mp.mfrom(q=mp.qfrom(m=gpass / gstop) ** (1 / n))
    return mp.sqrt(m) * nat


def band_edges(btype: str, p, rho: mp.mpf) -> list:
    """the analog edges where the equivalent lowpass is at rho"""
    if btype == "lowpass":
        return [rho * p]
    if btype == "highpass":
        return [p / rho]
    width, centre_sq = p[1] - p[0], p[0] * p[1]
    if btype == "bandpass":
        root = mp.sqrt(rho**2 * width**2 / 4 + centre_sq)
        return [root - rho * width / 2, root + rho * width / 2]
    root = mp.sqrt(width**2 + 4 * rho**2 * centre_sq)
    return [(root - width) / (2 * rho), (root + width) / (2 * rho)]


def exact_order(kind: str, wp, ws, rp: float, rs: float) -> tuple[int, list, list]:
    """
    (order, wn matching the passband, wn matching the stopband), each wn as a
    list of floats, after checking the raw order isn't within INTEGER_MARGIN of
    a whole number, where rounding could decide it
    """
    with mp.workdps(DPS):
        btype = band_type(wp, ws)
        p, s = prewarp(wp), prewarp(ws)
        if btype == "bandstop":
            p = symmetric_bandstop(p, s)
        nat = stopband_ratio(btype, p, s)
        gpass, gstop = pow10m1(mp.mpf(rp) / 10), pow10m1(mp.mpf(rs) / 10)

        raw = raw_order(kind, nat, gpass, gstop)
        order = int(mp.ceil(raw))
        case = (kind, wp, ws, rp, rs)
        assert abs(raw - mp.nint(raw)) > INTEGER_MARGIN, (case, raw)

        wns = []
        for match in ("passband", "stopband"):
            rho = prototype_wn(kind, match, order, nat, gpass, gstop)
            wns.append([float(unwarp(w)) for w in band_edges(btype, p, rho)])
        return order, wns[0], wns[1]
