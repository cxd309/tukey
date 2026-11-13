"""
The exact elliptic filter design, at 50 significant digits with mpmath

scipy.signal.ellip and MATLAB's ellip disagree by up to 1e-8, as MATLAB forms
1 - k1^2 for a small k1 and loses its digits. Neither shows which is right, so
the reference is computed exactly instead, from the same design (Orfanidis,
"Lecture Notes on Elliptic Filter Design", as SciPy's ellipap follows), but
solving each step with mpmath's own functions:

- the degree equation by the nome, with mpmath's qfrom and mfrom
- the inverse sc by the incomplete elliptic integral F(arctan(w) | m)
- K, sn, cn and dn by mpmath's ellipk and ellipfun

Each prototype is checked against its specification, independently of how it
was computed: exactly -ripple at 1 rad/s, within [-ripple, 0] dB in the
passband, exactly -attenuation at the stopband edge 1/sqrt(m), and at least
that far down beyond it.
"""

from typing import Final

import mpmath as mp

DPS: Final = 50  # working precision, in significant digits
SPEC_TOLERANCE_DB: Final = mp.mpf("1e-30")  # how far a prototype may miss its spec


def pow10m1(x: mp.mpf) -> mp.mpf:
    """10^x - 1"""
    return mp.expm1(x * mp.log(10))


def prototype(order: int, rp: mp.mpf, rs: mp.mpf) -> tuple[list, list, mp.mpf, mp.mpf]:
    """
    (zeros, poles, gain, m) of the analog elliptic lowpass prototype: rp dB of
    ripple up to 1 rad/s, at least rs dB down from 1/sqrt(m) rad/s
    """
    eps_sq = pow10m1(rp / 10)
    eps = mp.sqrt(eps_sq)
    if order == 1:
        p = -1 / eps
        return [], [mp.mpc(p)], -p, mp.mpf(0)

    m1 = eps_sq / pow10m1(rs / 10)  # the specification's modulus, squared
    # the degree equation, N·K'(m)/K(m) = K'(m1)/K(m1), through the nome:
    # q(m) = q(m1)^(1/N)
    m = mp.mfrom(q=mp.qfrom(m=m1) ** (mp.mpf(1) / order))
    capk, capk1 = mp.ellipk(m), mp.ellipk(m1)

    js = range(1 - order % 2, order, 2)
    sn = [mp.ellipfun("sn", j * capk / order, m=m) for j in js]
    cn = [mp.ellipfun("cn", j * capk / order, m=m) for j in js]
    dn = [mp.ellipfun("dn", j * capk / order, m=m) for j in js]

    # j = 0, only in an odd order, gives sn = 0: a zero at infinity, and the real
    # pole. mpmath's sn(0) is only 1e-57 or so, so j is tested, not sn.
    zeros = [1j / (mp.sqrt(m) * s) for j, s in zip(js, sn, strict=True) if j != 0]
    zeros += [mp.conj(z) for z in zeros]

    # sc(r | 1 - m1) = 1/eps, so r = F(arctan(1/eps) | 1 - m1)
    r = mp.ellipf(mp.atan(1 / eps), 1 - m1)
    v0 = capk * r / (order * capk1)
    sv = mp.ellipfun("sn", v0, m=1 - m)
    cv = mp.ellipfun("cn", v0, m=1 - m)
    dv = mp.ellipfun("dn", v0, m=1 - m)

    poles = [
        -(c * d * sv * cv + 1j * s * dv) / (1 - (d * sv) ** 2)
        for s, c, d in zip(sn, cn, dn, strict=True)
    ]
    # the real pole of an odd order (j = 0) has no conjugate to add
    complex_poles = [p for j, p in zip(js, poles, strict=True) if j != 0]
    real_poles = [mp.mpc(mp.re(p)) for j, p in zip(js, poles, strict=True) if j == 0]
    poles = real_poles + complex_poles + [mp.conj(p) for p in complex_poles]

    gain = mp.re(mp.fprod(-p for p in poles) / mp.fprod(-z for z in zeros))
    if order % 2 == 0:
        gain /= mp.sqrt(1 + eps_sq)
    return zeros, poles, gain, m


def analog_gain_db(zeros: list, poles: list, gain: mp.mpf, w: mp.mpf) -> mp.mpf:
    """the gain in dB of an analog filter at w rad/s"""
    s = mp.mpc(0, w)
    h = gain * mp.fprod(s - z for z in zeros) / mp.fprod(s - p for p in poles)
    return 20 * mp.log10(abs(h))


def check_prototype(order: int, rp: mp.mpf, rs: mp.mpf) -> None:
    """asserts the prototype meets its specification"""
    zeros, poles, gain, m = prototype(order, rp, rs)
    name = (order, rp, rs)

    def at(w):
        return analog_gain_db(zeros, poles, gain, w)

    assert abs(at(1) + rp) < SPEC_TOLERANCE_DB, name
    for i in range(201):
        g = at(mp.mpf(i) / 200)
        assert -rp - SPEC_TOLERANCE_DB < g < SPEC_TOLERANCE_DB, (name, i)
    if order == 1:  # no stopband ripple: it falls monotonically
        return
    stop = 1 / mp.sqrt(m)
    assert abs(at(stop) + rs) < SPEC_TOLERANCE_DB, name
    for i in range(1, 201):
        assert at(stop * (1 + mp.mpf(i) / 20)) < -rs + SPEC_TOLERANCE_DB, (name, i)


def split_roots(roots: list, wo: mp.mpf) -> list:
    """each root r as the pair r ± sqrt(r² - wo²), all + roots first"""
    disc = [mp.sqrt(r * r - wo * wo) for r in roots]
    return [r + d for r, d in zip(roots, disc, strict=True)] + [
        r - d for r, d in zip(roots, disc, strict=True)
    ]


def transform(zeros, poles, gain, wn, btype):
    """
    the prototype moved onto the band, as scipy.signal's lp2*_zpk, with the
    edges prewarped for the bilinear transform at fs = 2
    """
    degree = len(poles) - len(zeros)

    def warp(w):
        return 4 * mp.tan(mp.pi * mp.mpf(w) / 2)

    def neg_prod_ratio():
        return mp.re(mp.fprod(-z for z in zeros) / mp.fprod(-p for p in poles))

    if btype == "lowpass":
        wo = warp(wn)
        return [wo * z for z in zeros], [wo * p for p in poles], gain * wo**degree
    if btype == "highpass":
        wo = warp(wn)
        new_zeros = [wo / z for z in zeros] + [mp.mpc(0)] * degree
        return new_zeros, [wo / p for p in poles], gain * neg_prod_ratio()
    low, high = warp(wn[0]), warp(wn[1])
    wo, bw = mp.sqrt(low * high), high - low
    if btype == "bandpass":
        new_zeros = split_roots([bw / 2 * z for z in zeros], wo) + [mp.mpc(0)] * degree
        new_poles = split_roots([bw / 2 * p for p in poles], wo)
        return new_zeros, new_poles, gain * bw**degree
    new_zeros = split_roots([bw / 2 / z for z in zeros], wo)
    new_zeros += [mp.mpc(0, wo)] * degree + [mp.mpc(0, -wo)] * degree
    new_poles = split_roots([bw / 2 / p for p in poles], wo)
    return new_zeros, new_poles, gain * neg_prod_ratio()


def bilinear(zeros, poles, gain):
    """the analog filter in the z-plane, s = 4(z - 1)/(z + 1), at fs = 2"""
    degree = len(poles) - len(zeros)
    new_gain = gain * mp.re(
        mp.fprod(4 - z for z in zeros) / mp.fprod(4 - p for p in poles)
    )
    new_zeros = [(4 + z) / (4 - z) for z in zeros] + [mp.mpc(-1)] * degree
    return new_zeros, [(4 + p) / (4 - p) for p in poles], new_gain


def poly(roots: list) -> list:
    """the monic polynomial's coefficients, highest power first"""
    coeffs = [mp.mpc(1)]
    for r in roots:
        coeffs = [*coeffs, mp.mpc(0)]
        for i in range(len(coeffs) - 1, 0, -1):
            coeffs[i] -= r * coeffs[i - 1]
    return coeffs


def exact_ellip(order: int, rp: float, rs: float, wn, btype: str):
    """
    (z, p, k, b, a) of the digital elliptic filter, every value its exact
    design rounded once to float64, after checking the prototype's spec
    """
    with mp.workdps(DPS):
        rp_mp, rs_mp = mp.mpf(rp), mp.mpf(rs)
        check_prototype(order, rp_mp, rs_mp)
        zeros, poles, gain, _ = prototype(order, rp_mp, rs_mp)
        zeros, poles, gain = bilinear(*transform(zeros, poles, gain, wn, btype))
        b = [mp.re(gain * c) for c in poly(zeros)]
        a = [mp.re(c) for c in poly(poles)]
        return (
            [complex(z) for z in zeros],
            [complex(p) for p in poles],
            float(gain),
            [float(c) for c in b],
            [float(c) for c in a],
        )
