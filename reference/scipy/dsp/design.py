"""
Filter design: Butterworth, Chebyshev, elliptic and Bessel coefficients, order
selection, and zpk to second-order sections
"""

from collections.abc import Callable
from typing import Final

import numpy as np
import scipy as sc

from dsp.bessel import NORMS, exact_bessel
from dsp.elliptic import exact_ellip
from dsp.notch import exact_notch
from dsp.order import band_type, exact_order
from write import complex_json, write_vector

# (order, wn, btype) for scipy.signal.butter with a single cutoff
BUTTER_LOWHIGH_CASES: Final = [
    (1, 0.3, "lowpass"),  # simplest: one real pole, sanity check
    (2, 0.3, "lowpass"),  # first conjugate pair
    (4, 0.3, "lowpass"),  # multiple conjugate pairs
    (2, 0.5, "lowpass"),  # prewarp == fs2; a[1] is exactly 0
    (1, 0.3, "highpass"),  # zero at origin maps to z=+1
    (2, 0.3, "highpass"),  # conjugate pair through lp2hp gain
    (4, 0.3, "highpass"),  # multiple pairs through lp2hp
    (2, 0.5, "highpass"),  # mirror of lowpass 0.5; a[1] is 0
    (8, 0.01, "lowpass"),  # tiny b (~1e-14); guards tolerance
    (3, 0.3, "lowpass"),  # odd order: real pole at s=-1
    (5, 0.7, "highpass"),  # odd order; cutoff above 0.5
    (10, 0.3, "lowpass"),  # high order: poly conditioning
    (2, 0.99, "lowpass"),  # near Nyquist: steep tan in prewarp
    (2, 0.01, "highpass"),  # poles right next to zeros at z=+1
    (2, 0.001, "lowpass"),  # target application requirement
    (2, 0.003, "lowpass"),  # target application requirement
    (2, 0.0008, "highpass"),  # target application requirement
]

# (order, (low, high), btype) for scipy.signal.butter with two band edges
BUTTER_BAND_CASES: Final = [
    (2, (0.2, 0.5), "bandpass"),  # basic lp2bp: pole count doubles
    (4, (0.2, 0.5), "bandpass"),  # multiple pairs through lp2bp
    (2, (0.2, 0.5), "bandstop"),  # notch zeros land on unit circle
    (4, (0.2, 0.5), "bandstop"),  # multiple pairs through lp2bs
    (1, (0.2, 0.5), "bandpass"),  # sqrt of negative real: branch cut
    (3, (0.2, 0.5), "bandpass"),  # odd order: real + pair quadratics
    (3, (0.2, 0.5), "bandstop"),  # odd order: lp2bs divides by real p
    (2, (0.30, 0.31), "bandpass"),  # narrow band: clustered poles
    (2, (0.01, 0.99), "bandstop"),  # edges ~4000x apart; worst case
    (8, (0.2, 0.5), "bandpass"),  # 16th-order poly: stress test
    (2, (0.0032, 0.032), "bandpass"),  # target application requirement
    (1, (0.55, 0.7), "bandpass"),  # target application requirement
]

# (order, passband ripple dB, stopband attenuation dB, wn, btype), for
# scipy.signal.cheby1, which uses the ripple, and cheby2, the attenuation
CHEBY_CASES: Final = [
    (1, 1.0, 40.0, 0.3, "lowpass"),  # simplest: one real pole
    (2, 1.0, 40.0, 0.3, "lowpass"),  # even order: cheby1's DC gain is -ripple
    (3, 0.5, 30.0, 0.3, "lowpass"),  # odd: cheby1 0 dB at DC, cheby2 zero at inf
    (4, 1.0, 40.0, 0.3, "lowpass"),  # multiple conjugate pairs
    (4, 0.01, 80.0, 0.3, "lowpass"),  # tiny ripple, deep stopband
    (4, 6.0, 10.0, 0.3, "lowpass"),  # large ripple, shallow stopband
    (2, 1.0, 40.0, 0.5, "lowpass"),  # prewarp == fs2
    (10, 1.0, 40.0, 0.3, "lowpass"),  # high order: poly conditioning
    (8, 0.5, 60.0, 0.01, "lowpass"),  # poles packed against z = 1
    (2, 1.0, 40.0, 0.99, "lowpass"),  # near Nyquist: steep tan in prewarp
    (2, 1.0, 40.0, 0.001, "lowpass"),  # very low cutoff
    (4, 1.0, 40.0, 0.3, "highpass"),  # through lp2hp
    (5, 1.0, 40.0, 0.7, "highpass"),  # odd order; cutoff above 0.5
    (2, 1.0, 40.0, 0.0008, "highpass"),  # very low cutoff
    (1, 1.0, 40.0, (0.2, 0.5), "bandpass"),  # one prototype pole
    (2, 1.0, 40.0, (0.2, 0.5), "bandpass"),  # pole count doubles
    (3, 1.0, 40.0, (0.2, 0.5), "bandpass"),  # odd order
    (8, 1.0, 40.0, (0.2, 0.5), "bandpass"),  # 16th-order poly: stress test
    (2, 1.0, 40.0, (0.30, 0.31), "bandpass"),  # narrow band: clustered poles
    (2, 1.0, 40.0, (0.2, 0.5), "bandstop"),  # zeros on the unit circle
    (3, 1.0, 40.0, (0.2, 0.5), "bandstop"),  # odd order
    (4, 1.0, 40.0, (0.2, 0.5), "bandstop"),  # multiple pairs through lp2bs
    (2, 1.0, 40.0, (0.01, 0.99), "bandstop"),  # edges far apart
]

# CHEBY_CASES, then the sharp transitions where the elliptic maths is most
# delicate, for scipy.signal.ellip, which uses both the ripple and attenuation
ELLIP_CASES: Final = CHEBY_CASES + [
    (12, 0.5, 120.0, 0.3, "lowpass"),  # high order, very deep stopband
    (8, 0.01, 120.0, 0.2, "lowpass"),  # tiny ripple and deep stopband: modulus near 1
    (3, 3.0, 6.0, 0.3, "lowpass"),  # attenuation barely above ripple
    (6, 1.0, 60.0, (0.2, 0.25), "bandpass"),  # narrow band at high order
]

# (order, wn, btype) for scipy.signal.bessel, each written in every norm
BESSEL_CASES: Final = [
    (1, 0.3, "lowpass"),  # simplest: one real pole
    (2, 0.3, "lowpass"),  # first conjugate pair
    (3, 0.3, "lowpass"),  # odd order: real pole and a pair
    (4, 0.3, "lowpass"),  # multiple conjugate pairs
    (8, 0.3, "lowpass"),  # higher order
    (25, 0.3, "lowpass"),  # MATLAB's besselap's highest order
    (30, 0.3, "lowpass"),  # beyond it: SciPy only
    (10, 0.01, "lowpass"),  # poles packed against z = 1
    (2, 0.99, "lowpass"),  # near Nyquist: steep tan in prewarp
    (2, 0.001, "lowpass"),  # very low cutoff
    (4, 0.3, "highpass"),  # through lp2hp
    (5, 0.7, "highpass"),  # odd order; cutoff above 0.5
    (2, 0.0008, "highpass"),  # very low cutoff
    (1, (0.2, 0.5), "bandpass"),  # one prototype pole
    (4, (0.2, 0.5), "bandpass"),  # pole count doubles
    (2, (0.30, 0.31), "bandpass"),  # narrow band: clustered poles
    (3, (0.2, 0.5), "bandstop"),  # odd order through lp2bs
    (2, (0.01, 0.99), "bandstop"),  # edges far apart
]

# (w0, Q) for scipy.signal.iirnotch and iirpeak, each written for both
NOTCH_CASES: Final = [
    (0.3, 10.0),  # the documentation's shape
    (0.05, 35.0),  # mains hum: 50 Hz at 2 kHz
    (0.5, 30.0),  # cos(pi·w0) = 0: a[1] is exactly 0
    (0.1, 1.0),  # wide
    (0.3, 1000.0),  # very narrow: poles hard against the unit circle
    (0.002, 30.0),  # very low centre: poles near z = 1
    (0.95, 5.0),  # near Nyquist
    (0.25, 0.5),  # beta = 1: a[2] is exactly 0, the poles at ±j/sqrt(2)
    (0.6, 0.8),  # bandwidth past half Nyquist: beta > 1, real poles
]

# where each dB parameter sits in a CHEBY_CASES or ELLIP_CASES row
DB_COLUMNS: Final = {"ripple_db": 1, "attenuation_db": 2}

# (wp, ws, ripple dB, attenuation dB) for the order selections; the band type
# follows from the edges, as in scipy.signal's *ord
ORD_CASES: Final = [
    (0.2, 0.3, 1.0, 40.0),  # lowpass: the documentation's shape
    (0.3, 0.32, 0.5, 60.0),  # narrow transition: high order
    (0.1, 0.5, 3.0, 6.0),  # wide transition, shallow stopband: order 1 or 2
    (0.01, 0.015, 0.01, 80.0),  # low edges, tiny ripple, deep stopband
    (0.9, 0.95, 1.0, 40.0),  # near Nyquist: steep tan in the prewarp
    (0.3, 0.2, 1.0, 40.0),  # highpass
    (0.7, 0.6, 0.1, 60.0),  # highpass above 0.5
    (0.002, 0.001, 1.0, 30.0),  # highpass with very low edges
    ((0.2, 0.5), (0.1, 0.6), 1.0, 40.0),  # bandpass
    ((0.2, 0.3), (0.1, 0.5), 1.0, 40.0),  # bandpass, the upper stopband edge slack
    ((0.3, 0.31), (0.28, 0.33), 0.5, 50.0),  # narrow bandpass
    ((0.1, 0.6), (0.2, 0.5), 1.0, 40.0),  # bandstop: the lower passband edge moves
    ((0.1, 0.8), (0.2, 0.3), 1.0, 40.0),  # bandstop: the upper passband edge moves
    ((0.05, 0.9), (0.4, 0.5), 0.5, 60.0),  # bandstop, wide transitions
    ((0.3, 0.5), (0.35, 0.45), 3.0, 30.0),  # bandstop, narrow
]

# name -> (z, p, k) for scipy.signal.zpk2sos, beyond what Butterworth designs produce
ZPK2SOS_CASES = {
    "empty": ([], [], 2.0),  # no roots: a single gain-only section
    "one_real_pole": ([], [0.5], 1.0),  # first-order: padded with zeros at the origin
    "scipy_doc_mixed_odd": (  # SciPy's docstring example: odd order, real + complex
        [-1, -0.5 - 0.5j, -0.5 + 0.5j],
        [0.75, 0.8 + 0.1j, 0.8 - 0.1j],
        1.0,
    ),
    "complex_zero_real_poles": (  # pole real, nearest zero complex: pairing step 3.3
        [0.3 + 0.9j, 0.3 - 0.9j],
        [0.5, -0.4],
        1.0,
    ),
    "more_zeros_than_poles": (
        [1j, -1j, -1.0],
        [0.5],
        3.0,
    ),  # poles padded at the origin
}


def generate_butter() -> None:
    """
    b, a coefficients from scipy.signal.butter, analog=False.

    wn is normalised to Nyquist (1.0 == Nyquist), matching both scipy and
    MATLAB. Low/high-pass take a scalar wn; band-pass/stop take [low, high].
    """

    for order, wn, btype in BUTTER_LOWHIGH_CASES:
        b, a = sc.signal.butter(order, wn, btype=btype, analog=False)
        name = f"{btype}_order{order}_wn{wn}"
        write_vector(
            "butter",
            name,
            f"scipy.signal.butter({order}, {wn}, btype='{btype}')",
            {"order": order, "wn": wn, "btype": btype},
            {
                "b": b.tolist(),
                "a": a.tolist(),
                "sos": sc.signal.butter(order, wn, btype=btype, output="sos").tolist(),
            },
        )

    for order, wn, btype in BUTTER_BAND_CASES:
        b, a = sc.signal.butter(order, wn, btype=btype, analog=False)
        name = f"{btype}_order{order}_wn{wn[0]}-{wn[1]}"
        write_vector(
            "butter",
            name,
            f"scipy.signal.butter({order}, {list(wn)}, btype='{btype}')",
            {"order": order, "wn": list(wn), "btype": btype},
            {
                "b": b.tolist(),
                "a": a.tolist(),
                "sos": sc.signal.butter(order, wn, btype=btype, output="sos").tolist(),
            },
        )


def write_iir(
    category: str, design: Callable, cases: list, db_params: tuple[str, ...]
) -> None:
    """
    b, a and sos from a design taking dB arguments, one vector per case row,
    passing the row's db_params (see DB_COLUMNS) in order, by those names
    """
    for case in cases:
        order, wn, btype = case[0], case[3], case[4]
        dbs = [case[DB_COLUMNS[param]] for param in db_params]
        db_name = f"db{dbs[0]}" if len(dbs) == 1 else f"rp{dbs[0]}_rs{dbs[1]}"
        wn_name = f"{wn[0]}-{wn[1]}" if isinstance(wn, tuple) else wn
        wn_param = list(wn) if isinstance(wn, tuple) else wn
        db_args = ", ".join(str(db) for db in dbs)
        b, a = design(order, *dbs, wn, btype=btype)
        write_vector(
            category,
            f"{btype}_order{order}_{db_name}_wn{wn_name}",
            f"scipy.signal.{category}({order}, {db_args}, {wn_param}, btype='{btype}')",
            {
                "order": order,
                **dict(zip(db_params, dbs, strict=True)),
                "wn": wn_param,
                "btype": btype,
            },
            {
                "b": b.tolist(),
                "a": a.tolist(),
                "sos": design(order, *dbs, wn, btype=btype, output="sos").tolist(),
            },
        )


def generate_cheby1() -> None:
    """b, a and sos from scipy.signal.cheby1: equal ripple in the passband"""
    write_iir("cheby1", sc.signal.cheby1, CHEBY_CASES, ("ripple_db",))


def generate_cheby2() -> None:
    """
    b, a and sos from scipy.signal.cheby2: equal ripple in the stopband, whose
    start is wn
    """
    write_iir("cheby2", sc.signal.cheby2, CHEBY_CASES, ("attenuation_db",))


def generate_ellip() -> None:
    """
    The exact elliptic design (see dsp.elliptic): equal ripple in both bands,
    the passband's edge at wn as for cheby1. b and a are exact, rounded once;
    sos is scipy.signal.zpk2sos of the exact zeros, poles and gain, each rounded
    once. Asserts scipy.signal.ellip agrees to within its own accuracy.
    """
    for order, rp, rs, wn, btype in ELLIP_CASES:
        z, p, k, b, a = exact_ellip(order, rp, rs, wn, btype)
        scipy_b, scipy_a = sc.signal.ellip(order, rp, rs, wn, btype=btype)
        for exact, got in ((b, scipy_b), (a, scipy_a)):
            scale = np.max(np.abs(exact))
            case = (order, rp, rs, wn)
            assert np.allclose(got, exact, rtol=0, atol=1e-11 * scale), case
        wn_name = f"{wn[0]}-{wn[1]}" if isinstance(wn, tuple) else wn
        wn_param = list(wn) if isinstance(wn, tuple) else wn
        write_vector(
            "ellip",
            f"{btype}_order{order}_rp{rp}_rs{rs}_wn{wn_name}",
            f"the exact elliptic design, as scipy.signal.ellip({order}, {rp}, {rs}, "
            f"{wn_param}, btype='{btype}')",
            {
                "order": order,
                "ripple_db": rp,
                "attenuation_db": rs,
                "wn": wn_param,
                "btype": btype,
            },
            {"b": b, "a": a, "sos": sc.signal.zpk2sos(z, p, k).tolist()},
        )


def generate_bessel() -> None:
    """
    The exact Bessel design (see dsp.bessel) in each of scipy.signal.bessel's
    normalisations. b and a are exact, rounded once; sos is
    scipy.signal.zpk2sos of the exact zeros, poles and gain, each rounded once.
    Asserts scipy.signal.bessel agrees to within its own accuracy.
    """
    for order, wn, btype in BESSEL_CASES:
        for norm in NORMS:
            z, p, k, b, a = exact_bessel(order, wn, btype, norm)
            scipy_b, scipy_a = sc.signal.bessel(order, wn, btype=btype, norm=norm)
            for exact, got in ((b, scipy_b), (a, scipy_a)):
                scale = np.max(np.abs(exact))
                case = (order, wn, btype, norm)
                assert np.allclose(got, exact, rtol=0, atol=1e-11 * scale), case
            wn_name = f"{wn[0]}-{wn[1]}" if isinstance(wn, tuple) else wn
            wn_param = list(wn) if isinstance(wn, tuple) else wn
            write_vector(
                "bessel",
                f"{btype}_order{order}_{norm}_wn{wn_name}",
                f"the exact Bessel design, as scipy.signal.bessel({order}, "
                f"{wn_param}, btype='{btype}', norm='{norm}')",
                {"order": order, "wn": wn_param, "btype": btype, "norm": norm},
                {"b": b, "a": a, "sos": sc.signal.zpk2sos(z, p, k).tolist()},
            )


# how far scipy.signal's notch and peak may be from the exact design, relative
# to the largest coefficient: its peak's gain is 1 - g, which cancels for a
# narrow bandwidth, off by 6.1e-13 at w0 = 0.002, Q = 30, where beta/(1 + beta)
# keeps every digit
NOTCH_SCIPY_ACCURACY: Final = {"notch": 1e-15, "peak": 1e-12}


def write_notch(category: str, ftype: str, scipy_design: Callable) -> None:
    """
    The exact second-order notch or peak (see dsp.notch): b, a and the zeros,
    poles and gain, each exact rounded once. Asserts scipy.signal agrees to
    within its own accuracy.
    """
    for w0, q in NOTCH_CASES:
        z, p, k, b, a = exact_notch(w0, q, ftype)
        scipy_b, scipy_a = scipy_design(w0, q)
        for exact, got in ((b, scipy_b), (a, scipy_a)):
            atol = NOTCH_SCIPY_ACCURACY[ftype] * np.max(np.abs(exact))
            assert np.allclose(got, exact, rtol=0, atol=atol), (w0, q)
        write_vector(
            category,
            f"w0{w0}_q{q}",
            f"the exact design, as scipy.signal.{category}({w0}, {q})",
            {"w0": w0, "q": q},
            {"b": b, "a": a, "z": complex_json(z), "p": complex_json(p), "k": k},
        )


def generate_iirnotch() -> None:
    """a notch at w0, -3 dB at w0 ± w0/(2Q)"""
    write_notch("iirnotch", "notch", sc.signal.iirnotch)


def generate_iirpeak() -> None:
    """a peak at w0, -3 dB at w0 ± w0/(2Q)"""
    write_notch("iirpeak", "peak", sc.signal.iirpeak)


def write_ord(category: str, kind: str, scipy_ord: Callable) -> None:
    """
    The exact order and both natural frequencies (see dsp.order), one vector
    per ORD_CASES row. Asserts scipy_ord agrees on the order, and on the
    passband-matched wn to within its own accuracy: about 1e-5 for bandstop,
    whose edges SciPy finds by a numerical minimiser.
    """
    for wp, ws, rp, rs in ORD_CASES:
        order, wn_passband, wn_stopband = exact_order(kind, wp, ws, rp, rs)
        scipy_order, scipy_wn = scipy_ord(wp, ws, rp, rs)
        case = (category, wp, ws, rp, rs)
        assert scipy_order == order, (case, scipy_order, order)
        tolerance = 1e-5 if band_type(wp, ws) == "bandstop" else 1e-12
        assert np.allclose(scipy_wn, wn_passband, rtol=0, atol=tolerance), case

        def edges(w):
            return f"{w[0]}-{w[1]}" if isinstance(w, tuple) else w

        def param(w):
            return list(w) if isinstance(w, tuple) else w

        write_vector(
            category,
            f"{band_type(wp, ws)}_wp{edges(wp)}_ws{edges(ws)}_rp{rp}_rs{rs}",
            f"the exact order selection, as scipy.signal.{category}({param(wp)}, "
            f"{param(ws)}, {rp}, {rs}), with wn matching each band",
            {
                "wp": param(wp),
                "ws": param(ws),
                "ripple_db": rp,
                "attenuation_db": rs,
            },
            {
                "order": order,
                "wn_passband": wn_passband,
                "wn_stopband": wn_stopband,
            },
        )


def generate_buttord() -> None:
    """Butterworth order selection: Wn is the -3 dB point"""
    write_ord("buttord", "butter", sc.signal.buttord)


def generate_cheb1ord() -> None:
    """Chebyshev type I order selection: Wn is the passband edge"""
    write_ord("cheb1ord", "cheby1", sc.signal.cheb1ord)


def generate_cheb2ord() -> None:
    """Chebyshev type II order selection: Wn is the stopband edge"""
    write_ord("cheb2ord", "cheby2", sc.signal.cheb2ord)


def generate_ellipord() -> None:
    """Elliptic order selection: Wn is the passband edge"""
    write_ord("ellipord", "ellip", sc.signal.ellipord)


def generate_zpk2sos() -> None:
    """sos from scipy.signal.zpk2sos: SciPy's own designs plus hand-built cases."""
    cases = dict(ZPK2SOS_CASES)
    for order, wn, btype in BUTTER_LOWHIGH_CASES + BUTTER_BAND_CASES:
        z, p, k = sc.signal.butter(order, wn, btype=btype, output="zpk")
        wn_name = f"{wn[0]}-{wn[1]}" if isinstance(wn, tuple) else wn
        cases[f"butter_{btype}_order{order}_wn{wn_name}"] = (z, p, k)

    for name, (z, p, k) in cases.items():
        sos = sc.signal.zpk2sos(z, p, k)
        write_vector(
            "zpk2sos",
            name,
            f"scipy.signal.zpk2sos(z, p, k) for {name}",
            {"z": complex_json(z), "p": complex_json(p), "k": float(k)},
            {"sos": sos.tolist()},
        )
