"""Filter design: Butterworth coefficients and zpk to second-order sections"""

from typing import Final

import scipy as sc

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
