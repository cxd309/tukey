"""Filter design: Butterworth coefficients and zpk to second-order sections"""

import scipy as sc

from .cases import BUTTER_BAND_CASES, BUTTER_LOWHIGH_CASES, ZPK2SOS_CASES
from .write import complex_json, write_vector


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
