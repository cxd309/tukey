"""
Generate golden-vector JSON fixtures for tukey's Go tests.

Run with `uv run generate.py` from the testdata/ directory (or `uv run
testdata/generate.py` from the repo root). Requires the pinned scipy/numpy
in this directory's uv.lock — do not run against an unpinned environment.

Each JSON file records its own scipy/numpy version so a stale fixture is
easy to spot even if this script is re-run later against a newer pin.
"""

import json
from pathlib import Path
from typing import Final

import numpy as np
import scipy
from scipy import signal

OUTPUT_ROOT = Path(__file__).parent
VERSIONS = {"scipy_version": scipy.__version__, "numpy_version": np.__version__}

# -----------------------------------------------------------------------------
# Reference cases: everything tukey is checked against
# -----------------------------------------------------------------------------

SIGNAL_LENGTH: Final = 200

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

# name -> (b, a): filters shared by the filter and lfilter_zi vectors
FILTERS: Final = {
    "butter1_lp0.3": signal.butter(1, 0.3),  # shortest IIR
    "butter2_lp0.3": signal.butter(2, 0.3),  # typical call site
    "butter2_bp0.2-0.5": signal.butter(2, [0.2, 0.5], "bandpass"),
    "fir_avg4": ([0.25] * 4, [1.0]),  # a is length 1: no feedback
    "unnormalised": ([2.0, 1.0], [2.0, -0.5]),  # a[0] != 1
    "short_b": ([1.0], [1.0, -0.9]),  # len(b) < len(a)
}

# (name, padtype, padlen) variations of scipy.signal.filtfilt's edge handling,
# each run on butter2_lp0.3; padlen None means SciPy's 3 * max(len(a), len(b))
FILTFILT_PADDING_CASES = (
    ("even", "even", None),  # mirrored edges: slope flips
    ("constant", "constant", None),  # flat edges: slope drops to zero
    ("none", None, None),  # no extension: transients land on the data
    ("padlen0", "odd", 0),  # odd padtype but nothing added
    ("padlen_matlab", "odd", 6),  # MATLAB's 3 * (nfilt - 1)
    ("padlen50", "odd", 50),  # longer than the default 9
    ("padlen_max", "odd", SIGNAL_LENGTH - 1),  # longest SciPy allows
)


# -----------------------------------------------------------------------------
# Vector generators
# -----------------------------------------------------------------------------


def write_vector(
    category: str, name: str, description: str, params: dict, output: dict
) -> None:
    out_dir = OUTPUT_ROOT / category
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "description": description,
        **VERSIONS,
        "params": params,
        "output": output,
    }
    path = out_dir / f"{name}.json"
    path.write_text(json.dumps(payload, indent=2) + "\n")
    # print(f"wrote {path.relative_to(OUTPUT_ROOT)}")


def complex_json(values) -> dict:
    """Complex values as {"re":[...], "im":[...]}, since JSON has no complex type"""
    values = np.asarray(values, dtype=complex)
    return {"re": values.real.tolist(), "im": values.imag.tolist()}


def generate_butter() -> None:
    """
    b, a coefficients from scipy.signal.butter, analog=False.

    wn is normalised to Nyquist (1.0 == Nyquist), matching both scipy and
    MATLAB. Low/high-pass take a scalar wn; band-pass/stop take [low, high].
    """

    for order, wn, btype in BUTTER_LOWHIGH_CASES:
        b, a = signal.butter(order, wn, btype=btype, analog=False)
        name = f"{btype}_order{order}_wn{wn}"
        write_vector(
            "butter",
            name,
            f"scipy.signal.butter({order}, {wn}, btype='{btype}')",
            {"order": order, "wn": wn, "btype": btype},
            {
                "b": b.tolist(),
                "a": a.tolist(),
                "sos": signal.butter(order, wn, btype=btype, output="sos").tolist(),
            },
        )

    for order, wn, btype in BUTTER_BAND_CASES:
        b, a = signal.butter(order, wn, btype=btype, analog=False)
        name = f"{btype}_order{order}_wn{wn[0]}-{wn[1]}"
        write_vector(
            "butter",
            name,
            f"scipy.signal.butter({order}, {list(wn)}, btype='{btype}')",
            {"order": order, "wn": list(wn), "btype": btype},
            {
                "b": b.tolist(),
                "a": a.tolist(),
                "sos": signal.butter(order, wn, btype=btype, output="sos").tolist(),
            },
        )


def generate_zpk2sos() -> None:
    """sos from scipy.signal.zpk2sos: SciPy's own designs plus hand-built cases."""
    cases = dict(ZPK2SOS_CASES)
    for order, wn, btype in BUTTER_LOWHIGH_CASES + BUTTER_BAND_CASES:
        z, p, k = signal.butter(order, wn, btype=btype, output="zpk")
        wn_name = f"{wn[0]}-{wn[1]}" if isinstance(wn, tuple) else wn
        cases[f"butter_{btype}_order{order}_wn{wn_name}"] = (z, p, k)

    for name, (z, p, k) in cases.items():
        sos = signal.zpk2sos(z, p, k)
        write_vector(
            "zpk2sos",
            name,
            f"scipy.signal.zpk2sos(z, p, k) for {name}",
            {"z": complex_json(z), "p": complex_json(p), "k": float(k)},
            {"sos": sos.tolist()},
        )


def synthetic_signals(n: int) -> dict[str, np.ndarray]:
    t = np.arange(n) / n
    rng = np.random.default_rng(seed=0)  # fixed seed: reproducible noise
    impulse = np.zeros(n)
    impulse[0] = 1.0
    return {
        "impulse": impulse,  # output is the impulse response itself
        "step": np.r_[np.zeros(n // 2), np.ones(n - n // 2)],  # mid-signal edge
        "constant": np.full(n, 3.0),  # DC: filtfilt must preserve exactly
        "sines": np.sin(2 * np.pi * 5 * t) + 0.5 * np.sin(2 * np.pi * 40 * t),
        "chirp": signal.chirp(t, f0=1, t1=1, f1=80),  # sweeps through cutoff
        "noise": rng.standard_normal(n),  # broadband
    }


def generate_filter() -> None:
    """y from scipy.signal.lfilter(b, a, x) with the filter at rest."""

    for fname, (b, a) in FILTERS.items():
        for sname, x in synthetic_signals(SIGNAL_LENGTH).items():
            y = signal.lfilter(b, a, x)
            write_vector(
                "filter",
                f"{fname}_{sname}",
                f"scipy.signal.lfilter(b, a, x) for {fname}, {sname}",
                {"b": list(b), "a": list(a), "x": x.tolist()},
                {"y": y.tolist()},
            )


def generate_lfilter_zi() -> None:
    """zi from scipy.signal.lfilter_zi(b, a): the steady-state step response state."""
    for fname, (b, a) in FILTERS.items():
        zi = signal.lfilter_zi(b, a)
        write_vector(
            "lfilter_zi",
            fname,
            f"scipy.signal.lfilter_zi(b,a) for {fname}",
            {"b": list(b), "a": list(a)},
            {"zi": zi.tolist()},
        )


def generate_filtfilt() -> None:
    """y from scipy.signal.filtfilt(b, a, x) aross filters, signals and edge handling."""
    signals = synthetic_signals(SIGNAL_LENGTH)

    # default edge handling, for every filter and signal
    for fname, (b, a) in FILTERS.items():
        for sname, x in signals.items():
            write_filtfilt(f"{fname}_{sname}", b, a, x, "odd", None)
    # each edge-handling variation, on one typical filter
    b, a = FILTERS["butter2_lp0.3"]
    for case, padtype, padlen in FILTFILT_PADDING_CASES:
        for sname in ("sines", "step"):
            write_filtfilt(
                f"butter2_lp0.3_{sname}_{case}", b, a, signals[sname], padtype, padlen
            )
    # shortest signal the default padding allows: len(x) = padlen + 1 = 3*3 + 1
    write_filtfilt("butter2_1p0.3_shortest", b, a, signals["sines"][:10], "odd", None)


def write_filtfilt(
    name: str, b, a, x: np.ndarray, padtype: str | None, padlen: int | None
) -> None:
    y = signal.filtfilt(b, a, x, padtype=padtype, padlen=padlen)
    write_vector(
        "filtfilt",
        name,
        f"scipy.signal.filtfilt(b, a, x, padtype={padtype!r}, padlen={padlen!r})",
        {
            "b": list(b),
            "a": list(a),
            "x": x.tolist(),
            "padtype": padtype,
            "padlen": padlen,
        },
        {"y": y.tolist()},
    )


# -----------------------------------------------------------------------------
# Main function: run all the tests
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    generate_butter()
    print("wrote butter/")
    generate_zpk2sos()
    print("wrote zpk2sos/")
    generate_filter()
    print("wrote filter/")
    generate_lfilter_zi()
    print("wrote lfilter_zi/")
    generate_filtfilt()
    print("wrote filtfilt/")
