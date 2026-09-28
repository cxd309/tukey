"""
Generate golden-vector JSON fixtures for godsp's Go tests.

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
# Reference cases: everything godsp is checked against
# -----------------------------------------------------------------------------

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

# name -> (b, a): filters shared by the filter and lfilter_zi vectors
FILTERS: Final = {
    "butter1_lp0.3": signal.butter(1, 0.3),  # shortest IIR
    "butter2_lp0.3": signal.butter(2, 0.3),  # typical call site
    "butter2_bp0.2-0.5": signal.butter(2, [0.2, 0.5], "bandpass"),
    "fir_avg4": ([0.25] * 4, [1.0]),  # a is length 1: no feedback
    "unnormalised": ([2.0, 1.0], [2.0, -0.5]),  # a[0] != 1
    "short_b": ([1.0], [1.0, -0.9]),  # len(b) < len(a)
}

SIGNAL_LENGTH: Final = 200

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
            {"b": b.tolist(), "a": a.tolist()},
        )

    for order, wn, btype in BUTTER_BAND_CASES:
        b, a = signal.butter(order, wn, btype=btype, analog=False)
        name = f"{btype}_order{order}_wn{wn[0]}-{wn[1]}"
        write_vector(
            "butter",
            name,
            f"scipy.signal.butter({order}, {list(wn)}, btype='{btype}')",
            {"order": order, "wn": list(wn), "btype": btype},
            {"b": b.tolist(), "a": a.tolist()},
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


# -----------------------------------------------------------------------------
# Main function: run all the tests
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    generate_butter()
    print("wrote butter/")
    generate_filter()
    print("wrote filter/")
    generate_lfilter_zi()
    print("wrote lfilter_zi/")
