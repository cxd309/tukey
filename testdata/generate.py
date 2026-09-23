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

import numpy as np
import scipy
from scipy import signal

OUTPUT_ROOT = Path(__file__).parent
VERSIONS = {"scipy_version": scipy.__version__, "numpy_version": np.__version__}


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
    print(f"wrote {path.relative_to(OUTPUT_ROOT)}")


def generate_butter() -> None:
    """
    b, a coefficients from scipy.signal.butter, analog=False.

    wn is normalised to Nyquist (1.0 == Nyquist), matching both scipy and
    MATLAB. Low/high-pass take a scalar wn; band-pass/stop take [low, high].
    """
    lowhigh_cases = [
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
    for order, wn, btype in lowhigh_cases:
        b, a = signal.butter(order, wn, btype=btype, analog=False)
        name = f"{btype}_order{order}_wn{wn}"
        write_vector(
            "butter",
            name,
            f"scipy.signal.butter({order}, {wn}, btype='{btype}')",
            {"order": order, "wn": wn, "btype": btype},
            {"b": b.tolist(), "a": a.tolist()},
        )

    band_cases = [
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
    for order, wn, btype in band_cases:
        b, a = signal.butter(order, wn, btype=btype, analog=False)
        name = f"{btype}_order{order}_wn{wn[0]}-{wn[1]}"
        write_vector(
            "butter",
            name,
            f"scipy.signal.butter({order}, {list(wn)}, btype='{btype}')",
            {"order": order, "wn": list(wn), "btype": btype},
            {"b": b.tolist(), "a": a.tolist()},
        )


if __name__ == "__main__":
    generate_butter()
