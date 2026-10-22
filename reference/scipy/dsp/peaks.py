"""Peak finding: find_peaks with each option, peak_prominences and peak_widths"""

from typing import Final

import numpy as np
import scipy as sc

from dsp.cases import SIGNAL_LENGTH
from dsp.signals import synthetic_signals
from write import write_vector

# find_peaks option sets, each applied to the PEAK_SIGNALS
# "" is no options: every local maximum
FIND_PEAKS_OPTIONS: Final = {
    "": {},
    "distance5": {"distance": 5},
    "distance20": {
        "distance": 20
    },  # the target app uses distances of 80 and 400 on long signals
    "height": {"height": 0.5},
    "threshold": {"threshold": 0.05},
    "prominence": {"prominence": 0.5},
    "width": {"width": 3},
    "distance_prominence": {"distance": 10, "prominence": 0.2},
}

# small hand-built signals for the edge cases, checked with no options
# and with a distance, which is where MATLAB and SciPy are most likely to differ
EDGE_SIGNALS: Final = {
    "plateaus": [
        0.0,
        1,
        1,
        0,
        2,
        2,
        2,
        0,
        1,
        0,
        3,
        3,
        0,
    ],  # flat tops of 2 and 3 samples
    "ties": [0.0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0],  # equal peaks closer than the distance
    "ends_higher": [3.0, 1, 2, 1, 3],  # the ends are never peaks
    "constant": [1.0] * 6,  # no peaks
    "empty": [],
    "one": [1.0],
    "two": [1.0, 2],
    "three": [1.0, 2, 1],
    # a tall peak of low prominence (4.9 at 5) next to a shorter one of high
    # prominence (4.5 at 7): with distance 3 and prominence 0.3, applying the
    # distance rule first loses both, applying prominence first keeps 7
    "distance_vs_prominence": [0.0, 5, 4.8, 4.8, 4.8, 4.9, 4.0, 4.5, 0],
    # a tall narrow spike (5 at 2, width 1) next to a shorter wide hump (3 at
    # 6, width 3): with distance 5 and width 2, applying the distance rule
    # first loses both, applying width first keeps 6
    "distance_vs_width": [0.0, 0, 5, 0, 1, 2, 3, 2, 1, 0, 0],
}


def peak_signals() -> dict[str, np.ndarray]:
    """the signals every find_peaks option is checked on"""
    signals = synthetic_signals(SIGNAL_LENGTH)
    rng = np.random.default_rng(seed=3)
    return {
        "sines": signals["sines"],
        "chirp": signals["chirp"],
        "noise": signals["noise"],
        "abs_noise": np.abs(rng.standard_normal(500)),  # as the app: peaks of |x|
    }


def generate_find_peaks() -> None:
    """scipy.signal.find_peaks on each signal with each option set, and the edge cases"""
    cases = [
        (f"{sname}_{label}" if label else sname, x, kwargs)
        for sname, x in peak_signals().items()
        for label, kwargs in FIND_PEAKS_OPTIONS.items()
    ]
    for name, values in EDGE_SIGNALS.items():
        x = np.array(values)
        cases.append((name, x, {}))
        cases.append((f"{name}_distance3", x, {"distance": 3}))
    # each pair of options together, which decides the order they're applied in
    cases.append(
        (
            "distance_vs_prominence_both",
            np.array(EDGE_SIGNALS["distance_vs_prominence"]),
            {"distance": 3, "prominence": 0.3},
        )
    )
    cases.append(
        (
            "distance_vs_width_both",
            np.array(EDGE_SIGNALS["distance_vs_width"]),
            {"distance": 5, "width": 2},
        )
    )

    for name, x, kwargs in cases:
        peaks, _ = sc.signal.find_peaks(x, **kwargs)
        args = "".join(f", {k}={v!r}" for k, v in kwargs.items())
        write_vector(
            "find_peaks",
            name,
            f"scipy.signal.find_peaks(x{args})",
            {"x": x.tolist(), **kwargs},
            {"peaks": peaks.tolist(), "heights": x[peaks].tolist()},
        )


def generate_peak_prominences() -> None:
    """scipy.signal.peak_prominences for every peak find_peaks finds"""
    for name, x in peak_signals().items():
        peaks, _ = sc.signal.find_peaks(x)
        prominences, left, right = sc.signal.peak_prominences(x, peaks)
        write_vector(
            "peak_prominences",
            name,
            "scipy.signal.peak_prominences(x, find_peaks(x)[0])",
            {"x": x.tolist(), "peaks": peaks.tolist()},
            {
                "prominences": prominences.tolist(),
                "left_bases": left.tolist(),
                "right_bases": right.tolist(),
            },
        )


def generate_peak_widths() -> None:
    """scipy.signal.peak_widths for every peak find_peaks finds, at half and full prominence"""
    for name, x in peak_signals().items():
        peaks, _ = sc.signal.find_peaks(x)
        for rel_height in (0.5, 1.0):
            widths, heights, left, right = sc.signal.peak_widths(
                x, peaks, rel_height=rel_height
            )
            write_vector(
                "peak_widths",
                f"{name}_rel{rel_height}",
                f"scipy.signal.peak_widths(x, find_peaks(x)[0], rel_height={rel_height})",
                {"x": x.tolist(), "peaks": peaks.tolist(), "rel_height": rel_height},
                {
                    "widths": widths.tolist(),
                    "width_heights": heights.tolist(),
                    "left_ips": left.tolist(),
                    "right_ips": right.tolist(),
                },
            )
