"""Filtering with (b, a) and second-order sections: at rest, from a state, and zero-phase"""

import numpy as np
import scipy as sc

from dsp.cases import FILTERS, SIGNAL_LENGTH, SOS_FILTERS
from dsp.signals import synthetic_signals
from write import write_vector

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


# filtfilt's padding variations, minus MATLAB's padlen, which was computed for 3 taps
SOSFILTFILT_PADDING_CASES = tuple(
    c for c in FILTFILT_PADDING_CASES if c[0] != "padlen_matlab"
)


def generate_filter() -> None:
    """y from scipy.signal.lfilter(b, a, x) with the filter at rest."""

    for fname, (b, a) in FILTERS.items():
        for sname, x in synthetic_signals(SIGNAL_LENGTH).items():
            y = sc.signal.lfilter(b, a, x)
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
        zi = sc.signal.lfilter_zi(b, a)
        write_vector(
            "lfilter_zi",
            fname,
            f"scipy.signal.lfilter_zi(b,a) for {fname}",
            {"b": list(b), "a": list(a)},
            {"zi": zi.tolist()},
        )


def generate_lfilter_state() -> None:
    """lfilter from an arbitrary starting state, to check zi in and zf out"""
    x = synthetic_signals(SIGNAL_LENGTH)["sines"]
    for fname, (b, a) in FILTERS.items():
        zi = 0.5 * sc.signal.lfilter_zi(b, a)
        y, zf = sc.signal.lfilter(b, a, x, zi=zi)
        write_vector(
            "lfilter_state",
            fname,
            f"scipy.signal.lfilter(b, a, x, zi=0.5*lfilter_zi(b, a)) for {fname}, sines",
            {"b": list(b), "a": list(a), "x": x.tolist(), "zi": zi.tolist()},
            {"y": y.tolist(), "zf": zf.tolist()},
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
    write_filtfilt("butter2_lp0.3_shortest", b, a, signals["sines"][:10], "odd", None)


def write_filtfilt(
    name: str, b, a, x: np.ndarray, padtype: str | None, padlen: int | None
) -> None:
    y = sc.signal.filtfilt(b, a, x, padtype=padtype, padlen=padlen)
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


def generate_sos_filtering() -> None:
    """sosfilt, sosfilt_zi and sosfiltfilt on filters given as second-order sections."""
    signals = synthetic_signals(SIGNAL_LENGTH)

    for fname, sos in SOS_FILTERS.items():
        zi = sc.signal.sosfilt_zi(sos)
        write_vector(
            "sosfilt_zi",
            fname,
            f"scipy.signal.sosfilt_zi(sos) for {fname}",
            {"sos": sos.tolist()},
            {"zi": zi.tolist()},
        )

        for sname, x in signals.items():
            write_vector(
                "sosfilt",
                f"{fname}_{sname}",
                f"scipy.signal.sosfilt(sos, x) for {fname}, {sname}",
                {"sos": sos.tolist(), "x": x.tolist(), "zi": None},
                {"y": sc.signal.sosfilt(sos, x).tolist(), "zf": None},
            )
            write_sosfiltfilt(f"{fname}_{sname}", sos, x, "odd", None)

        # an arbitrary starting state, to exercise zi in and zf out
        zi0 = 0.5 * zi
        y, zf = sc.signal.sosfilt(sos, signals["sines"], zi=zi0)
        write_vector(
            "sosfilt",
            f"{fname}_sines_with_state",
            f"scipy.signal.sosfilt(sos, x, zi=0.5*sosfilt_zi(sos)) for {fname}, sines",
            {"sos": sos.tolist(), "x": signals["sines"].tolist(), "zi": zi0.tolist()},
            {"y": y.tolist(), "zf": zf.tolist()},
        )

    # each edge-handling variation, on the filter with a first-order section
    sos = SOS_FILTERS["butter3_lp0.3"]
    for case, padtype, padlen in SOSFILTFILT_PADDING_CASES:
        for sname in ("sines", "step"):
            write_sosfiltfilt(
                f"butter3_lp0.3_{sname}_{case}", sos, signals[sname], padtype, padlen
            )

    # shortest signal the default padding allows: 2 sections, one first-order,
    # so ntaps = 2*2 + 1 - 1 = 4 and padlen = 12
    write_sosfiltfilt("butter3_lp0.3_shortest", sos, signals["sines"][:13], "odd", None)


def write_sosfiltfilt(
    name: str, sos: np.ndarray, x: np.ndarray, padtype: str | None, padlen: int | None
) -> None:
    y = sc.signal.sosfiltfilt(sos, x, padtype=padtype, padlen=padlen)
    write_vector(
        "sosfiltfilt",
        name,
        f"scipy.signal.sosfiltfilt(sos, x, padtype={padtype!r}, padlen={padlen!r})",
        {"sos": sos.tolist(), "x": x.tolist(), "padtype": padtype, "padlen": padlen},
        {"y": y.tolist()},
    )
