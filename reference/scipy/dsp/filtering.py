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


# MATLAB's filtfilt, from its source (filtfilt.m, R2025a): every section is
# filtered forward and backward in turn, each odd-extended by the same edge
# length, max(1, 3·order), where the order is summed over the sections, each
# counting its longer polynomial's last nonzero coefficient. That's SciPy's
# filtfilt per section with padlen set to it, so the reference is built from
# SciPy's, and MATLAB's own vectors check the definition.


def effective_order(b, a) -> int:
    """a filter's order as MATLAB counts it: the longer polynomial's last
    nonzero coefficient's index, trailing zeros dropped"""

    def last_nonzero(c) -> int:
        nonzero = np.flatnonzero(np.asarray(c, dtype=float))
        return int(nonzero[-1]) if nonzero.size else 0

    return max(last_nonzero(b), last_nonzero(a))


def matlab_padlen(sections) -> int:
    """MATLAB filtfilt's edge length for a filter of (b, a) sections"""
    return max(1, 3 * sum(effective_order(b, a) for b, a in sections))


def scipy_ntaps(sos: np.ndarray) -> int:
    """scipy.signal.sosfiltfilt's tap count, for its default padding of 3·ntaps"""
    first_order = min(int(np.sum(sos[:, 2] == 0)), int(np.sum(sos[:, 5] == 0)))
    return 2 * len(sos) + 1 - first_order


def filtfilt_by_section(
    sos: np.ndarray, x: np.ndarray, padtype: str | None, padlen: int
) -> np.ndarray:
    """each section forward and backward in turn, each padded by padlen"""
    y = x
    for section in sos:
        y = sc.signal.filtfilt(
            section[:3], section[3:], y, padtype=padtype, padlen=padlen
        )
    return y


# sections whose order MATLAB counts differently from SciPy's tap count,
# name -> (sos, SciPy's ntaps - 1, MATLAB's order), each asserted
ORDER_COUNT_FILTERS = {
    # b2 = 0 in one section and a2 = 0 in the other: SciPy counts one
    # first-order section, MATLAB none: order 3 against 4
    "mixed_zeros": (
        np.array([[0.2, 0.3, 0.0, 1.0, -0.5, 0.2], [0.4, 0.1, 0.05, 1.0, -0.3, 0.0]]),
        3,
        4,
    ),
    # a gain-only section is order 0 to MATLAB, 1 to SciPy: 2 against 3
    "gain_section": (
        np.array(
            [[2.0, 0.0, 0.0, 1.0, 0.0, 0.0], list(SOS_FILTERS["butter3_lp0.3"][1])]
        ),
        3,
        2,
    ),
    # a pure gain: order 0, so MATLAB's edge is its minimum, 1 sample
    "gain_only": (np.array([[2.0, 0.0, 0.0, 1.0, 0.0, 0.0]]), 1, 0),
}


def generate_sosfiltfilt_by_section() -> None:
    """SOSFiltFilt's BySection: SciPy's filtfilt on each section in turn"""
    signals = synthetic_signals(SIGNAL_LENGTH)

    def write(name, sos, x, padtype, padlen, matlab):
        if matlab:
            n = matlab_padlen([(s[:3], s[3:]) for s in sos])
        elif padlen is None:
            n = 3 * scipy_ntaps(sos)
        else:
            n = padlen
        write_vector(
            "sosfiltfilt_by_section",
            name,
            f"scipy.signal.filtfilt on each section in turn, padtype={padtype!r}, padlen={n}",
            {
                "sos": sos.tolist(),
                "x": x.tolist(),
                "padtype": padtype,
                "padlen": padlen,
                "matlab_padlen": matlab,
            },
            {"y": filtfilt_by_section(sos, x, padtype, n).tolist()},
        )

    # MATLAB's filtfilt(sos, 1, x), for every filter and signal
    for fname, sos in SOS_FILTERS.items():
        for sname, x in signals.items():
            write(f"{fname}_{sname}_matlab", sos, x, "odd", None, True)
        # SciPy's default padding, 3·ntaps, by section
        write(f"{fname}_sines", sos, signals["sines"], "odd", None, False)
    # each edge-handling variation, on the filter with a first-order section
    sos = SOS_FILTERS["butter3_lp0.3"]
    for case, padtype, padlen in SOSFILTFILT_PADDING_CASES:
        for sname in ("sines", "step"):
            write(
                f"butter3_lp0.3_{sname}_{case}",
                sos,
                signals[sname],
                padtype,
                padlen,
                False,
            )
    # shortest signal MATLAB allows. Its first section has the real pole and
    # two zeros, order 2 to MATLAB, where SciPy counts it first-order: order 4,
    # an edge of 12, so 13 samples
    n = matlab_padlen([(s[:3], s[3:]) for s in sos])
    assert n == 12, n
    write(
        "butter3_lp0.3_shortest_matlab",
        sos,
        signals["sines"][: n + 1],
        "odd",
        None,
        True,
    )
    # sections whose order MATLAB counts its own way
    for fname, (sos, scipy_order, order) in ORDER_COUNT_FILTERS.items():
        assert scipy_ntaps(sos) - 1 == scipy_order, fname
        assert sum(effective_order(s[:3], s[3:]) for s in sos) == order, fname
        for sname in ("sines", "step"):
            write(f"{fname}_{sname}_matlab", sos, signals[sname], "odd", None, True)


# (b, a) filters whose order MATLAB counts by dropping trailing zeros
TRAILING_ZERO_FILTERS = {
    "trailing_b": ([0.3, 0.2, 0.0], [1.0, -0.5, 0.0]),  # order 1, not 2
    "trailing_a": ([0.25, 0.5, 0.25], [1.0, -0.4, 0.0]),  # order 2: b is longer
    "fir_trailing": ([0.25, 0.5, 0.25, 0.0], [1.0]),  # order 2, not 3
    "gain": ([2.0], [1.0]),  # order 0: an edge of 1 sample
    "butter2_lp0.3": FILTERS["butter2_lp0.3"],  # the usual case, order 2
}


def generate_filtfilt_matlab() -> None:
    """FiltFilt with MATLABPadLen: SciPy's filtfilt with MATLAB's edge length"""
    signals = synthetic_signals(SIGNAL_LENGTH)
    for fname, (b, a) in TRAILING_ZERO_FILTERS.items():
        n = matlab_padlen([(b, a)])
        for sname in ("sines", "step"):
            x = signals[sname]
            write_vector(
                "filtfilt_matlab",
                f"{fname}_{sname}",
                f"scipy.signal.filtfilt(b, a, x, padlen={n}), MATLAB's edge length",
                {"b": list(b), "a": list(a), "x": x.tolist()},
                {"y": sc.signal.filtfilt(b, a, x, padlen=n).tolist()},
            )
