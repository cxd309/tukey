"""Frequency response of (b, a) filters and second-order sections"""

import scipy as sc

from dsp.cases import FILTERS, SOS_FILTERS
from write import complex_json, write_vector


def generate_freqz() -> None:
    """freqz on every filter in FILTERS, and freqz_sos on every filter in SOS_FILTERS."""
    for fname, (b, a) in FILTERS.items():
        write_freqz("freqz", f"{fname}_n512", b, a, 512)
    # point counts at the edges, on one typical filter
    b, a = FILTERS["butter2_lp0.3"]
    for n in (0, 1, 7):  # empty, DC only, odd
        write_freqz("freqz", f"butter2_lp0.3_n{n}", b, a, n)

    for fname, sos in SOS_FILTERS.items():
        w, h = sc.signal.freqz_sos(sos, worN=512)
        write_vector(
            "freqz_sos",
            f"{fname}_n512",
            f"scipy.signal.freqz_sos(sos, worN=512) for {fname}",
            {"sos": sos.tolist(), "n": 512},
            {"w": w.tolist(), "h": complex_json(h)},
        )


def write_freqz(category: str, name: str, b, a, n: int) -> None:
    w, h = sc.signal.freqz(b, a, worN=n)
    write_vector(
        category,
        name,
        f"scipy.signal.freqz(b, a, worN={n}) for {name}",
        {"b": list(b), "a": list(a), "n": n},
        {"w": w.tolist(), "h": complex_json(h)},
    )
