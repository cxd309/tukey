"""Frequency response and group delay of (b, a) filters and second-order sections"""

from typing import Final

import numpy as np
import scipy as sc

from dsp.cases import FILTERS, SOS_FILTERS
from dsp.delay import exact_group_delay
from write import complex_json, write_vector


def symmetric_firwin() -> list:
    """firwin(31, 0.3)'s taps made exactly symmetric, so its delay is exactly 15"""
    h = sc.signal.firwin(31, 0.3)
    return np.concatenate([h[:16], h[:15][::-1]]).tolist()


# FILTERS, then filters whose numerator has zeros on or near the unit circle,
# where the group delay's two small polynomials lose digits: scipy.signal's
# group_delay is off by 3% on butter8_lp0.01 near Nyquist
GROUP_DELAY_FILTERS: Final = {
    **FILTERS,
    "butter8_lp0.01": sc.signal.butter(8, 0.01),  # eightfold zero at z = -1
    "ellip4_lp0.3": sc.signal.ellip(4, 1, 40, 0.3),  # stopband zeros on the circle
    "notch0.3_q30": sc.signal.iirnotch(0.3, 30),  # zeros on the circle at 0.3
    "fir_hamming31_lp0.3": (symmetric_firwin(), [1.0]),  # linear phase: 15 everywhere
    "difference": ([1.0, -1.0], [1.0]),  # a zero exactly at DC, a grid point
}


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


def generate_group_delay() -> None:
    """
    The exact group delay (see dsp.delay) of every filter in
    GROUP_DELAY_FILTERS at 512 points, and of every filter in SOS_FILTERS, as
    the sum of its sections' delays. Each vector lists the points where a zero
    or pole is on the unit circle, where the delay is its limit from both
    sides. SciPy isn't asserted against: off by up to 3% away from those
    points, it's the formula these vectors check against, not a reference.
    """
    for fname, (b, a) in GROUP_DELAY_FILTERS.items():
        write_group_delay(f"{fname}_n512", list(b), list(a), 512)
    b, a = FILTERS["butter2_lp0.3"]
    for n in (0, 1, 7):  # empty, DC only, odd
        write_group_delay(f"butter2_lp0.3_n{n}", list(b), list(a), n)

    for fname, sos in SOS_FILTERS.items():
        sections = [(row[:3].tolist(), row[3:].tolist()) for row in sos]
        w, delay, singular = exact_group_delay(sections, 512)
        write_vector(
            "group_delay_sos",
            f"{fname}_n512",
            f"the exact group delay of {fname}'s sections at 512 points, as "
            "scipy.signal.group_delay of their product",
            {"sos": sos.tolist(), "n": 512},
            {"w": w, "delay": delay, "singular": singular},
        )


def write_group_delay(name: str, b: list, a: list, n: int) -> None:
    w, delay, singular = exact_group_delay([(b, a)], n)
    write_vector(
        "group_delay",
        name,
        f"the exact group delay, as scipy.signal.group_delay((b, a), w={n}) for {name}",
        {"b": b, "a": a, "n": n},
        {"w": w, "delay": delay, "singular": singular},
    )
