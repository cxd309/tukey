"""Cases shared across the dsp generators: signal length and filters under test"""

from typing import Final

import scipy as sc

SIGNAL_LENGTH: Final = 200

# name -> sos: filters as second-order sections, for sosfilt, sosfilt_zi, sosfiltfilt
SOS_FILTERS = {
    "butter3_lp0.3": sc.signal.butter(
        3, 0.3, output="sos"
    ),  # odd order: a first-order section
    "butter4_bs0.2-0.5": sc.signal.butter(
        4, [0.2, 0.5], "bandstop", output="sos"
    ),  # 4 sections
    "butter8_lp0.01": sc.signal.butter(
        8, 0.01, output="sos"
    ),  # overall gain ~1e-14 in section 0
    "butter10_bp0.3-0.31": sc.signal.butter(
        10, [0.3, 0.31], "bandpass", output="sos"
    ),  # order 20: (b, a) fails here
}

# name -> (b, a): filters shared by the filter and lfilter_zi vectors
FILTERS: Final = {
    "butter1_lp0.3": sc.signal.butter(1, 0.3),  # shortest IIR
    "butter2_lp0.3": sc.signal.butter(2, 0.3),  # typical call site
    "butter2_bp0.2-0.5": sc.signal.butter(2, [0.2, 0.5], "bandpass"),
    "fir_avg4": ([0.25] * 4, [1.0]),  # a is length 1: no feedback
    "unnormalised": ([2.0, 1.0], [2.0, -0.5]),  # a[0] != 1
    "short_b": ([1.0], [1.0, -0.9]),  # len(b) < len(a)
}
