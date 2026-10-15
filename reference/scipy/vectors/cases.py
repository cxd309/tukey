"""Every case tukey is checked against: the tables the generators iterate over"""

from typing import Final

import scipy as sc

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
    "butter1_lp0.3": sc.signal.butter(1, 0.3),  # shortest IIR
    "butter2_lp0.3": sc.signal.butter(2, 0.3),  # typical call site
    "butter2_bp0.2-0.5": sc.signal.butter(2, [0.2, 0.5], "bandpass"),
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

# filtfilt's padding variations, minus MATLAB's padlen, which was computed for 3 taps
SOSFILTFILT_PADDING_CASES = tuple(
    c for c in FILTFILT_PADDING_CASES if c[0] != "padlen_matlab"
)

# lengths the FFT must handle, chosen to exercise each algorithm path
FFT_LENGTHS = (
    1,
    2,
    4,
    8,
    16,
    128,
    1024,
    4096,  # powers of two: radix-2
    3,
    5,
    9,
    15,
    30,
    60,
    100,
    1000,  # smooth: mixed radix over 2, 3, 5
    6,
    12,
    210,  # more mixed factors; 210 = 2*3*5*7 brings in a 7
    7,
    17,
    97,
    1009,  # primes: the generic butterfly or Bluestein
)
