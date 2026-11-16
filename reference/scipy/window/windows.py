"""scipy.signal.windows, symetric and period, at lengths covering the edge cases"""

from typing import Final

import scipy as sc

from write import write_vector

# empty and single-sample (SciPy returns [] and [1]),
# the smallest lengths where the shape shows
# then even/odd pairs where symmetric and periodic differ
LENGTHS: Final = (0, 1, 2, 3, 8, 9, 64, 255)

# HFT90D flat-top coefficients, from scipy.signal.windows.general_cosine's docs
HFT90D: Final = [1, 1.942604, 1.340318, 0.440811, 0.043097]

# window -> {case label: keyword arguments}, "" for windows with no parameters
WINDOWS: Final = {
    "boxcar": {"": {}},
    "triang": {"": {}},
    "bartlett": {"": {}},
    "hann": {"": {}},
    "hamming": {"": {}},
    "blackman": {"": {}},
    "blackmanharris": {"": {}},
    "nuttall": {"": {}},
    "flattop": {"": {}},
    "cosine": {"": {}},
    "barthann": {"": {}},
    "bohman": {"": {}},
    "parzen": {"": {}},
    "lanczos": {"": {}},
    # alpha 0 is boxcar, 1 is hann
    "tukey": {f"alpha{a}": {"alpha": a} for a in (0.0, 0.25, 0.5, 1.0)},
    # beta 0 is boxcar; 8.6 approximates blackman
    "kaiser": {f"beta{b}": {"beta": b} for b in (0.0, 5.0, 8.6, 14.0)},
    # even length and symmetric only: see unsupported
    "kaiser_bessel_derived": {f"beta{b}": {"beta": b} for b in (4.0, 12.0)},
    "gaussian": {f"std{s}": {"std": s} for s in (1.0, 2.5, 7.0)},
    # p 1 is gaussian, p 1.5, sig 7 is SciPy's docs example, p 4 is nearly flat-topped
    "general_gaussian": {
        f"p{p}_sig{sig}": {"p": p, "sig": sig}
        for p, sig in ((1.0, 2.5), (1.5, 7.0), (4.0, 20.0))
    },
    # [0.5, 0.5] is hann
    "general_cosine": {"hann": {"a": [0.5, 0.5]}, "hft90d": {"a": HFT90D}},
    # 0.5 is hann, 0.54 is hamming
    "general_hamming": {f"alpha{a}": {"alpha": a} for a in (0.5, 0.54, 0.75)},
    # center defaults to the middle, moving it is periodic only, see unsupported
    "exponential": {
        **{f"tau{t}": {"tau": t} for t in (1.0, 3.0, 10.0)},
        "center0.0_tau3.0": {"center": 0.0, "tau": 3.0},
    },
    # MATLAB's sll is negative dB, and it has no norm options:
    # both are generated so MATLAB's vectors show which it matches
    "taylor": {
        f"nbar{nbar}_sll{sll}_{'norm' if norm else 'raw'}": {
            "nbar": nbar,
            "sll": sll,
            "norm": norm,
        }
        for nbar, sll in ((4, 30), (9, 50))
        for norm in (True, False)
    },
    # SciPy warns below 45dB
    "chebwin": {f"at{at}": {"at": at} for at in (50, 100)},
}


def generate_windows() -> None:
    """Every window in WINDOWS, at every length, symmetric and periodic, where scipy allows"""
    for name, cases in WINDOWS.items():
        window = getattr(sc.signal.windows, name)
        for label, kwargs in cases.items():
            for n in LENGTHS:
                for sym in (True, False):
                    if unsupported(name, kwargs, n, sym):
                        continue
                    shape = "sym" if sym else "periodic"
                    args = "".join(f"{k}={v!r}, " for k, v in kwargs.items())
                    write_vector(
                        f"window/{name}",
                        "_".join(filter(None, (label, f"n{n}", shape))),
                        f"scipy.signal.windows.{name}({n}, {args}sym={sym})",
                        {"n": n, "sym": sym, **kwargs},
                        {"w": window(n, **kwargs, sym=sym).tolist()},
                    )


def unsupported(name: str, kwargs: dict, n: int, sym: bool) -> bool:
    """Whether SciPy rejects this case with a ValueError, so there's nothing to check"""
    if name == "kaiser_bessel_derived":
        return not sym or n % 2 == 1
    if name == "exponential":
        return sym and "center" in kwargs
    return False
