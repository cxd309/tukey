"""
The exact window-method FIR design, at 50 significant digits with mpmath

scipy.signal.firwin's definition: for numtaps taps, centred on
m = k - (numtaps - 1)/2, the ideal response's impulse response is

    Σ over passbands (right·sinc(right·m) - left·sinc(left·m))

with sinc(x) = sin(pi·x)/(pi·x), multiplied by the window. Scaled, it's
divided by its amplitude Σ h·cos(pi·m·f) at f, the first passband's centre,
or 0 or 1 if it starts at DC or ends at Nyquist, so the response there is
exactly 1.

The window's samples are the input, as given, and only the design is exact:
SciPy's symmetric windows are themselves symmetric only to a few ULP, so the
taps are too. A scaled design is checked to be exactly 1 at its scale
frequency.
"""

from typing import Final

import mpmath as mp

from dsp.elliptic import DPS

SPEC_TOLERANCE: Final = mp.mpf("1e-40")  # how far a design may miss its spec


def passbands(cutoff: list, btype: str) -> list:
    """(left, right) of each passband, normalised to Nyquist"""
    edges = [mp.mpf(c) for c in cutoff]
    if btype == "lowpass":
        return [(mp.mpf(0), edges[0])]
    if btype == "highpass":
        return [(edges[0], mp.mpf(1))]
    if btype == "bandpass":
        return [(edges[0], edges[1])]
    return [(mp.mpf(0), edges[0]), (edges[1], mp.mpf(1))]


def scale_frequency(bands: list) -> mp.mpf:
    """where a scaled design's response is 1: the first passband's centre,
    or DC or Nyquist if it reaches either"""
    left, right = bands[0]
    if left == 0:
        return mp.mpf(0)
    if right == 1:
        return mp.mpf(1)
    return (left + right) / 2


def amplitude(h: list, m: list, f: mp.mpf) -> mp.mpf:
    """the zero-phase amplitude at f, as firwin scales by"""
    return mp.fsum(hk * mp.cospi(mk * f) for hk, mk in zip(h, m, strict=True))


def exact_firwin(window: list, cutoff: list, btype: str, scale: bool) -> list:
    """the taps, exact rounded once, after checking the design"""
    with mp.workdps(DPS):
        n = len(window)
        m = [k - mp.mpf(n - 1) / 2 for k in range(n)]
        bands = passbands(cutoff, btype)
        h = [
            mp.mpf(w)
            * mp.fsum(
                r * mp.sinc(mp.pi * r * mk) - left * mp.sinc(mp.pi * left * mk)
                for left, r in bands
            )
            for w, mk in zip(window, m, strict=True)
        ]
        if scale:
            f = scale_frequency(bands)
            s = amplitude(h, m, f)
            h = [hk / s for hk in h]
            assert abs(amplitude(h, m, f) - 1) < SPEC_TOLERANCE, (cutoff, btype)
        return [float(hk) for hk in h]
