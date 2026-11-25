"""
The exact Gaussian pulse, at 50 significant digits with mpmath

scipy.signal.gausspulse is exp(-a·t²)·cos(2·pi·fc·t), with its quadrature part
exp(-a·t²)·sin(2·pi·fc·t) and its envelope exp(-a·t²).
a sets the bandwidth: the spectrum, a Gaussian about fc, falls to bwr dB at
fc·(1 ± bw/2), which gives a = 5·(pi·fc·bw)²/(-bwr·ln(10)).
The cutoff is the time where the envelope falls to tpr dB,
tc = ln(10)·sqrt(tpr·bwr)/(10·pi·fc·bw).

Each is checked against its specification, independently of the formulas:
the spectrum at fc·bw/2 from fc is bwr dB below its peak, and the envelope at
tc is tpr dB below its peak.
"""

from typing import Final

import mpmath as mp
import numpy as np
import scipy as sc

from write import write_vector

DPS: Final = 50  # working precision, in significant digits
SPEC_TOLERANCE: Final = mp.mpf("1e-40")  # how far a pulse may miss its spec, relative

# how far scipy may be from exact: the pulse absolutely, at most 1.4e-15,
# and the cutoff relatively, at most 2.2e-16
SCIPY_ACCURACY: Final = 1e-14
SCIPY_CUTOFF_ACCURACY: Final = 1e-15

# the sample times, and the pulse over them, as (t, fc, bw, bwr)
# unit: scipy's own example, 5 Hz over -1 to 1 seconds
# defaults: scipy's default 1 kHz over 2 ms either side
# ultrasound: MATLAB's example, 50 kHz at 60% bandwidth, 100 ns apart
# narrow: a 10% bandwidth at -3 dB, many cycles under the envelope
# wide: a 200% bandwidth at -20 dB, about one cycle
# tails: 1 kHz out to 40 ms, where the envelope underflows to subnormals and then 0
# empty: no samples
GRIDS: Final = {
    "unit": (np.linspace(-1, 1, 2001), 5.0, 0.5, -6.0),
    "defaults": (np.arange(-1000, 1001) * 2e-6, 1000.0, 0.5, -6.0),
    "ultrasound": (np.arange(-300, 301) * 1e-7, 50e3, 0.6, -6.0),
    "narrow": (np.linspace(-1, 1, 2001), 5.0, 0.1, -3.0),
    "wide": (np.linspace(-1, 1, 2001), 5.0, 2.0, -20.0),
    "tails": (np.linspace(-0.04, 0.04, 801), 1000.0, 0.5, -6.0),
    "empty": (np.array([]), 1000.0, 0.5, -6.0),
}

# the cutoffs, as (fc, bw, bwr, tpr)
CUTOFFS: Final = [
    (1000.0, 0.5, -6.0, -60.0),  # scipy's defaults
    (50e3, 0.6, -6.0, -40.0),  # MATLAB's example
    (5.0, 0.1, -3.0, -100.0),
    (5.0, 2.0, -20.0, -20.0),
    (1e6, 0.5, -6.0, -300.0),
    (1.0, 0.5, -6.0, -6.0),  # the cutoff at the bandwidth level
]


def envelope_rate(fc: float, bw: float, bwr: float) -> mp.mpf:
    """a, where the envelope is exp(-a·t²)"""
    return 5 * (mp.pi * fc * bw) ** 2 / (-bwr * mp.log(10))


def exact_cutoff(fc: float, bw: float, bwr: float, tpr: float) -> mp.mpf:
    """the time where the envelope falls to tpr dB"""
    return mp.log(10) * mp.sqrt(mp.mpf(tpr) * bwr) / (10 * mp.pi * fc * bw)


def decibels(db: float) -> mp.mpf:
    """db as an amplitude ratio"""
    return mp.mpf(10) ** (mp.mpf(db) / 20)


def check_bandwidth(fc: float, bw: float, bwr: float) -> None:
    """asserts the spectrum exp(-pi²·f²/a) about fc is bwr dB down at fc·bw/2 from it"""
    a = envelope_rate(fc, bw, bwr)
    level = mp.exp(-((mp.pi * fc * bw / 2) ** 2) / a)
    assert abs(level - decibels(bwr)) <= SPEC_TOLERANCE * decibels(bwr), (
        fc,
        bw,
        bwr,
        level,
    )


def check_cutoff(fc: float, bw: float, bwr: float, tpr: float) -> None:
    """asserts the envelope is tpr dB down at the cutoff"""
    a = envelope_rate(fc, bw, bwr)
    tc = exact_cutoff(fc, bw, bwr, tpr)
    level = mp.exp(-a * tc**2)
    assert abs(level - decibels(tpr)) <= SPEC_TOLERANCE * decibels(tpr), (
        fc,
        bw,
        bwr,
        tpr,
        level,
    )


def exact_pulse(t, fc: float, bw: float, bwr: float) -> dict:
    """the in-phase part, quadrature part and envelope at each time in t"""
    a = envelope_rate(fc, bw, bwr)
    in_phase, quadrature, envelope = [], [], []
    for ti in t:
        ti = mp.mpf(float(ti))
        e = mp.exp(-a * ti**2)
        in_phase.append(float(e * mp.cos(2 * mp.pi * fc * ti)))
        quadrature.append(float(e * mp.sin(2 * mp.pi * fc * ti)))
        envelope.append(float(e))
    return {"in_phase": in_phase, "quadrature": quadrature, "envelope": envelope}


def check_scipy_pulse(t, fc: float, bw: float, bwr: float, output: dict) -> None:
    """asserts scipy's pulse is within SCIPY_ACCURACY of the exact pulse"""
    got = sc.signal.gausspulse(t, fc, bw, bwr, retquad=True, retenv=True)
    for name, y in zip(("in_phase", "quadrature", "envelope"), got, strict=True):
        worst = max(
            (abs(yi - xi) for yi, xi in zip(y, output[name], strict=True)),
            default=0,
        )
        assert worst <= SCIPY_ACCURACY, (fc, bw, bwr, name, worst)


def check_scipy_cutoff(fc: float, bw: float, bwr: float, tpr: float, tc: float) -> None:
    """asserts scipy's cutoff is within SCIPY_CUTOFF_ACCURACY of the exact cutoff"""
    got = sc.signal.gausspulse("cutoff", fc, bw, bwr, tpr)
    assert abs(got - tc) <= SCIPY_CUTOFF_ACCURACY * tc, (fc, bw, bwr, tpr, got, tc)


def generate_gausspulse() -> None:
    """the exact pulse on every grid, and the exact cutoff for every case"""
    with mp.workdps(DPS):
        for grid, (t, fc, bw, bwr) in GRIDS.items():
            check_bandwidth(fc, bw, bwr)
            output = exact_pulse(t, fc, bw, bwr)
            check_scipy_pulse(t, fc, bw, bwr, output)
            write_vector(
                "waveform/gausspulse",
                grid,
                f"exact, mpmath at {DPS} digits",
                {"t": t.tolist(), "fc": fc, "bw": bw, "bwr": bwr},
                output,
            )
        for fc, bw, bwr, tpr in CUTOFFS:
            check_bandwidth(fc, bw, bwr)
            check_cutoff(fc, bw, bwr, tpr)
            tc = float(exact_cutoff(fc, bw, bwr, tpr))
            check_scipy_cutoff(fc, bw, bwr, tpr, tc)
            write_vector(
                "waveform/gausspulse_cutoff",
                f"fc{fc:g}_bw{bw:g}_bwr{bwr:g}_tpr{tpr:g}",
                f"exact, mpmath at {DPS} digits",
                {"fc": fc, "bw": bw, "bwr": bwr, "tpr": tpr},
                {"tc": [tc]},
            )
