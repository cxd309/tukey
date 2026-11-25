"""
scipy.signal.sawtooth, with exact values from mpmath

sawtooth(t, width) rises from -1 to 1 over the first width of each period of 2·pi,
then falls back to -1 over the rest.
As for square, the period is the float64 2·pi, which scipy reduces t by,
so t = 2·pi·k built with the float64 pi starts a cycle.
The reduction r = t mod 2·pi is exact in mpmath, and so is the wave from it.

scipy is within SCIPY_ACCURACY of the exact wave, except at width 1 for a tiny
negative t, where t % 2·pi rounds up to 2·pi and scipy divides 0 by 0, giving NaN
where the exact wave is just below 1.
"""

from typing import Final

import mpmath as mp
import numpy as np
import scipy as sc

from write import write_vector

DPS: Final = 50  # working precision, in significant digits
SCIPY_ACCURACY: Final = 1e-14  # how far scipy may be from exact, absolute
PERIOD: Final = 2 * np.pi  # the float64 2·pi, exact as an mpf

# the sample times
# steps: hundredths of a period from -2 to 8 periods, many on or next to a corner
# pi_multiples: k·pi from -20·pi to 20·pi, every one a corner at width 0.5
# large: around 10^6, where t has fewer digits left for the fraction of a cycle
# random: uniform from -50 to 50
# near_zero: either side of 0, where a tiny negative t reduces to just below 2·pi
# empty: no samples
GRIDS: Final = {
    "steps": np.arange(-200, 801) * (PERIOD / 100),
    "pi_multiples": np.arange(-20, 21) * np.pi,
    "large": 1e6 + np.arange(1001) * 0.01,
    "random": np.random.default_rng(seed=0).uniform(-50, 50, 1001),
    "near_zero": np.array(
        [-1e-300, -1e-20, -1e-16, -1e-8, -0.0, 0.0, 1e-300, 1e-20, 1e-16, 1e-8]
    ),
    "empty": np.array([]),
}

# 1 is a rising ramp, 0 a falling ramp and 0.5 a triangle
DUTIES: Final = (1.0, 0.5, 0.25, 0.9, 0.0)


def exact_sawtooth(t: float, duty: float) -> mp.mpf:
    """sawtooth at t, rising for the first duty of each float64 2·pi"""
    with mp.workdps(DPS):
        period = mp.mpf(PERIOD)
        t = mp.mpf(t)
        # r = t - k·period exactly, as a tiny negative t leaves r further below 2·pi than DPS digits reach
        k = mp.floor(t / period)
        r = mp.fsub(t, mp.fmul(k, period, exact=True), exact=True)
        # t/period rounded to DPS digits may have floored to the next whole period either way
        if r < 0:
            r = mp.fadd(r, period, exact=True)
        elif r >= period:
            r = mp.fsub(r, period, exact=True)
        if r < duty * period:
            return 2 * r / (duty * period) - 1
        return 1 - 2 * (r - duty * period) / ((1 - duty) * period)


def is_scipy_nan(t: float, duty: float) -> bool:
    """whether scipy gives NaN at t, which it does at width 1 where t % 2·pi rounds up to 2·pi"""
    return duty == 1 and float(np.float64(t) % PERIOD) == PERIOD


def check_scipy(t, duty: float, x) -> None:
    """asserts scipy's sawtooth is within SCIPY_ACCURACY of x, the exact wave"""
    y = sc.signal.sawtooth(t, duty)
    for ti, xi, yi in zip(t, x, y, strict=True):
        if is_scipy_nan(ti, duty):
            assert np.isnan(yi), (ti, duty, yi)
            continue
        assert abs(yi - xi) <= SCIPY_ACCURACY, (ti, duty, xi, yi)


def generate_sawtooth() -> None:
    """the exact sawtooth on every grid at every duty"""
    for grid, t in GRIDS.items():
        for duty in DUTIES:
            x = [float(exact_sawtooth(float(ti), duty)) for ti in t]
            check_scipy(t, duty, x)
            write_vector(
                "waveform/sawtooth",
                f"{grid}_duty{duty:g}",
                f"exact, mpmath at {DPS} digits",
                {"t": t.tolist(), "duty": duty},
                {"x": x},
            )
