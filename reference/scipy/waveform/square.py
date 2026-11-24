"""
scipy.signal.square, checked against the exact square wave away from its edges

square(t, duty) is 1 where t mod 2·pi is below duty·2·pi, and -1 elsewhere.
scipy reduces t by the float64 2·pi, which is 2.4e-16 short of the true 2·pi,
where the exact square wave reduces it by the true 2·pi.
They differ only at samples within about t·4e-17 of an edge, and there scipy's
is the one a caller means: t = 2·pi·k built with the float64 pi is the start of
a cycle in scipy, and just before its end in the exact wave.
So the vectors are scipy's own output, and tukey follows it.

Each sample is checked against the exact square wave, computed with mpmath,
wherever it's further than EDGE_DISTANCE cycles from an edge.
"""

from typing import Final

import mpmath as mp
import numpy as np
import scipy as sc

from write import write_vector

DPS: Final = 50  # working precision, in significant digits
EDGE_DISTANCE: Final = mp.mpf(
    "1e-9"
)  # how close to an edge, in cycles, scipy may differ

# the sample times
# steps: hundredths of a period from -2 to 8 periods, many on or next to an edge
# pi_multiples: k·pi from -20·pi to 20·pi, every one an edge at duty 0.5
# large: around 10^6, where t has fewer digits left for the fraction of a cycle
# random: uniform from -50 to 50
# empty: no samples
GRIDS: Final = {
    "steps": np.arange(-200, 801) * (2 * np.pi / 100),
    "pi_multiples": np.arange(-20, 21) * np.pi,
    "large": 1e6 + np.arange(1001) * 0.01,
    "random": np.random.default_rng(seed=0).uniform(-50, 50, 1001),
    "empty": np.array([]),
}

# 0 is -1 everywhere and 1 is 1 everywhere
DUTIES: Final = (0.5, 0.25, 0.9, 0.0, 1.0)


def check_away_from_edges(t, duty: float, x) -> None:
    """asserts x is the exact square wave at every t further than EDGE_DISTANCE from an edge"""
    with mp.workdps(DPS):
        for ti, xi in zip(t, x, strict=True):
            fraction = mp.frac(mp.mpf(float(ti)) / (2 * mp.pi))
            distance = min(fraction, 1 - fraction, abs(fraction - duty))
            if distance > EDGE_DISTANCE:
                want = 1.0 if fraction < duty else -1.0
                assert xi == want, (ti, duty, xi)


def generate_square() -> None:
    """scipy.signal.square on every grid at every duty"""
    for grid, t in GRIDS.items():
        for duty in DUTIES:
            x = sc.signal.square(t, duty)
            check_away_from_edges(t, duty, x)
            write_vector(
                "waveform/square",
                f"{grid}_duty{duty:g}",
                f"scipy.signal.square(t, {duty:g})",
                {"t": t.tolist(), "duty": duty},
                {"x": x.tolist()},
            )
