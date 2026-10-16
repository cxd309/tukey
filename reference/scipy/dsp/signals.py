"""Synthetic input signals shared by the filtering vectors"""

import numpy as np
import scipy as sc


def synthetic_signals(n: int) -> dict[str, np.ndarray]:
    t = np.arange(n) / n
    rng = np.random.default_rng(seed=0)  # fixed seed: reproducible noise
    impulse = np.zeros(n)
    impulse[0] = 1.0
    return {
        "impulse": impulse,  # output is the impulse response itself
        "step": np.r_[np.zeros(n // 2), np.ones(n - n // 2)],  # mid-signal edge
        "constant": np.full(n, 3.0),  # DC: filtfilt must preserve exactly
        "sines": np.sin(2 * np.pi * 5 * t) + 0.5 * np.sin(2 * np.pi * 40 * t),
        "chirp": sc.signal.chirp(t, f0=1, t1=1, f1=80),  # sweeps through cutoff
        "noise": rng.standard_normal(n),  # broadband
    }
