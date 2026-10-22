"""
Spectral analysis: the STFT, at FFT bins and at chosen frequencies, the
spectrogram, Welch's method and the periodogram

Every case passes its window as an array, so the segment length is the window's
and both sources transform exactly the same thing, with per-segment detrending
off unless a case is testing it. Time-frequency outputs are written
[segment][frequency], tukey's layout.
"""

import math
from fractions import Fraction
from typing import Final

import numpy as np
import scipy as sc

from write import complex_json, write_vector

FS: Final = 1000.0  # Hz, for every case but the band zoom
LENGTH: Final = 4096


def spectral_signals() -> dict[str, np.ndarray]:
    """two tones in noise, and a chirp, LENGTH samples at FS"""
    t = np.arange(LENGTH) / FS
    rng = np.random.default_rng(seed=11)
    tones = (
        np.sin(2 * np.pi * 50 * t)
        + 0.5 * np.sin(2 * np.pi * 123.4 * t)  # between bins at every length used
        + 0.1 * rng.standard_normal(LENGTH)
    )
    chirp = sc.signal.chirp(t, f0=10, t1=t[-1], f1=400)
    return {"tones": tones, "chirp": chirp}


def spectral_windows() -> dict[str, np.ndarray]:
    """a periodic window, a symmetric one of odd length, and a rectangular one"""
    return {
        "hann256": sc.signal.windows.hann(256, sym=False),
        "hamming255": sc.signal.windows.hamming(255),
        "boxcar128": sc.signal.windows.boxcar(128),
    }


def segment_starts(n: int, length: int, overlap: int) -> range:
    """the start of every whole segment: a partial last one is dropped"""
    return range(0, n - length + 1, length - overlap)


def plain_stft(
    x: np.ndarray, window: np.ndarray, overlap: int, nfft: int, two_sided: bool
) -> np.ndarray:
    """
    the textbook STFT, [segment][frequency]:
    each windowed segment's FFT, with its phase measured from the segment's start and no scaling
    """
    length = len(window)
    transform = np.fft.fft if two_sided else np.fft.rfft
    return np.array(
        [
            transform(window * x[s : s + length], nfft)
            for s in segment_starts(len(x), length, overlap)
        ]
    )


def twiddles(f: float, fs: float, length: int) -> tuple[np.ndarray, np.ndarray]:
    """
    cos and sin of 2πfn/fs for n = 0 to length-1, with fn/fs reduced to a
    fraction of a turn exactly, in rationals, before it's rounded. Computing
    2πfn/fs in floats instead rounds an angle of up to hundreds of radians,
    which put the sums 1e-14 from the true values
    """
    turns = [float((Fraction(f) * n / Fraction(fs)) % 1) for n in range(length)]
    phase = 2 * np.pi * np.array(turns)
    return np.cos(phase), np.sin(phase)


def exact_dft(segment: np.ndarray, cos: np.ndarray, sin: np.ndarray) -> complex:
    """
    sum w[n]·x[n]·e^(-j2πfn/fs) for one windowed segment, given its twiddles,
    each part summed exactly with math.fsum
    """
    re = math.fsum((segment * cos).tolist())
    im = -math.fsum((segment * sin).tolist())
    return complex(re, im)


def generate_stft() -> None:
    """scipy.signal.stft without padding, scaled back to the plain transform"""
    signals, windows = spectral_signals(), spectral_windows()
    tones, chirp = signals["tones"], signals["chirp"]
    cases = [
        # name, x, window, overlap, nfft (None for the window's length), two-sided
        ("tones_hann256", tones, "hann256", 128, None, False),
        ("tones_hann256_overlap0", tones, "hann256", 0, None, False),
        ("chirp_hamming255_overlap200", chirp, "hamming255", 200, None, False),
        ("tones_hann256_nfft512", tones, "hann256", 128, 512, False),
        ("tones_hamming255_nfft301", tones, "hamming255", 127, 301, False),  # odd
        ("tones_boxcar128", tones, "boxcar128", 64, None, False),
        (
            "partial_hann256",
            tones[:4000],
            "hann256",
            128,
            None,
            False,
        ),  # a partial last segment
        ("one_segment", tones[:256], "hann256", 128, None, False),
        ("tones_hann256_twosided", tones, "hann256", 128, None, True),
    ]
    for name, x, wname, overlap, nfft, two_sided in cases:
        w = windows[wname]
        f, t, z = sc.signal.stft(
            x,
            fs=FS,
            window=w,
            nperseg=len(w),
            noverlap=overlap,
            nfft=nfft,
            detrend=False,
            return_onesided=not two_sided,
            boundary=None,
            padded=False,
            scaling="spectrum",
        )
        values = z.T * w.sum()  # 'spectrum' divides by the window's sum
        plain = plain_stft(x, w, overlap, nfft or len(w), two_sided)
        assert np.allclose(values, plain, rtol=1e-12, atol=1e-12 * abs(plain).max()), (
            f"{name}: scipy's stft isn't the plain transform scaled by 1/sum(w)"
        )
        write_vector(
            "stft",
            name,
            "scipy.signal.stft(boundary=None, padded=False, detrend=False) * sum(window)",
            {
                "x": x.tolist(),
                "fs": FS,
                "window": w.tolist(),
                "overlap": overlap,
                "nfft": nfft,
                "two_sided": two_sided,
            },
            {
                "frequencies": f.tolist(),
                "times": t.tolist(),
                "values": complex_json(values),
            },
        )


def generate_stft_frequencies() -> None:
    """the STFT at chosen frequencies, each summed exactly"""
    signals, windows = spectral_signals(), spectral_windows()
    tones = signals["tones"]
    # a band of interest at fine resolution: 0.1 Hz steps, where a 1250-sample
    # window at 250 Hz has FFT bins every 0.2 Hz
    zoom_fs = 250.0
    zoom_t = np.arange(5000) / zoom_fs
    zoom_rng = np.random.default_rng(seed=12)
    zoom_x = (
        np.sin(2 * np.pi * 3.2 * zoom_t)
        + 0.3 * np.sin(2 * np.pi * 7.75 * zoom_t)
        + 0.1 * zoom_rng.standard_normal(len(zoom_t))
    )
    cases = [
        # name, x, fs, window, overlap, frequencies
        (
            "tones_between_bins",
            tones,
            FS,
            windows["hann256"],
            128,
            np.array([0.0, 12.3, 50.0, 123.4, 499.9, 500.0]),  # DC and Nyquist included
        ),
        (
            "tones_fft_bins",  # exactly the FFT's bins: must agree with stft
            tones,
            FS,
            windows["hann256"],
            128,
            np.arange(129) * FS / 256,
        ),
        (
            "band_zoom",
            zoom_x,
            zoom_fs,
            sc.signal.windows.hamming(1250),
            250,
            np.arange(5, 201) / 10,  # 0.5 to 20 Hz
        ),
    ]
    for name, x, fs, w, overlap, freqs in cases:
        length = len(w)
        starts = segment_starts(len(x), length, overlap)
        turns = [twiddles(f, fs, length) for f in freqs]  # the same for every segment
        values = np.array(
            [
                [exact_dft(w * x[s : s + length], cos, sin) for cos, sin in turns]
                for s in starts
            ]
        )
        times = np.array([(s + length / 2) / fs for s in starts])  # as scipy's stft
        if name == "tones_fft_bins":
            plain = plain_stft(x, w, overlap, length, False)
            assert np.allclose(
                values, plain, rtol=1e-12, atol=1e-12 * abs(plain).max()
            ), "the exact sums disagree with the FFT at its own bins"
        write_vector(
            "stft_frequencies",
            name,
            "sum w[n]·x[n]·exp(-j2πfn/fs) per segment, summed exactly with math.fsum",
            {
                "x": x.tolist(),
                "fs": fs,
                "window": w.tolist(),
                "overlap": overlap,
                "frequencies": freqs.tolist(),
            },
            {
                "frequencies": freqs.tolist(),
                "times": times.tolist(),
                "values": complex_json(values),
            },
        )


def generate_spectrogram() -> None:
    """scipy.signal.spectrogram's power, density and spectrum scaled"""
    signals, windows = spectral_signals(), spectral_windows()
    tones, chirp = signals["tones"], signals["chirp"]
    cases = [
        # name, x, window, overlap, nfft, scaling
        ("tones_hann256_density", tones, "hann256", 128, None, "density"),
        ("tones_hann256_spectrum", tones, "hann256", 128, None, "spectrum"),
        ("chirp_hamming255_nfft301", chirp, "hamming255", 200, 301, "density"),
        ("tones_boxcar128_overlap0", tones, "boxcar128", 0, None, "spectrum"),
    ]
    for name, x, wname, overlap, nfft, scaling in cases:
        w = windows[wname]
        f, t, sxx = sc.signal.spectrogram(
            x,
            fs=FS,
            window=w,
            nperseg=len(w),
            noverlap=overlap,
            nfft=nfft,
            detrend=False,
            scaling=scaling,
            mode="psd",
        )
        write_vector(
            "spectrogram",
            name,
            f"scipy.signal.spectrogram(detrend=False, scaling={scaling!r}, mode='psd')",
            {
                "x": x.tolist(),
                "fs": FS,
                "window": w.tolist(),
                "overlap": overlap,
                "nfft": nfft,
                "scaling": scaling,
            },
            {"frequencies": f.tolist(), "times": t.tolist(), "values": sxx.T.tolist()},
        )


def generate_welch() -> None:
    """scipy.signal.welch: scalings, averages, and per-segment detrending"""
    signals, windows = spectral_signals(), spectral_windows()
    tones, chirp = signals["tones"], signals["chirp"]
    offset = tones + 5.0  # a large mean, so detrending each segment matters
    cases = [
        # name, x, window, overlap, nfft, scaling, average, detrend
        ("tones_hann256", tones, "hann256", 128, None, "density", "mean", "none"),
        (
            "tones_hann256_spectrum",
            tones,
            "hann256",
            128,
            None,
            "spectrum",
            "mean",
            "none",
        ),
        (
            "tones_hann256_median",
            tones,
            "hann256",
            128,
            None,
            "density",
            "median",
            "none",
        ),
        (
            "chirp_hamming255_nfft512",
            chirp,
            "hamming255",
            200,
            512,
            "density",
            "mean",
            "none",
        ),
        (
            "partial_boxcar128",
            tones[:4000],
            "boxcar128",
            64,
            None,
            "density",
            "mean",
            "none",
        ),
        ("offset_hann256", offset, "hann256", 128, None, "density", "mean", "none"),
        (
            "offset_hann256_constant",
            offset,
            "hann256",
            128,
            None,
            "density",
            "mean",
            "constant",
        ),
        (
            "offset_hann256_linear",
            offset,
            "hann256",
            128,
            None,
            "density",
            "mean",
            "linear",
        ),
    ]
    for name, x, wname, overlap, nfft, scaling, average, detrend in cases:
        w = windows[wname]
        f, pxx = sc.signal.welch(
            x,
            fs=FS,
            window=w,
            nperseg=len(w),
            noverlap=overlap,
            nfft=nfft,
            detrend=False if detrend == "none" else detrend,
            scaling=scaling,
            average=average,
        )
        write_vector(
            "welch",
            name,
            f"scipy.signal.welch(detrend={detrend!r}, scaling={scaling!r}, average={average!r})",
            {
                "x": x.tolist(),
                "fs": FS,
                "window": w.tolist(),
                "overlap": overlap,
                "nfft": nfft,
                "scaling": scaling,
                "average": average,
                "detrend": detrend,
            },
            {"frequencies": f.tolist(), "power": pxx.tolist()},
        )


def generate_periodogram() -> None:
    """scipy.signal.periodogram of the whole signal"""
    tones = spectral_signals()["tones"]
    cases = [
        # name, window, nfft, scaling
        ("tones_boxcar", sc.signal.windows.boxcar(LENGTH), None, "density"),
        (
            "tones_hann_spectrum",
            sc.signal.windows.hann(LENGTH, sym=False),
            None,
            "spectrum",
        ),
        ("tones_boxcar_nfft8192", sc.signal.windows.boxcar(LENGTH), 8192, "density"),
    ]
    for name, w, nfft, scaling in cases:
        f, pxx = sc.signal.periodogram(
            tones, fs=FS, window=w, nfft=nfft, detrend=False, scaling=scaling
        )
        write_vector(
            "periodogram",
            name,
            f"scipy.signal.periodogram(detrend=False, scaling={scaling!r})",
            {
                "x": tones.tolist(),
                "fs": FS,
                "window": w.tolist(),
                "nfft": nfft,
                "scaling": scaling,
            },
            {"frequencies": f.tolist(), "power": pxx.tolist()},
        )
