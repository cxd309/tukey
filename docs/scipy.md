<a href="../README.md"><img src="../images/logo.svg" alt="tukey" width="120" align="right"></a>

# Coming from SciPy

How to port existing applications from SciPy, or how to use tukey for engineers familiar with SciPy and NumPy.

## Packages

tukey is three packages, each usable on its own, mirroring SciPy's:

| Package        | For                                                                         | SciPy                                                         |
| -------------- | --------------------------------------------------------------------------- | ------------------------------------------------------------- |
| `tukey/fft`    | [Fourier transforms](#fft)                                                  | `scipy.fft`                                                   |
| `tukey/window` | [window functions](#window)                                                 | `scipy.signal.windows`                                        |
| `tukey/dsp`    | [filter design, filtering, signals, spectral analysis and reductions](#dsp) | `scipy.signal`, and NumPy's `sum`, `percentile`, `diff`, etc. |

## The basics

- Filter design takes a band, `dsp.Lowpass(Wn)`, `dsp.Highpass(Wn)`, `dsp.Bandpass(W1, W2)` or `dsp.Bandstop(W1, W2)`,
  rather than `btype`, and returns zeros, poles and gain (ZPK form), the same as `output='zpk'`.
  Call `BA()` for `(b, a)` or `SOS()` for sections.
- Frequencies for a design are fractions of Nyquist, the same as SciPy without `fs=`.
  `dsp.HzToNyquist(f, fs)` converts from Hz.
- Keyword arguments are passed as options, e.g. `padtype='even'` is `dsp.PaddingEven()` and `norm='ortho'` is `fft.NormOrtho()`.
- A function that can fail returns an error, where SciPy raises an exception.
- Functions take one-dimensional `[]float64`, so there is no `axis`, work a 2-D array row by row.

See [Accuracy](accuracy.md) for output accuracy comparison against each function.

## fft

| SciPy                                   | tukey                             |
| --------------------------------------- | --------------------------------- |
| `scipy.fft.fft(x)`                      | `y, err := fft.Forward(x)`        |
| `scipy.fft.fft(x, n)`                   | `fft.Forward(x, fft.N(n))`        |
| `scipy.fft.fft(x, norm='ortho')`        | `fft.Forward(x, fft.NormOrtho())` |
| `scipy.fft.ifft(X)`                     | `x, err := fft.Inverse(X)`        |
| `scipy.fft.rfft(x)`                     | `y, err := fft.ForwardReal(x)`    |
| `scipy.fft.irfft(X, n)`                 | `x, err := fft.InverseReal(X, n)` |
| `scipy.fft.fftfreq(n, d)`               | `f, err := fft.Freq(n, d)`        |
| `scipy.fft.rfftfreq(n, d)`              | `f, err := fft.RealFreq(n, d)`    |
| `scipy.fft.fftshift(X)`                 | `fft.Shift(X)`                    |
| `scipy.fft.ifftshift(X)`                | `fft.InverseShift(X)`             |
| `scipy.fft.next_fast_len(n, real=True)` | `fft.FastLength(n)`               |

### Differences to note:

- An empty input is an error, the same as SciPy.
- `fft.FastLength` allows factors of 2, 3 and 5, the same as `next_fast_len(n, real=True)`,
  where SciPy's default also allows 7 and 11.

## window

Every `scipy.signal.windows` window has the same name in Go's style, e.g. `blackmanharris` is `window.BlackmanHarris`.
`sym=False` is `window.Periodic()`.

| SciPy                                    | tukey                                                                  |
| ---------------------------------------- | ---------------------------------------------------------------------- |
| `hann(n)`                                | `w, err := window.Hann(n)`                                             |
| `hann(n, sym=False)`                     | `window.Hann(n, window.Periodic())`                                    |
| `kaiser(n, beta)`                        | `window.Kaiser(n, beta)`                                               |
| `taylor(n, nbar, sll, norm=False)`       | `window.Taylor(n, nbar, sll, window.Unnormalised())`                   |
| `exponential(n, center, tau, sym=False)` | `window.Exponential(n, tau, window.Periodic(), window.Center(center))` |
| `get_window('hann', n)`                  | `window.Hann(n)`, calling the window directly                          |

### Differences to note:

- Symmetric windows are exactly symmetric, where SciPy's are to within a few ULP.
- `dpss` isn't supported.

## dsp

### Filter design

| SciPy                                   | tukey                                                                           |
| --------------------------------------- | ------------------------------------------------------------------------------- |
| `butter(n, Wn, output='zpk')`           | `f, err := dsp.Butter(n, dsp.Lowpass(Wn))`, then `f.Zeros`, `f.Poles`, `f.Gain` |
| `butter(n, Wn)`                         | `f, err := dsp.Butter(n, dsp.Lowpass(Wn))`, then `b, a, err := f.BA()`          |
| `butter(n, Wn, output='sos')`           | `f, err := dsp.Butter(n, dsp.Lowpass(Wn))`, then `sos, err := f.SOS()`          |
| `butter(n, Wn, btype='highpass')`       | `dsp.Butter(n, dsp.Highpass(Wn))`                                               |
| `butter(n, [W1, W2], btype='bandpass')` | `dsp.Butter(n, dsp.Bandpass(W1, W2))`                                           |
| `butter(n, [W1, W2], btype='bandstop')` | `dsp.Butter(n, dsp.Bandstop(W1, W2))`                                           |
| `butter(n, f, fs=fs)`                   | `wn, err := dsp.HzToNyquist(f, fs)`, then `dsp.Butter(n, dsp.Lowpass(wn))`      |
| `cheby1(n, rp, Wn)`                     | `dsp.Cheby1(n, rp, dsp.Lowpass(Wn))`                                            |
| `cheby2(n, rs, Wn)`                     | `dsp.Cheby2(n, rs, dsp.Lowpass(Wn))`                                            |
| `ellip(n, rp, rs, Wn)`                  | `dsp.Ellip(n, rp, rs, dsp.Lowpass(Wn))`                                         |
| `bessel(n, Wn)`                         | `dsp.Bessel(n, dsp.Lowpass(Wn), dsp.BesselNormPhase())`                         |
| `bessel(n, Wn, norm='mag')`             | `dsp.Bessel(n, dsp.Lowpass(Wn))`                                                |
| `buttord(wp, ws, gpass, gstop)`         | `n, wn, err := dsp.ButterOrd(dsp.Lowpass(wp), dsp.Lowpass(ws), gpass, gstop)`   |
| `cheb1ord`, `cheb2ord`, `ellipord`      | `dsp.Cheby1Ord`, `dsp.Cheby2Ord`, `dsp.EllipOrd`, the same way                  |
| `iirnotch(w0, Q)`                       | `f, err := dsp.IIRNotch(w0, Q)`                                                 |
| `iirpeak(w0, Q)`                        | `f, err := dsp.IIRPeak(w0, Q)`                                                  |
| `firwin(n, Wn, window='hamming')`       | `w, err := window.Hamming(n)`, then `h, err := dsp.FIRWin(w, dsp.Lowpass(Wn))`  |
| `firwin(n, Wn, pass_zero=False)`        | `dsp.FIRWin(w, dsp.Highpass(Wn))`                                               |
| `firwin(n, Wn, scale=False)`            | `dsp.FIRWin(w, dsp.Lowpass(Wn), dsp.FIRNoScale())`                              |
| `zpk2tf(z, p, k)`                       | `b, a, err := f.BA()`                                                           |
| `zpk2sos(z, p, k)`                      | `sos, err := f.SOS()`                                                           |

#### Differences to note:

- `Bessel`'s band edge is the -3 dB point by default, the same as `Butter`'s,
  where SciPy's default is `norm='phase'`, pass `BesselNormPhase()` as an option to match this.
- `iirnotch` and `iirpeak` return `(b, a)`, where tukey returns zeros, poles and gain, so call `BA()`.
  `IIRPeak`'s gain is exact, where `iirpeak` is 6.1e-13 off for a narrow bandwidth.
- `FIRWin` takes the window itself, whose length is the number of taps, where `firwin` takes the number and a window name.
  A highpass or bandstop needs an odd number of taps.
- The order selections' bandstop edges are exact, where SciPy's come from a numerical minimiser, good to about 1e-5.
- Analog designs, `analog=True`, aren't supported.

### Filtering and response

| SciPy                               | tukey                                                              |
| ----------------------------------- | ------------------------------------------------------------------ |
| `lfilter(b, a, x)`                  | `y, err := dsp.Filter(b, a, x)`                                    |
| `lfilter(b, a, x, zi=zi)`           | `y, zf, err := dsp.FilterState(b, a, x, zi)`                       |
| `lfilter_zi(b, a)`                  | `zi, err := dsp.FilterZi(b, a)`                                    |
| `sosfilt(sos, x)`                   | `y, err := dsp.SOSFilter(sos, x)`                                  |
| `sosfilt(sos, x, zi=zi)`            | `y, zf, err := dsp.SOSFilterState(sos, x, zi)`                     |
| `sosfilt_zi(sos)`                   | `zi, err := dsp.SOSFilterZi(sos)`                                  |
| `filtfilt(b, a, x)`                 | `y, err := dsp.FiltFilt(b, a, x)`                                  |
| `filtfilt(b, a, x, padtype='even')` | `dsp.FiltFilt(b, a, x, dsp.PaddingEven())`                         |
| `filtfilt(b, a, x, padlen=n)`       | `dsp.FiltFilt(b, a, x, dsp.PadLen(n))`                             |
| `sosfiltfilt(sos, x)`               | `y, err := dsp.SOSFiltFilt(sos, x)`                                |
| `freqz(b, a, worN=n)`               | `w, h, err := dsp.Freqz(b, a, n)`                                  |
| `freqz(b, a, worN=n, fs=fs)`        | `dsp.Freqz(b, a, n)`, then `f, err := dsp.RadPerSampleToHz(w, fs)` |
| `freqz_sos(sos, worN=n)`            | `w, h, err := dsp.SOSFreqz(sos, n)`                                |
| `group_delay((b, a), w=n)`          | `w, gd, err := dsp.GroupDelay(b, a, n)`                            |

#### Differences to note:

- `filtfilt`'s `method='gust'` isn't supported.
- `GroupDelay` is exact, including the limit at a zero on the unit circle, where `group_delay` gives 0 there.

### Signals

| SciPy or NumPy                                        | tukey                                                             |
| ----------------------------------------------------- | ----------------------------------------------------------------- |
| `detrend(x)`                                          | `y := dsp.Detrend(x, dsp.DetrendLinear())`                        |
| `detrend(x, type='constant')`                         | `y := dsp.Detrend(x, dsp.DetrendConstant())`                      |
| `interp1d(x, y, bounds_error=False)(xi)`              | `yi, err := dsp.Interp1(x, y, xi)`                                |
| `medfilt(x, n)`                                       | `y, err := dsp.MedFilt(x, n, dsp.MedFiltZeroPad())`               |
| `savgol_filter(x, n, polyorder)`                      | `y, err := dsp.SavGolFilter(x, n, polyorder)`                     |
| `savgol_filter(x, n, polyorder, deriv=d, delta=dx)`   | `dsp.SavGolFilter(x, n, polyorder, dsp.SavGolDeriv(d, dx))`       |
| `convolve(a, v)`                                      | `c, err := dsp.Convolve(a, v)`                                    |
| `convolve(a, v, mode='same')`                         | `dsp.Convolve(a, v, dsp.ConvolveSame())`                          |
| `convolve(a, v, method='direct')`                     | `dsp.Convolve(a, v, dsp.ConvolveDirect())`                        |
| `fftconvolve(a, v)`                                   | `dsp.Convolve(a, v, dsp.ConvolveFFT())`                           |
| `correlate(a, v)`, `correlation_lags(len(a), len(v))` | `c, err := dsp.Correlate(a, v)`, then `c.Values` and `c.Lags`     |
| `find_peaks(x)`                                       | `peaks, err := dsp.FindPeaks(x)`                                  |
| `find_peaks(x, height=h, threshold=t)`                | `dsp.FindPeaks(x, dsp.MinPeakHeight(h), dsp.PeakThreshold(t))`    |
| `find_peaks(x, distance=d)`                           | `dsp.FindPeaks(x, dsp.MinPeakDistance(d))`                        |
| `find_peaks(x, prominence=p, width=w)`                | `dsp.FindPeaks(x, dsp.MinPeakProminence(p), dsp.MinPeakWidth(w))` |
| `peak_prominences(x, peaks)`                          | `dsp.PeakProminences(x, peaks)`                                   |
| `peak_widths(x, peaks, rel_height=r)`                 | `dsp.PeakWidths(x, peaks, r)`                                     |
| `hilbert(x)`                                          | `y := dsp.Hilbert(x)`                                             |
| `envelope(x)`                                         | `up, lo, err := dsp.Envelope(x, dsp.EnvelopeAnalytic())`          |

#### Differences to note:

- `MedFilt` truncates its window at the ends by default, where `medfilt` pads with zeros, `MedFiltZeroPad()`.
  A window holding NaN gives NaN, where `medfilt`'s result depends on how NaN sorts.
- `SavGolFilter` is exact at any degree, where `savgol_filter` loses digits for long windows and high degrees.
  Only `mode='interp'`, SciPy's default, is supported.
- `Interp1` gives NaN outside x's range, the same as `interp1d(bounds_error=False)`, where `numpy.interp` repeats the end values.
- `Correlate` returns the lags with the values.
  For `mode='same'` with an odd `len(a)` and an even `len(v)`, `correlation_lags` labels every lag one high, and tukey doesn't.
- `FindPeaks` takes a minimum for each limit, where `find_peaks` also takes a `(min, max)` range.
  It returns each peak's index and height.
  `wlen`, `rel_height` and `plateau_size` aren't supported.
- `Envelope` keeps the Nyquist bin, which `scipy.signal.envelope` drops, and adds an RMS envelope, `EnvelopeRMS(n)`.

### Spectral analysis

| SciPy                                                                               | tukey                                                            |
| ----------------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| `stft(x, fs, window=w, nperseg=len(w), boundary=None, padded=False, detrend=False)` | `tf, err := dsp.STFT(x, fs, w)`, which isn't divided by `sum(w)` |
| `spectrogram(x, fs, window=w, nperseg=len(w), detrend=False)`                       | `tf, err := dsp.Spectrogram(x, fs, w)`                           |
| `welch(x, fs, window=w, nperseg=len(w), detrend=False)`                             | `psd, err := dsp.Welch(x, fs, w)`                                |
| `welch(x, fs, ...)`, with its default `detrend='constant'`                          | `dsp.Welch(x, fs, w, dsp.SegmentDetrend(dsp.DetrendConstant()))` |
| `periodogram(x, fs, window=w, detrend=False)`                                       | `psd, err := dsp.Periodogram(x, fs, w)`                          |
| `noverlap=n`, `nfft=n`                                                              | `dsp.Overlap(n)`, `dsp.NFFT(n)`                                  |
| `scaling='spectrum'`                                                                | `dsp.ScalingSpectrum()`                                          |
| `average='median'`                                                                  | `dsp.AverageMedian()`                                            |
| `return_onesided=False`                                                             | `dsp.TwoSided()`                                                 |

#### Differences to note:

- There's no default window or segment length, where SciPy defaults to a 256-sample Hann window.
  Pass the window, its length is the segment length.
- There's no detrending by default, where `welch`, `spectrogram` and `periodogram` subtract each segment's mean.
- `STFT` isn't divided by the window's sum, so multiply SciPy's `stft` by `sum(w)` to compare.
  It doesn't pad the signal, so a partial last segment is dropped.
- `Values` are `[segment][frequency]`, the transpose of SciPy's.

### Reductions and calculus

| NumPy or SciPy                              | tukey                                         |
| ------------------------------------------- | --------------------------------------------- |
| `numpy.sum(x)`                              | `s := dsp.Sum(x)`                             |
| `numpy.mean(x)`                             | `m, err := dsp.Mean(x)`                       |
| `numpy.sqrt(numpy.mean(x**2))`              | `r, err := dsp.RMS(x)`                        |
| `numpy.percentile(x, p)`                    | `v, err := dsp.Percentile(x, p)`              |
| `numpy.percentile(x, p, method='hazen')`    | `dsp.Percentile(x, p, dsp.PercentileHazen())` |
| `numpy.diff(x)`                             | `d := dsp.Diff(x)`                            |
| `numpy.gradient(x, dx)`                     | `g, err := dsp.Grad(x, dx)`                   |
| `numpy.cumsum(x)`                           | `s := dsp.CumSum(x)`                          |
| `trapezoid(y, dx=dx)`                       | `q, err := dsp.Trapz(y, dx)`                  |
| `cumulative_trapezoid(y, dx=dx, initial=0)` | `q, err := dsp.CumTrapz(y, dx)`               |

#### Differences to note:

- Sums are compensated, so they're exact on cases where NumPy's pairwise sum rounds.
- `Mean`, `RMS` and `Percentile` of an empty signal return an error.
- `CumTrapz` starts at 0, the same as `cumulative_trapezoid(..., initial=0)`.
