<a href="../README.md"><img src="../images/logo.svg" alt="tukey" width="120" align="right"></a>

# Coming from MATLAB

How to port existing applications from MATLAB or how to use tukey for engineers familiar with matlab. Some of these functions are inside the DSP System Toolbox or Signal Processing Toolbox add-ons.

## Packages

tukey is three packages, each usable on its own:

| Package        | For                                                                         | MATLAB                                            |
| -------------- | --------------------------------------------------------------------------- | ------------------------------------------------- |
| `tukey/fft`    | [Fourier transforms](#fft)                                                  | `fft`, `ifft` and their helpers                   |
| `tukey/window` | [window functions](#window)                                                 | `hann`, `kaiser` and the other windows            |
| `tukey/dsp`    | [filter design, filtering, signals, spectral analysis and reductions](#dsp) | `butter`, `filtfilt`, `findpeaks`, `pwelch`, etc. |

## The basics

- Filter design takes a band, `dsp.Lowpass(Wn)`, `dsp.Highpass(Wn)`, `dsp.Bandpass(W1, W2)` or `dsp.Bandstop(W1, W2)`,
  rather than a type string, and returns in zeros, poles and gain (ZPK) form.
  Call `BA()` for `[b, a]` form or `SOS()` for sections.
- Frequencies for a design are fractions of Nyquist, the same as MATLAB's `Wn`.
  `dsp.HzToNyquist(f, fs)` converts from Hz.
- Every function that can fail returns an error, (standard Go functionality).
- Indices are 0-based, so `findpeaks`' `locs` and `findsignal`'s `istart` are one less in tukey.
- Functions take vectors, as `[]float64`, not matrices, so work a matrix column by column.
- Where MATLAB and SciPy differ, tukey follows SciPy by default, and an option gives MATLAB's behaviour.
  The tables below include those options.

See [Accuracy](accuracy.md) for output accuracy comparison against each function.

## fft

| MATLAB                 | tukey                         |
| ---------------------- | ----------------------------- |
| `fft(x)`               | `fft.Forward(x)`              |
| `fft(x, n)`            | `fft.Forward(x, fft.N(n))`    |
| `ifft(X)`              | `fft.Inverse(X)`              |
| `fft(x)`, first half   | `fft.ForwardReal(x)`          |
| `ifft(X, 'symmetric')` | `fft.InverseReal(X, n)`       |
| `(0:n-1)*fs/n`         | `f, err := fft.Freq(n, 1/fs)` |
| `fftshift(X)`          | `fft.Shift(X)`                |
| `ifftshift(X)`         | `fft.InverseShift(X)`         |

### Differences to note:

- `fft.Freq` lists the upper half of the bins as negative frequencies, where `(0:n-1)*fs/n` lists them above Nyquist.
- `fft.InverseReal` needs its length, as n/2 + 1 bins come from both an even and an odd length.

## window

Most windows have the same name.
`window.Periodic()` option is MATLAB's `'periodic'`.

| MATLAB                    | tukey                                                 |
| ------------------------- | ----------------------------------------------------- |
| `hann(n, 'periodic')`     | `window.Hann(n, window.Periodic())`                   |
| `rectwin(n)`              | `window.Boxcar(n)`                                    |
| `nuttallwin(n)`           | `window.Nuttall(n)`                                   |
| `flattopwin(n)`           | `window.FlatTop(n)`                                   |
| `tukeywin(n, r)`          | `window.Tukey(n, r)`                                  |
| `gausswin(n, alpha)`      | `window.Gaussian(n, (n-1)/(2*alpha))`                 |
| `taylorwin(n, nbar, sll)` | `window.Taylor(n, nbar, -sll, window.Unnormalised())` |
| `hanning(n)`              | `window.Hann(n + 2)`, less its first and last sample  |

### Differences to note:

- `gausswin` takes the reciprocal of a width relative to half the window, where `window.Gaussian` takes the standard deviation in samples.
- `taylorwin` takes its sidelobe level as negative dB and doesn't normalise, so pass `-sll` and `window.Unnormalised()`.
- `hanning` is the Hann window without its zero ends.
- At n = 1 every window is `[1]`, where `taylorwin` evaluates its formula.

## dsp

### Filter design

| MATLAB                              | tukey                                                                                        |
| ----------------------------------- | -------------------------------------------------------------------------------------------- |
| `[b, a] = butter(n, Wn)`            | `f, err := dsp.Butter(n, dsp.Lowpass(Wn))`, then `b, a, err := f.BA()`                       |
| `butter(n, Wn, 'high')`             | `dsp.Butter(n, dsp.Highpass(Wn))`                                                            |
| `butter(n, [W1 W2])`                | `dsp.Butter(n, dsp.Bandpass(W1, W2))`                                                        |
| `butter(n, [W1 W2], 'stop')`        | `dsp.Butter(n, dsp.Bandstop(W1, W2))`                                                        |
| `[z, p, k] = butter(n, Wn)`         | `f, err := dsp.Butter(n, dsp.Lowpass(Wn))`, then `f.Zeros`, `f.Poles`, `f.Gain`              |
| `[sos, g] = zp2sos(z, p, k)`        | `sos, err := f.SOS()`, with `g` already in the first section                                 |
| `cheby1(n, Rp, Wn)`                 | `dsp.Cheby1(n, Rp, dsp.Lowpass(Wn))`                                                         |
| `cheby2(n, Rs, Wn)`                 | `dsp.Cheby2(n, Rs, dsp.Lowpass(Wn))`                                                         |
| `ellip(n, Rp, Rs, Wn)`              | `dsp.Ellip(n, Rp, Rs, dsp.Lowpass(Wn))`                                                      |
| `besself(n, Wo)`, then `bilinear`   | `dsp.Bessel(n, dsp.Lowpass(Wn), dsp.BesselNormPhase())`, digital directly                    |
| `[n, Wn] = buttord(Wp, Ws, Rp, Rs)` | `n, wn, err := dsp.ButterOrd(dsp.Lowpass(Wp), dsp.Lowpass(Ws), Rp, Rs, dsp.MatchStopband())` |
| `cheb1ord(Wp, Ws, Rp, Rs)`          | `dsp.Cheby1Ord(pass, stop, Rp, Rs)`                                                          |
| `cheb2ord(Wp, Ws, Rp, Rs)`          | `dsp.Cheby2Ord(pass, stop, Rp, Rs, dsp.MatchStopband())`                                     |
| `ellipord(Wp, Ws, Rp, Rs)`          | `dsp.EllipOrd(pass, stop, Rp, Rs)`                                                           |
| `[b, a] = iirnotch(w0, bw)`         | `f, err := dsp.IIRNotch(w0, w0/bw)`                                                          |
| `[b, a] = iirpeak(w0, bw)`          | `f, err := dsp.IIRPeak(w0, w0/bw)`                                                           |
| `b = fir1(n, Wn)`                   | `w, err := window.Hamming(n + 1)`, then `h, err := dsp.FIRWin(w, dsp.Lowpass(Wn))`           |
| `b = fir1(n, Wn, w)`                | `h, err := dsp.FIRWin(w, dsp.Lowpass(Wn))`                                                   |

#### Differences to note:

- `f.SOS()` pairs and orders sections the same as SciPy's `zpk2sos`, with the gain in the first section.
  MATLAB's `zp2sos` gives the same filter in different sections, so compare responses rather than coefficients.
- `Bessel`'s band edge is the -3 dB point by default, the same as `Butter`'s.
  MATLAB has no digital Bessel design, and `besself` is the phase normalisation, `BesselNormPhase()`.
- The order selections meet the passband edge by default.
  MATLAB's `buttord` and `cheb2ord` meet the stopband edge, so pass `MatchStopband()`.
  For a bandstop, MATLAB's `cheb1ord` and `ellipord` return the passband edges unmoved,
  which for `cheb1ord` can miss the stopband specification, where tukey moves them.
  For a bandpass, MATLAB's `cheb2ord` returns both stopband edges, where tukey meets the tighter one.
  The orders always match.
- `IIRNotch` and `IIRPeak` take the quality factor q, where MATLAB takes the bandwidth, w0/q.
- `FIRWin` takes the window itself, whose length is the number of taps, where `fir1` takes the order, one less.
  A highpass or bandstop needs an odd number of taps.

### Filtering and response

| MATLAB                          | tukey                                                                    |
| ------------------------------- | ------------------------------------------------------------------------ |
| `y = filter(b, a, x)`           | `y, err := dsp.Filter(b, a, x)`                                          |
| `[y, zf] = filter(b, a, x, zi)` | `y, zf, err := dsp.FilterState(b, a, x, zi)`                             |
| `y = sosfilt(sos, x)`           | `y, err := dsp.SOSFilter(sos, x)`                                        |
| `y = filtfilt(b, a, x)`         | `y, err := dsp.FiltFilt(b, a, x, dsp.MATLABPadLen())`                    |
| `y = filtfilt(sos, g, x)`       | `y, err := dsp.SOSFiltFilt(sos, x, dsp.BySection(), dsp.MATLABPadLen())` |
| `[h, w] = freqz(b, a, n)`       | `w, h, err := dsp.Freqz(b, a, n)`                                        |
| `[h, w] = freqz(sos, n)`        | `w, h, err := dsp.SOSFreqz(sos, n)`                                      |
| `[gd, w] = grpdelay(b, a, n)`   | `w, gd, err := dsp.GroupDelay(b, a, n)`                                  |
| `[gd, w] = grpdelay(sos, n)`    | `w, gd, err := dsp.SOSGroupDelay(sos, n)`                                |

#### Differences to note:

- `FiltFilt` extends the signal past each end before filtering, so the filter's start-up transients fall on the extension rather than the signal.

  By default the extension uses the same convention as SciPy, 3 × the number of coefficients, e.g. 15 samples for a 4th-order filter.

  MATLAB instead uses 3 × the filter's order, i.e. 12 samples for the same 4th-order files, to match this use the `MATLABPadLen()` option.

  The two results differ slightly near the ends, and agree elsewhere.
- `SOSFiltFilt` filters the whole cascade forward then backward, the same as SciPy.
  `BySection()` option filters each section in turn, the same as MATLAB.
  They give the same filter, but different edge transients, which for a slow filter reach well into the signal.
- `Freqz` and `GroupDelay` return the frequencies first, then the response.
  `dsp.RadPerSampleToHz(w, fs)` converts the frequencies to Hz.
- `Freqz` with n= 0 returns empty slices, where MATLAB errors.
- `GroupDelay` is exact near poles and zeros crowded against the unit circle, where `grpdelay` loses digits.

### Signals

| MATLAB                                             | tukey                                                                                            |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| `y = detrend(x)`                                   | `y := dsp.Detrend(x, dsp.DetrendLinear())`                                                       |
| `y = detrend(x, 0)`                                | `y := dsp.Detrend(x, dsp.DetrendConstant())`                                                     |
| `yi = interp1(x, y, xi)`                           | `yi, err := dsp.Interp1(x, y, xi)`                                                               |
| `fillmissing(x, 'linear')`                         | `y := dsp.FillMissing(x, dsp.FillLinear())`                                                      |
| `fillmissing(x, 'constant', v)`                    | `y := dsp.FillMissing(x, dsp.FillConstant(v))`                                                   |
| `y = medfilt1(x, n)`                               | `y, err := dsp.MedFilt(x, n, dsp.MedFiltZeroPad())`                                              |
| `y = medfilt1(x, n, 'truncate')`                   | `y, err := dsp.MedFilt(x, n)`                                                                    |
| `y = sgolayfilt(x, order, n)`                      | `y, err := dsp.SavGolFilter(x, n, order)`                                                        |
| `c = conv(a, v)`                                   | `c, err := dsp.Convolve(a, v)`                                                                   |
| `c = conv(a, v, 'same')`                           | `c, err := dsp.Convolve(a, v, dsp.ConvolveSame())`                                               |
| `[r, lags] = xcorr(a, v)`                          | `c, err := dsp.Correlate(a, v)`, then `c.Values` and `c.Lags`                                    |
| `[pks, locs] = findpeaks(x)`                       | `peaks, err := dsp.FindPeaks(x, dsp.MATLABPeaks())`                                              |
| `findpeaks(x, 'MinPeakDistance', d)`               | `dsp.FindPeaks(x, dsp.MATLABPeaks(), dsp.MinPeakDistance(d+1))`                                  |
| `findpeaks(x, 'MinPeakHeight', h)`                 | `dsp.FindPeaks(x, dsp.MATLABPeaks(), dsp.MinPeakHeight(h))`                                      |
| `findpeaks(x, 'Threshold', t)`                     | `dsp.FindPeaks(x, dsp.MATLABPeaks(), dsp.PeakThreshold(t))`                                      |
| `findpeaks(x, 'MinPeakProminence', p)`             | `dsp.FindPeaks(x, dsp.MATLABPeaks(), dsp.MinPeakProminence(p))`                                  |
| `findpeaks(x, 'MinPeakWidth', w)`                  | `dsp.FindPeaks(x, dsp.MATLABPeaks(), dsp.MinPeakWidth(w))`                                       |
| `[~, ~, w, p] = findpeaks(x)`                      | `dsp.PeakWidths(x, indices, 0.5)` and `dsp.PeakProminences(x, indices)`                          |
| `y = hilbert(x)`                                   | `y := dsp.Hilbert(x)`                                                                            |
| `[up, lo] = envelope(x)`                           | `up, lo, err := dsp.Envelope(x, dsp.EnvelopeAnalytic())`                                         |
| `[up, lo] = envelope(x, n, 'rms')`                 | `up, lo, err := dsp.Envelope(x, dsp.EnvelopeRMS(n))`                                             |
| `[istart, istop, dist] = findsignal(data, signal)` | `match, err := dsp.FindSignal(data, signal)`, then `match.Start`, `match.Stop`, `match.Distance` |

#### Differences to note:

- `ConvolveSame` starts one output earlier than `conv(a, v, 'same')` for an even-length `v`, the same as SciPy.
- `ConvolveValid` with `v` longer than `a` gives `len(v) - len(a) + 1` outputs, where `conv` gives none.
- `Correlate` gives only the lags where `a` and `v` overlap, where `xcorr` pads the shorter with zeros.
  `xcorr`'s normalisations aren't supported yet.
- `MedFilt` truncates its window at the ends by default, where `medfilt1` pads with zeros, `MedFiltZeroPad()`.
- `SavGolFilter` takes `(x, n, order)`, where `sgolayfilt` takes `(x, order, n)`,
  and `SavGolDeriv` gives derivatives, which `sgolayfilt` doesn't.
- `FindPeaks` follows SciPy's rules unless given `MATLABPeaks()`.
  MATLAB's `'MinPeakDistance', d` keeps peaks more than d apart, so it's `MinPeakDistance(d+1)`.
  `findpeaks` errors for fewer than 3 samples, where tukey returns no peaks.
  `'SortStr'`, `'NPeaks'` and `'WidthReference'` aren't supported yet.
- `Envelope` supports MATLAB's analytic and RMS envelopes, not `'peak'` or an FIR Hilbert filter.
- `FindSignal` gives the exact distance of the closest window, where `findsignal` rounds,
  so its perfect matches can be up to 2^-29 rather than 0.
  `SearchFast()` rounds the same way as MATLAB.
  Only `findsignal`'s default options are supported so far.
- `Interp1` needs `x` strictly increasing.

### Spectral analysis

| MATLAB                                              | tukey                                                                                                |
| --------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `[S, F, T] = spectrogram(x, w, noverlap, nfft, fs)` | `tf, err := dsp.STFT(x, fs, w, dsp.Overlap(noverlap), dsp.NFFT(nfft))`                               |
| `spectrogram(x, w, noverlap, f, fs)`                | `tf, err := dsp.STFT(x, fs, w, dsp.Overlap(noverlap), dsp.AtFrequencies(f))`                         |
| `[~, F, T, P] = spectrogram(..., 'power')`          | `tf, err := dsp.Spectrogram(x, fs, w, dsp.Overlap(noverlap), dsp.NFFT(nfft), dsp.ScalingSpectrum())` |
| `[pxx, f] = pwelch(x, w, noverlap, nfft, fs)`       | `psd, err := dsp.Welch(x, fs, w, dsp.Overlap(noverlap), dsp.NFFT(nfft))`                             |
| `[pxx, f] = periodogram(x, w, nfft, fs)`            | `psd, err := dsp.Periodogram(x, fs, w, dsp.NFFT(nfft))`                                              |

Where MATLAB is given an integer n for its window, pass `window.Hamming(n)`, which matches MATLAB exactly.

#### Differences to note:

- There's no default window or segment length, as MATLAB and SciPy each choose their own, and differently.
  Pass the window, and its length is the segment length.
- The defaults are half-window overlap, an FFT the window's length, and no detrending.
- `Values` are `[segment][frequency]`, the transpose of MATLAB's `S` and `P`.
- Two-sided frequencies run from 0 up then negative, the same as `fft.Freq`, where MATLAB's run up to fs.
- `AverageMedian()` and `SegmentDetrend` add a median average and per-segment detrending, which `pwelch` doesn't have.

### Reductions and calculus

| MATLAB                | tukey                                                   |
| --------------------- | ------------------------------------------------------- |
| `s = sum(x)`          | `s := dsp.Sum(x)`                                       |
| `m = mean(x)`         | `m, err := dsp.Mean(x)`                                 |
| `r = rms(x)`          | `r, err := dsp.RMS(x)`                                  |
| `v = prctile(x, p)`   | `v, err := dsp.Percentile(x, p, dsp.PercentileHazen())` |
| `d = diff(x)`         | `d := dsp.Diff(x)`                                      |
| `d = diff(x, 2)`      | `d := dsp.Diff(dsp.Diff(x))`                            |
| `g = gradient(x, dx)` | `g, err := dsp.Grad(x, dx)`                             |
| `s = cumsum(x)`       | `s := dsp.CumSum(x)`                                    |
| `q = trapz(dx, y)`    | `q, err := dsp.Trapz(y, dx)`                            |
| `q = cumtrapz(dx, y)` | `q, err := dsp.CumTrapz(y, dx)`                         |

#### Differences to note:

- Sums are compensated, so they're closer to the exact answer than MATLAB's plain sums,
  e.g. 1e8, then 20,000 values of 0.001, then -1e8, sums to 20, where `sum` gives 20.00001.
- `Mean` and `RMS` of an empty signal return an error, where MATLAB gives NaN.
- `Percentile` is linear by default, and `PercentileHazen()` is `prctile`.
  Any NaN gives NaN, where `prctile` ignores it, so use `FillMissing` first.
  It takes one `p` per call, where `prctile` takes a vector.
- `CumTrapz` starts at 0 and has the same length as `y`, the same as `cumtrapz`.
