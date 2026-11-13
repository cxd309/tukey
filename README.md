# tukey

Offline Digital Signal Processing (DSP) for GO, for engineers, checked against SciPy and MATLAB.

```go
import "github.com/cxd309/tukey/dsp"

f, err := dsp.Butter(2, dsp.Lowpass(0.3))
b, a, err := f.BA()
y, err := dsp.FiltFilt(b, a, x)

// higher orders and narrow bands: second-order sections
f, err = dsp.Butter(8, dsp.Bandpass(0.30, 0.32))
sos, err := f.SOS()
y, err = dsp.SOSFiltFilt(sos, x)
```

```go
// the lowest order that meets a specification, then the design
order, band, err := dsp.EllipOrd(dsp.Lowpass(0.2), dsp.Lowpass(0.25), 0.5, 60)
f, err := dsp.Ellip(order, 0.5, 60, band)

y, err := dsp.MedFilt(x, 5)         // spikes removed, steps kept
y, err = dsp.SavGolFilter(x, 21, 3) // noise smoothed, peaks kept
```

```go
import "github.com/cxd309/tukey/fft"

spectrum, err := fft.ForwardReal(x)  // any length, O(n log n)
freqs := fft.RealFreq(len(x), 1/fs)  // each bin's frequency in Hz
```

```go
import "github.com/cxd309/tukey/window"

w, err := window.Hann(len(x), window.Periodic()) // periodic for spectral analysis
w, err = window.Kaiser(64, 8.6)                  // symmetric, for filter design
```

```go
w, err := window.Hann(1024, window.Periodic())  // the window sets the segment length
psd, err := dsp.Welch(x, fs, w)                 // power spectral density, half-overlapped segments
tf, err := dsp.Spectrogram(x, fs, w)            // power over time: tf.Values[segment][frequency]
```

Named after [John Tukey](https://en.wikipedia.org/wiki/John_Tukey), who co-invented the FFT and considered by some as the father of Data Science.

## Why

This library was designed to support the development of appliations reliant on DSP in GO.

A survey of existing tools found some coverage but was specifically lacking zero-phase filtering which was essential for an application that we needed to develop.

Therefore the decision was made to build a library to support the application being developed and to make it public source so that others might not encounter the same speedbump in development.

The specific application was porting from an existing MATLAB web app which heavily relied on both IIR filter design and `filtfilt`, the forward-backward filtering that removes phase-distortion.

Phase distortion matters whenever the location of a feature in a signal matters. A causal filter delays each frequency by a different amount, which moves features: in the example below, a pulse at sample 50 comes at sample 54 while `FiltFilt` smooths it and leaves it at 50.

Python development with SciPy was considered but the I preferred GO tooling for web app development and therefore decided to develop this package instead. I have used SciPy to compare outputs against, as a comparison test.

## Compatibility contract

Every function matches its SciPy equivalent to within a measured tolerance. Each figure is the worst on either amd64 or arm64, which round slightly differently because arm64 fuses multiply-adds.

| Function          | SciPy equivalent                        | Worst measured | Tested to |
| ----------------- | --------------------------------------- | -------------- | --------- |
| `Butter`          | `butter`, as `(b, a)`                   | 7.6e-15        | 1e-13     |
| `Butter` + SOS    | `butter(output='sos')`                  | 7.3e-15        | 1e-13     |
| `Cheby1`          | `cheby1`, as `(b, a)`                   | 6.8e-15        | 1e-13     |
| `Cheby1` + SOS    | `cheby1(output='sos')`                  | 6.8e-15        | 1e-13     |
| `Cheby2`          | `cheby2`, as `(b, a)`                   | 5.7e-14        | 1e-12     |
| `Cheby2` + SOS    | `cheby2(output='sos')`                  | 5.7e-14        | 1e-12     |
| `Ellip`           | `ellip`, as `(b, a)`, exactly           | 3.9e-15        | 1e-13     |
| `Ellip` + SOS     | `ellip(output='sos')`, exactly          | 3.0e-15        | 1e-13     |
| `Bessel`          | `bessel`, as `(b, a)`, exactly          | 6.7e-15        | 1e-13     |
| `Bessel` + SOS    | `bessel(output='sos')`, exactly         | 6.7e-15        | 1e-13     |
| `ButterOrd`       | `buttord`, exactly                      | 4.2e-16        | 1e-14     |
| `Cheby1Ord`       | `cheb1ord`, exactly                     | 2.9e-16        | 1e-14     |
| `Cheby2Ord`       | `cheb2ord`, exactly                     | 3.9e-16        | 1e-14     |
| `EllipOrd`        | `ellipord`, exactly                     | 7.5e-16        | 1e-14     |
| `IIRNotch`        | `iirnotch`, exactly                     | 4.7e-16        | 1e-14     |
| `IIRPeak`         | `iirpeak`, exactly                      | 4.7e-16        | 1e-14     |
| `FIRWin`          | `firwin`, exactly                       | 3.5e-16        | 1e-14     |
| `ZPK.SOS`         | `zpk2sos`                               | 1.1e-16        | 1e-14     |
| `Filter`          | `lfilter`                               | 6.4e-16        | 1e-14     |
| `FilterState`     | `lfilter` with `zi`                     | 2.3e-16        | 1e-14     |
| `FilterZi`        | `lfilter_zi`                            | 3.3e-16        | 1e-14     |
| `FiltFilt`        | `filtfilt`                              | 8.8e-16        | 1e-14     |
| `SOSFilter`       | `sosfilt`, with and without `zi`        | 4.1e-14        | 1e-12     |
| `SOSFilterZi`     | `sosfilt_zi`                            | 9.8e-17        | 1e-14     |
| `SOSFiltFilt`     | `sosfiltfilt`                           | 1.7e-13        | 1e-11     |
| `Freqz`           | `freqz`                                 | 1.4e-15        | 1e-13     |
| `SOSFreqz`        | `freqz_sos`                             | 3.7e-13        | 1e-11     |
| `GroupDelay`      | `group_delay`, exactly                  | 7.2e-15        | 1e-13     |
| `SOSGroupDelay`   | `group_delay` of the sections, exactly  | 1.0e-14        | 1e-13     |
| `Detrend`         | `detrend`                               | 5.6e-16        | 1e-14     |
| `Interp1`         | `interp1d(bounds_error=False)`          | 1.3e-16        | 1e-14     |
| `RMS`             | `sqrt(mean(x**2))`, from NumPy          | 1.4e-16        | 1e-14     |
| `MedFilt`         | `medfilt`, exactly, and truncated       | exact          | 1e-14     |
| `SavGolFilter`    | `savgol_filter(mode='interp')`, exactly | 3.6e-16        | 1e-14     |
| `Convolve`        | `convolve`, exactly                     | exact          | 1e-14     |
| `Convolve` + FFT  | `convolve`, exactly                     | 5.5e-16        | 1e-14     |
| `Correlate`       | `correlate`, exactly                    | exact          | 1e-14     |
| `Correlate` + FFT | `correlate`, exactly                    | 4.7e-16        | 1e-14     |
| `FindPeaks`       | `find_peaks`                            | exact          | exact     |
| `PeakProminences` | `peak_prominences`                      | exact          | 1e-14     |
| `PeakWidths`      | `peak_widths`                           | exact          | 1e-14     |
| `Hilbert`         | `hilbert`                               | 5.7e-16        | 1e-14     |
| `Envelope`        | MATLAB's `envelope`, from `hilbert`     | 1.0e-15        | 1e-14     |
| `STFT`            | `stft`, unpadded, × `sum(window)`       | 4.2e-16        | 1e-14     |
| `Spectrogram`     | `spectrogram(mode='psd')`               | 7.3e-16        | 1e-14     |
| `Welch`           | `welch`                                 | 1.1e-15        | 1e-13     |
| `Periodogram`     | `periodogram`                           | 1.1e-15        | 1e-13     |
| `Sum`             | `sum`, exactly (`math.fsum`)            | exact          | 1e-14     |
| `Mean`            | `mean`, exactly                         | exact          | 1e-14     |
| `Diff`            | `numpy.diff`                            | exact          | 1e-14     |
| `Grad`            | `numpy.gradient`                        | exact          | 1e-14     |
| `CumSum`          | `cumsum`, exactly                       | exact          | 1e-14     |
| `Trapz`           | `trapezoid`, exactly                    | exact          | 1e-14     |
| `CumTrapz`        | `cumulative_trapezoid`, exactly         | exact          | 1e-14     |
| `Percentile`      | `numpy.percentile`                      | exact          | 1e-14     |

For scale, 1e-16 is ~ 1 ULP of a `float64`. The larger SOS figures all come from one deliberately difficult test filter (order 8, cutoff 0.01), whose poles sit so close to z = 1 that single rounding differences between Go and SciPy are amplified. Every other filter agrees within a few ULP.

Where SciPy's own answer is measurably off, its reference is instead the exact answer, computed with mpmath at 50 digits or in rational arithmetic and rounded once, marked "exactly" above. SciPy itself is off by 1.5e-13 for `ellip`, 2.7e-14 for `bessel`, 6.1e-13 for `iirpeak` at w0 = 0.002, Q = 30, as it forms the gain as 1 − 1/(1 + β), and 3.0e-15 for `firwin`. Its four order selections move a bandstop's edges by a numerical minimiser, good to about 1e-5. Its `group_delay` is 0 at a zero on the unit circle and 3% off for `butter(8, 0.01)`. Its `savgol_filter` fits raw powers of the positions by `lstsq`, 1e-12 off at n = 21 and wrong by 21% at n = 101, degree 10. And its `medfilt` only zero-pads, with NaN handled however it happens to sort.

Two very low cutoff filters, order 2 at 0.001 low-pass and 0.0008 high-pass, are checked separately. They're the original MATLAB application's own designs, and their poles sit even closer to z = 1, so rounding differences between implementations are magnified by roughly the square of 1/(1 − pole):

| Function      | Worst against SciPy | Tested to | Worst against MATLAB | Tested to |
| ------------- | ------------------- | --------- | -------------------- | --------- |
| `Filter`      | 9.5e-13             | 1e-11     | 6.2e-13              | 1e-11     |
| `FilterState` | 9.9e-14             | 1e-12     | 1.1e-13              | 1e-11     |
| `FiltFilt`    | 4.4e-12             | 1e-10     | 1.5e-11              | 1e-9      |
| `Freqz`       | 1.5e-12             | 1e-10     | 1.1e-11              | 1e-9      |

SciPy and MATLAB differ from each other on these filters by as much. At DC the 0.001 low-pass's gain is exactly 1, and `Freqz` gives 1 + 1.6e-12 where MATLAB gives 1 + 1.3e-11, so much of the difference is MATLAB's.

Functions MATLAB also has are checked against vectors generated by MATLAB R2025a (Signal Processing Toolbox 25.1), on the same inputs:

| Function                      | MATLAB equivalent                    | Worst measured | Tested to |
| ----------------------------- | ------------------------------------ | -------------- | --------- |
| `Butter`                      | `[b, a] = butter`                    | 1.2e-13        | 1e-11     |
| `Cheby1`                      | `[b, a] = cheby1`                    | 6.8e-13        | 1e-11     |
| `Cheby2`                      | `[b, a] = cheby2`                    | 1.8e-14        | 1e-12     |
| `Ellip`                       | `[b, a] = ellip`                     | 9.9e-9         | 1e-7      |
| `Bessel` + `BesselNormPhase`  | `besself`, then `bilinear`           | 1.6e-13        | 1e-11     |
| `ButterOrd` + `MatchStopband` | `buttord`                            | 3.9e-16        | 1e-14     |
| `Cheby1Ord`                   | `cheb1ord`                           | 1.8e-16        | 1e-14     |
| `Cheby2Ord` + `MatchStopband` | `cheb2ord`                           | 2.2e-16        | 1e-14     |
| `EllipOrd`                    | `ellipord`                           | 1.8e-16        | 1e-14     |
| `IIRNotch`                    | `iirnotch(w0, w0/q)`                 | 3.8e-16        | 1e-14     |
| `IIRPeak`                     | `iirpeak(w0, w0/q)`                  | 6.1e-13        | 1e-11     |
| `FIRWin`                      | `fir1(numel(w) - 1, wn, w)`          | 3.2e-15        | 1e-13     |
| `Filter`                      | `filter`                             | 6.4e-16        | 1e-14     |
| `FilterState`                 | `[y, zf] = filter(b, a, x, zi)`      | 3.0e-16        | 1e-14     |
| `FiltFilt` + MATLABPadLen     | `filtfilt(b, a, x)`                  | 8.8e-16        | 1e-14     |
| `SOSFilter`                   | `sosfilt`                            | 4.5e-14        | 1e-12     |
| `Freqz`                       | `freqz(b, a, n)`                     | 1.9e-15        | 1e-13     |
| `SOSFreqz`                    | `freqz(sos, n)`                      | 1.1e-12        | 1e-10     |
| `GroupDelay`                  | `[gd, w] = grpdelay(b, a, n)`        | 3.2e-11        | 1e-9      |
| `SOSGroupDelay`               | `grpdelay(sos, n)`                   | 5.7e-9         | 1e-7      |
| `Detrend`                     | `detrend`                            | 1.2e-15        | 1e-13     |
| `Interp1`                     | `interp1`                            | 2.2e-16        | 1e-14     |
| `FillMissing`                 | `fillmissing`                        | exact          | 1e-14     |
| `RMS`                         | `rms`                                | 1.2e-16        | 1e-14     |
| `MedFilt`                     | `medfilt1(x, n, 'truncate')`         | 2.8e-17        | 1e-14     |
| `SavGolFilter`                | `sgolayfilt(x, order, n)`            | 8.7e-16        | 1e-14     |
| `Convolve`                    | `conv`                               | 1.6e-15        | 1e-13     |
| `Convolve` + FFT              | `conv`                               | 1.8e-15        | 1e-13     |
| `Correlate`                   | `[r, lags] = xcorr(a, v)`            | 3.8e-16        | 1e-14     |
| `Correlate` + FFT             | `[r, lags] = xcorr(a, v)`            | 6.3e-16        | 1e-14     |
| `FindPeaks` + `MATLABPeaks`   | `findpeaks`                          | exact          | exact     |
| `PeakProminences`             | `[~, ~, ~, p] = findpeaks`           | exact          | 1e-14     |
| `PeakWidths`                  | `[~, ~, w] = findpeaks`              | 2.7e-15        | 1e-13     |
| `Hilbert`                     | `hilbert`                            | 5.8e-16        | 1e-14     |
| `Envelope`                    | `envelope`, `envelope(x, n, 'rms')`  | 4.5e-15        | 1e-13     |
| `FindSignal`                  | `findsignal(data, signal)`           | 5.9e-16        | 1e-14     |
| `FindSignal` + `SearchFast`   | `findsignal(data, signal)`           | 1.5e-15        | 1e-13     |
| `STFT`                        | `[S, F, T] = spectrogram`            | 4.2e-16        | 1e-14     |
| `STFT` + `AtFrequencies`      | `spectrogram(x, w, noverlap, f, fs)` | 2.4e-13        | 1e-11     |
| `Spectrogram`                 | `[~, F, T, P] = spectrogram`         | 7.3e-16        | 1e-14     |
| `Welch`                       | `pwelch`                             | 1.3e-15        | 1e-13     |
| `Periodogram`                 | `periodogram`                        | 5.0e-16        | 1e-14     |
| `Sum`                         | `sum`                                | 9.0e-14        | 1e-12     |
| `Mean`                        | `mean`                               | 9.0e-14        | 1e-12     |
| `Diff`                        | `diff`                               | exact          | 1e-14     |
| `Grad`                        | `gradient(x, dx)`                    | 2.1e-14        | 1e-12     |
| `CumSum`                      | `cumsum`                             | 3.6e-13        | 1e-11     |
| `Trapz`                       | `trapz(dx, y)`                       | 9.0e-14        | 1e-12     |
| `CumTrapz`                    | `cumtrapz(dx, y)`                    | 3.6e-13        | 1e-11     |
| `Percentile` + Hazen          | `prctile`                            | 4.0e-17        | 1e-14     |

The largest come from the same difficult filter, where MATLAB and SciPy differ from each other by as much. SciPy has no `fillmissing`, so `FillMissing` is checked against MATLAB's documented definition as well as MATLAB itself. MATLAB's RMS envelope is 4.5e-15 from an exactly summed answer, where `Envelope`'s compensated sum is within 1e-15 of it, so that gap is MATLAB's. SciPy has no `findsignal` either, so `FindSignal` is also checked against MATLAB's definition summed exactly, which it matches exactly. MATLAB's distance is measured relative to the window's and signal's sums of squares, the scale it rounds at: its perfect matches are up to 2^-29, not 0. `STFT` with `AtFrequencies` is within 1 ULP of exact sums, which SciPy has no equivalent for. MATLAB's `spectrogram` at a frequency vector is 2.4e-13 from them, even at the FFT's own bins, where its FFT path agrees within 3e-16. The sums are compensated, and match exact sums on every SciPy-side case. MATLAB sums plainly, so its sums are measured relative to the sum of the magnitudes summed, the scale it rounds at: its running sum of 20,000 tenths drifts by 3.6e-13 of that, and 1e8 among 20,000 thousandths, then cancelled, sums to 20.00001 where `Sum` gives 20. See [Porting from MATLAB](#porting-from-matlab) for the functions not in this table.

Of the filter designs, MATLAB's own rounding sets the larger figures. Its `ellip` forms 1 − k1² for the specification's small modulus k1 and loses its digits, 9.9e-9 from the exact design, where SciPy and tukey keep them. Its `iirpeak` forms the gain as SciPy's does, 6.1e-13 off, and its `cheby1` differs from SciPy's by as much as from tukey's. MATLAB has no digital Bessel design, so `Bessel` is checked against `besself` at the prewarped edges through `bilinear`, in the phase normalisation and up to order 25, `besselap`'s highest, which is itself 1.7e-13 from the exact design. `iirnotch` and `iirpeak` need the DSP System Toolbox. The order selections' band is compared only where MATLAB follows one of tukey's conventions (see below), and its bandstop `buttord` edges come from the same minimiser as SciPy's, 5.7e-5 off.

MATLAB's `grpdelay` loses digits near poles and zeros crowded against the unit circle: 3.2e-11 at worst on `(b, a)`, 2.2e-5 on the two very low cutoff filters, 5.7e-9 on `butter(8, 0.01)`'s sections, and 7.3% off on its `(b, a)`, which isn't compared, where `GroupDelay` is exact. MATLAB's `xcorr` is by FFT, 2.4e-4 off the exact answer on small outputs beside large ones, though within 6.3e-16 of the largest, the scale these comparisons use. Its `medfilt1` is 1 ULP off at some even-count midpoints near the truncated ends, and rounds a midpoint of two subnormals twice. Its `sgolayfilt` loses digits at a 101-sample window of degree 10, 2.5e-14 off, compared separately.

The `fft` package matches `scipy.fft` the same way, at 23 lengths covering powers of two, smooth lengths and primes:

| Function            | SciPy equivalent        | Worst measured | Tested to |
| ------------------- | ----------------------- | -------------- | --------- |
| `fft.Forward`       | `fft`                   | 9.3e-16        | 1e-14     |
| `fft.Inverse`       | `ifft`                  | 8.5e-16        | 1e-14     |
| `fft.ForwardReal`   | `rfft`                  | 5.8e-16        | 1e-14     |
| `fft.InverseReal`   | `irfft`                 | 5.0e-16        | 1e-14     |
| `fft.Norm`, `fft.N` | `norm=`, `n=`           | 6.9e-16        | 1e-14     |
| `fft.Freq`          | `fftfreq`, `rfftfreq`   | exact          | exact     |
| `fft.Shift`         | `fftshift`, `ifftshift` | exact          | exact     |

And against MATLAB, whose `fft` is FFTW, an independent implementation:

| Function          | MATLAB equivalent         | Worst measured | Tested to |
| ----------------- | ------------------------- | -------------- | --------- |
| `fft.Forward`     | `fft`                     | 9.5e-16        | 1e-14     |
| `fft.Inverse`     | `ifft`                    | 8.3e-16        | 1e-14     |
| `fft.ForwardReal` | `fft`, first n/2+1 bins   | 6.7e-16        | 1e-14     |
| `fft.InverseReal` | `ifft(x, n, 'symmetric')` | 5.3e-16        | 1e-14     |
| `fft.N`           | `fft(x, n)`, `ifft(x, n)` | 7.2e-16        | 1e-14     |
| `fft.Shift`       | `fftshift`, `ifftshift`   | exact          | exact     |

MATLAB has no normalisation option or `fftfreq`, so `fft.Norm` and `fft.Freq` are checked against SciPy only.

The `window` package has every `scipy.signal.windows` window except `dpss`, plus MATLAB's `hanning`. Each is checked against SciPy, and against MATLAB where MATLAB has it, on identical inputs at n = 0, 1, 2, 3, 8, 9, 64 and 255, symmetric and periodic:

| Function                     | SciPy                   | MATLAB           | SciPy worst | MATLAB worst | Tested to |
| ---------------------------- | ----------------------- | ---------------- | ----------- | ------------ | --------- |
| `window.Hann`                | `hann`                  | `hann`           | 4.4e-16     | 1.7e-16      | 1e-14     |
| `window.Hamming`             | `hamming`               | `hamming`        | 4.4e-16     | 2.2e-16      | 1e-14     |
| `window.Blackman`            | `blackman`              | `blackman`       | 5.6e-16     | 3.3e-16      | 1e-14     |
| `window.BlackmanHarris`      | `blackmanharris`        | `blackmanharris` | 6.7e-16     | 7.8e-16      | 1e-14     |
| `window.Nuttall`             | `nuttall`               | `nuttallwin`     | 6.7e-16     | 8.9e-16      | 1e-14     |
| `window.FlatTop`             | `flattop`               | `flattopwin`     | 7.5e-16     | 3.3e-16      | 1e-14     |
| `window.GeneralHamming`      | `general_hamming`       |                  | 4.4e-16     |              | 1e-14     |
| `window.GeneralCosine`       | `general_cosine`        |                  | 1.7e-15     |              | 1e-13     |
| `window.Boxcar`              | `boxcar`                | `rectwin`        | exact       | exact        | exact     |
| `window.Triang`              | `triang`                | `triang`         | exact       | exact        | exact     |
| `window.Bartlett`            | `bartlett`              | `bartlett`       | 1.1e-16     | exact        | 1e-14     |
| `window.Parzen`              | `parzen`                | `parzenwin`      | 2.2e-16     | 1.1e-16      | 1e-14     |
| `window.Bohman`              | `bohman`                | `bohmanwin`      | 5.6e-16     | 2.2e-16      | 1e-14     |
| `window.Barthann`            | `barthann`              | `barthannwin`    | 3.3e-16     | 3.3e-16      | 1e-14     |
| `window.Tukey`               | `tukey`                 | `tukeywin`       | 6.7e-16     | 6.7e-16      | 1e-14     |
| `window.Hanning`             | `hann(n+2)[1:-1]`       | `hanning`        | 3.3e-16     | 1.7e-16      | 1e-14     |
| `window.Cosine`              | `cosine`                |                  | 4.6e-16     |              | 1e-14     |
| `window.Lanczos`             | `lanczos`               |                  | 2.2e-16     |              | 1e-14     |
| `window.Gaussian`            | `gaussian`              | `gausswin`       | 1.1e-16     | 1.7e-16      | 1e-14     |
| `window.GeneralGaussian`     | `general_gaussian`      |                  | 1.1e-16     |              | 1e-14     |
| `window.Exponential`         | `exponential`           |                  | 1.1e-16     |              | 1e-14     |
| `window.Kaiser`              | `kaiser`                | `kaiser`         | 1.7e-15     | 2.0e-15      | 1e-13     |
| `window.KaiserBesselDerived` | `kaiser_bessel_derived` |                  | exact       |              | 1e-13     |
| `window.Taylor`              | `taylor`                | `taylorwin`      | 5.6e-16     | 3.4e-16      | 1e-14     |
| `window.Chebwin`             | `chebwin`               | `chebwin`        | 1.5e-12     | 3.7e-13      | 1e-11     |

Kaiser's error grows with beta, as the Bessel function's sensitivity to its argument does. Chebwin's comes from evaluating a Chebyshev polynomial of degree n−1 near ±1, where an ULP of difference in its argument is magnified; SciPy and MATLAB differ from each other by as much. MATLAB's windows without a `'periodic'` option are checked symmetric only.

These tables are generated by:

```sh
go run ./cmd/tolerances
```

## Porting from MATLAB

Function calls map one-to-one. Where MATLAB's behaviour differs from SciPy's, SciPy is chosen and the default with the option to implement MATLAB behaviour.

| MATLAB                                              | tukey                                                                                                    |
| --------------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| `[b, a] = butter(n, Wn)`                            | `f, err := dsp.Butter(n, dsp.Lowpass(Wn))`, `b, a, err := f.BA()`                                        |
| `butter(n, Wn, 'high')`                             | `dsp.Butter(n, dsp.Highpass(Wn))`                                                                        |
| `butter(n, [W1 W2])`                                | `dsp.Butter(n, dsp.Bandpass(W1, W2))`                                                                    |
| `butter(n, [W1 W2], 'stop')`                        | `dsp.Butter(n, dsp.Bandstop(W1, W2))`                                                                    |
| `[z, p, k] = butter(n, Wn)`                         | `f, err := dsp.Butter(n, dsp.Lowpass(Wn))`: `f.Zeros`, `f.Poles`, `f.Gain`                               |
| `cheby1(n, Rp, Wn)`, `cheby2(n, Rs, Wn)`            | `dsp.Cheby1(n, Rp, dsp.Lowpass(Wn))`, `dsp.Cheby2(n, Rs, dsp.Lowpass(Wn))`                               |
| `ellip(n, Rp, Rs, Wn)`                              | `dsp.Ellip(n, Rp, Rs, dsp.Lowpass(Wn))`                                                                  |
| `[z, p, k] = besself(n, Wo)`, then `bilinear`       | `dsp.Bessel(n, dsp.Lowpass(Wn), dsp.BesselNormalised(dsp.BesselNormPhase))`, digital directly            |
| `[n, Wn] = buttord(Wp, Ws, Rp, Rs)`                 | `n, wn, err := dsp.ButterOrd(dsp.Lowpass(Wp), dsp.Lowpass(Ws), Rp, Rs, dsp.Matching(dsp.MatchStopband))` |
| `cheb1ord(…)`, `ellipord(…)`                        | `dsp.Cheby1Ord(pass, stop, Rp, Rs)`, `dsp.EllipOrd(pass, stop, Rp, Rs)`                                  |
| `cheb2ord(…)`                                       | `dsp.Cheby2Ord(pass, stop, Rp, Rs, dsp.Matching(dsp.MatchStopband))`                                     |
| `[b, a] = iirnotch(w0, bw)`, `iirpeak(w0, bw)`      | `f, err := dsp.IIRNotch(w0, w0/bw)`, `dsp.IIRPeak(w0, w0/bw)`: q is w0/bw                                |
| `b = fir1(n, Wn)`                                   | `h, err := dsp.FIRWin(window.Hamming(n+1), dsp.Lowpass(Wn))`: the window's length is the number of taps  |
| `b = fir1(n, Wn, w)`                                | `h, err := dsp.FIRWin(w, dsp.Lowpass(Wn))`                                                               |
| `y = filter(b, a, x)`                               | `y, err := dsp.Filter(b, a, x)`                                                                          |
| `y = filtfilt(b, a, x)`                             | `y, err := dsp.FiltFilt(b, a, x, dsp.MATLABPadLen())`                                                    |
| `[h, w] = freqz(b, a, n)`                           | `w, h, err := dsp.Freqz(b, a, n)` (note the order)                                                       |
| `[gd, w] = grpdelay(b, a, n)`, `grpdelay(sos, n)`   | `w, gd, err := dsp.GroupDelay(b, a, n)`, `dsp.SOSGroupDelay(sos, n)` (note the order)                    |
| `y = sosfilt(sos, x)`                               | `y, err := dsp.SOSFilter(sos, x)`                                                                        |
| `y = detrend(x)`                                    | `y, err := dsp.Detrend(x, dsp.DetrendLinear)`                                                            |
| `y = detrend(x, 'constant')`                        | `y, err := dsp.Detrend(x, dsp.DetrendConstant)`                                                          |
| `yi = interp1(x, y, xi)`                            | `yi, err := dsp.Interp1(x, y, xi)`                                                                       |
| `fillmissing(x, 'linear')`                          | `y := dsp.FillMissing(x, dsp.FillLinear())`                                                              |
| `fillmissing(x, 'constant', v)`                     | `y := dsp.FillMissing(x, dsp.FillConstant(v))`                                                           |
| `y = medfilt1(x, n, 'truncate')`                    | `y, err := dsp.MedFilt(x, n)`                                                                            |
| `y = medfilt1(x, n)`                                | `y, err := dsp.MedFilt(x, n, dsp.MedFiltEdges(dsp.MedFiltZeroPad))`                                      |
| `y = sgolayfilt(x, order, n)`                       | `y, err := dsp.SavGolFilter(x, n, order)` (note the order)                                               |
| `c = conv(a, v)`, `conv(a, v, 'same')`              | `c, err := dsp.Convolve(a, v)`, `dsp.Convolve(a, v, dsp.ConvolveSame())`                                 |
| `[r, lags] = xcorr(a, v)`                           | `c, err := dsp.Correlate(a, v)`: `c.Values`, `c.Lags`                                                    |
| `r = rms(x)`                                        | `r, err := dsp.RMS(x)`                                                                                   |
| `[pks, locs] = findpeaks(x)`                        | `peaks, err := dsp.FindPeaks(x, dsp.MATLABPeaks())`                                                      |
| `findpeaks(x, 'MinPeakDistance', d)`                | `dsp.FindPeaks(x, dsp.MATLABPeaks(), dsp.MinPeakDistance(d+1))`                                          |
| `'MinPeakHeight'`, `'Threshold'`                    | `dsp.MinPeakHeight(h)`, `dsp.PeakThreshold(t)`                                                           |
| `'MinPeakProminence'`, `'MinPeakWidth'`             | `dsp.MinPeakProminence(p)`, `dsp.MinPeakWidth(w)`                                                        |
| `[~, ~, w, p] = findpeaks(x)`                       | `dsp.PeakWidths(x, indices, 0.5)`, `dsp.PeakProminences(x, indices)`                                     |
| `y = hilbert(x)`                                    | `y, err := dsp.Hilbert(x)`                                                                               |
| `[up, lo] = envelope(x)`                            | `up, lo, err := dsp.Envelope(x, dsp.EnvelopeAnalytic())`                                                 |
| `[up, lo] = envelope(x, n, 'rms')`                  | `up, lo, err := dsp.Envelope(x, dsp.EnvelopeRMS(n))`                                                     |
| `[istart, istop, dist] = findsignal(data, signal)`  | `match, err := dsp.FindSignal(data, signal)`: `match.Start`, `match.Stop`, `match.Distance`              |
| `[S, F, T] = spectrogram(x, w, noverlap, nfft, fs)` | `tf, err := dsp.STFT(x, fs, w, dsp.Overlap(noverlap), dsp.NFFT(nfft))`: `tf.Values[segment][frequency]`  |
| `spectrogram(x, w, noverlap, f, fs)`                | `dsp.STFT(x, fs, w, dsp.Overlap(noverlap), dsp.AtFrequencies(f))`                                        |
| `[~, F, T, P] = spectrogram(…, 'power')`            | `dsp.Spectrogram(x, fs, w, dsp.Overlap(noverlap), dsp.NFFT(nfft), dsp.Scaling(dsp.ScalingSpectrum))`     |
| `[pxx, f] = pwelch(x, w, noverlap, nfft, fs)`       | `psd, err := dsp.Welch(x, fs, w, dsp.Overlap(noverlap), dsp.NFFT(nfft))`: `psd.Frequencies`, `psd.Power` |
| `[pxx, f] = periodogram(x, w, nfft, fs)`            | `psd, err := dsp.Periodogram(x, fs, w, dsp.NFFT(nfft))`                                                  |
| `s = sum(x)`, `m = mean(x)`                         | `s := dsp.Sum(x)`, `m, err := dsp.Mean(x)`                                                               |
| `d = diff(x)`, `diff(x, 2)`                         | `d := dsp.Diff(x)`, `dsp.Diff(dsp.Diff(x))`                                                              |
| `g = gradient(x, dx)`                               | `g, err := dsp.Grad(x, dx)`                                                                              |
| `s = cumsum(x)`                                     | `s := dsp.CumSum(x)`                                                                                     |
| `q = trapz(dx, y)`, `cumtrapz(dx, y)`               | `q, err := dsp.Trapz(y, dx)`, `dsp.CumTrapz(y, dx)` (note the order)                                     |
| `v = prctile(x, p)`                                 | `v, err := dsp.Percentile(x, p, dsp.PercentileBy(dsp.PercentileHazen))`, one `p` per call                |

Frequencies are normalised to Nyquist in both: for a frequency `f`Hz at sample rate `fs`, use `f / (fs/2)`. `Freqz` returns frequencies in rad/sample, as MATLAB's `freqz` does.

`MATLABPadLen` uses MATLAB's `filtfilt` padding length, `3*(nfilt-1)`, instead of SciPy's `3*nfilt`. Only samples near the ends of the output are affected.

Where tukey and MATLAB differ, verified against MATLAB R2025a:

- **Second-order sections:** `f.SOS()` pairs and orders sections as SciPy's `zpk2sos` does, with the overall gain in the first section. MATLAB's `zp2sos` gives a first-order section one zero rather than two, and breaks ties in pole radius differently: the same filter, in different sections. `[sos, g] = zp2sos(z, p, k)` also returns the gain `g` separately
- **`SOSFiltFilt` with `MATLABPadLen`** does not match MATLAB's `filtfilt(sos, g, x)`, which filters forward then backward one section at a time, each padded by 3 × the total filter order. Expect differences near the ends of the output; a MATLAB-compatible mode is on the roadmap
- **Filter design** takes a `Band`, `dsp.Lowpass(Wn)` and so on, rather than a type string, and returns zeros, poles and gain, for `BA()` or `SOS()`. `Bessel`'s band edge is the −3 dB point by default, as `Butter`'s. MATLAB has no digital Bessel design: `besself` is analog, in the phase normalisation, `BesselNormalised(BesselNormPhase)`
- **Order selection** takes the passband and stopband as `Band`s. MATLAB's `buttord` and `cheb2ord` meet the stopband edge exactly, so pass `Matching(MatchStopband)`. For a bandstop, `cheb1ord` and `ellipord` return the passband edges unmoved, which for `cheb1ord` can miss the stopband's specification, where tukey moves them as SciPy does. For a bandpass, `cheb2ord` returns both stopband edges, where tukey meets the tighter one. Every order matches
- **`IIRNotch` and `IIRPeak`** take the quality factor q, as SciPy, where `iirnotch` and `iirpeak` take the bandwidth bw = w0/q. `IIRPeak` forms its gain exactly, where `iirpeak` is 6.1e-13 off for a narrow bandwidth
- **`FIRWin`** takes the window itself, whose length is the number of taps, where `fir1` takes the order, one less, and defaults to `hamming(n+1)`. A highpass or bandstop needs an odd number of taps, and an even one returns `ErrInvalidLength`
- **`GroupDelay`** returns the frequencies first, as `Freqz`, and is exact where `grpdelay` loses digits near poles crowded against the unit circle. At a zero on the unit circle both give the limit from either side
- **`Convolve`'s `ConvolveSame`** starts one output earlier than `conv(a, v, 'same')` for an even-length `v`, as SciPy's does, and `ConvolveValid` with `v` longer than `a` gives `len(v) − len(a) + 1` outputs, as NumPy, where `conv` gives none. `ConvolveAuto` sums by FFT for all but short inputs: set `ConvolveBy(ConvolveDirect)` for inputs with high dynamic range
- **`Correlate`** gives the lags where `a` and `v` overlap, `−(len(v) − 1)` to `len(a) − 1`, where `xcorr` pads the shorter with zeros and gives every lag up to the longer's length. The lags run the same way. `xcorr`'s normalisations aren't supported yet
- **`MedFilt`** truncates the window at the ends by default, as `medfilt1(x, n, 'truncate')`, where MATLAB's and SciPy's default zero padding pulls the ends towards 0. n must be odd
- **`SavGolFilter`** takes `(x, n, polyOrder)`, where `sgolayfilt` takes `(x, order, framelen)`, and gives derivatives with `SavGolDeriv`, which `sgolayfilt` doesn't. n must be odd and no longer than `x`
- **`Freqz` with n = 0** returns empty slices, as SciPy does, where MATLAB's `freqz` errors
- **`RMS` of an empty signal** returns `ErrInvalidLength`, where MATLAB's `rms` gives NaN
- **`FindPeaks`** returns 0-based indices, and follows SciPy's rules unless given `MATLABPeaks()`. MATLAB's `MinPeakDistance` d keeps peaks more than d apart, so it's `MinPeakDistance(d+1)`. MATLAB's `findpeaks` errors for fewer than 3 samples, and for a distance of n − 1 or more, where tukey returns the peaks as usual. `'SortStr'`, `'NPeaks'` and `'WidthReference'` aren't supported yet
- **`Envelope`** supports MATLAB's default analytic envelope and `'rms'`, not `envelope(x, fl, 'analytic')` with an FIR Hilbert filter, or `'peak'`
- **`FindSignal`** returns 0-based indices, and by default the exact distance of the closest window, where MATLAB's `findsignal` rounds at the scale of the sums of squares, so its perfect matches can be up to 2⁻²⁹ rather than 0. `SignalSearch(SearchFast)` rounds as MATLAB does. Only the default options are supported so far: `'Normalization'`, `'Metric'`, `'TimeAlignment'`, `'MaxDistance'` and `'MaxNumSegments'` aren't yet
- **Spectral functions** take the window as a required argument, and its length is the segment length. Where MATLAB is given an integer n for its window, pass `window.Hamming(n)`, the symmetric Hamming window, which matches MATLAB exactly (verified against R2025a). Without a window, MATLAB and SciPy each choose a window, segment length, overlap and FFT length, and differently, so tukey makes that choice explicit. Its defaults are half-window overlap, an NFFT of the window's length, and no detrending. `Values` are `[segment][frequency]`, the transpose of MATLAB's `S` and `P`, and two-sided frequencies run from 0 to Nyquist and then negative, as `fft.Freq`, where MATLAB labels the upper half above Nyquist. `pwelch` has no median average or per-segment detrending, which `Average(AverageMedian)` and `SegmentDetrend` add
- **Sums** in `Sum`, `Mean`, `CumSum`, `Trapz` and `CumTrapz` are compensated, where MATLAB's are plain, so on long signals, or small values among large ones, tukey is closer to the exact answer: 1e8 among 20,000 thousandths, then cancelled, sums to 20 where MATLAB's `sum` gives 20.00001. `Mean` of an empty signal returns `ErrInvalidLength`, where MATLAB's `mean` gives NaN
- **`Percentile`** interpolates linearly by default, as NumPy does. MATLAB's `prctile` is `PercentileBy(PercentileHazen)`. Any NaN in `x` gives NaN, where `prctile` ignores NaN: use `FillMissing` first to fill the gaps. `Percentile` takes one `p` per call, where `prctile` takes a vector

`Interp1` returns NaN outside the sampled range, as MATLAB's `interp1` and SciPy's `interp1d(bounds_error=False)` do. NumPy's `np.interp` differs there: it repeats the end values. `Interp1` takes vectors, not the matrices MATLAB also accepts, and needs `x` strictly increasing. `FillMissing` likewise works on one vector, so fill a matrix column by column.

`Trapz`, `CumTrapz` and `Grad` take MATLAB's and NumPy's established short names. SciPy's are `trapezoid(y, dx=dx)` and `cumulative_trapezoid(y, dx=dx, initial=0)`, and NumPy's `gradient(x, dx)`. `CumTrapz` starts at 0, the same length as `y`, as MATLAB's `cumtrapz` does: SciPy's `cumulative_trapezoid` leaves the 0 out unless given `initial=0`.

The `fft` package's names read after the package rather than repeating it:

| MATLAB                 | SciPy                        | tukey                                   |
| ---------------------- | ---------------------------- | --------------------------------------- |
| `fft(x)`               | `scipy.fft.fft(x)`           | `fft.Forward(x)`                        |
| `fft(x, n)`            | `scipy.fft.fft(x, n)`        | `fft.Forward(x, fft.N(n))`              |
| `ifft(X)`              | `scipy.fft.ifft(X)`          | `fft.Inverse(X)`                        |
| `fft(x)`, first half   | `scipy.fft.rfft(x)`          | `fft.ForwardReal(x)`                    |
| `ifft(X, 'symmetric')` | `scipy.fft.irfft(X, n)`      | `fft.InverseReal(X, n)` (n is required) |
| `(0:n-1)*fs/n`         | `scipy.fft.fftfreq(n, 1/fs)` | `fft.Freq(n, 1/fs)`                     |
| `fftshift(X)`          | `scipy.fft.fftshift(X)`      | `fft.Shift(X)`                          |
| `ifftshift(X)`         | `scipy.fft.ifftshift(X)`     | `fft.InverseShift(X)`                   |

`fft.Freq` lists the upper half of the bins as negative frequencies, same as SciPy.
MATLAB's usual `(0:n-1)*fs/n` lists them as frequencies above Nyquist. `fft.InverseReal` needs its length because n/2 + 1 bins come from both an even and an odd length.

Windows take their shape as an option, `window.Periodic()`, for MATLAB's `'periodic'` and SciPy's `sym=False`. Most map by name, these differ:

| MATLAB                    | SciPy                               | tukey                                                 |
| ------------------------- | ----------------------------------- | ----------------------------------------------------- |
| `hann(n, 'periodic')`     | `hann(n, sym=False)`                | `window.Hann(n, window.Periodic())`                   |
| `rectwin(n)`              | `boxcar(n)`                         | `window.Boxcar(n)`                                    |
| `nuttallwin(n)`           | `nuttall(n)`                        | `window.Nuttall(n)`                                   |
| `flattopwin(n)`           | `flattop(n)`                        | `window.FlatTop(n)`                                   |
| `tukeywin(n, r)`          | `tukey(n, r)`                       | `window.Tukey(n, r)`                                  |
| `gausswin(n, alpha)`      | `gaussian(n, (n-1)/(2*alpha))`      | `window.Gaussian(n, (n-1)/(2*alpha))`                 |
| `taylorwin(n, nbar, sll)` | `taylor(n, nbar, -sll, norm=False)` | `window.Taylor(n, nbar, -sll, window.Unnormalised())` |
| `hanning(n)`              |                                     | `window.Hanning(n)`                                   |

`gausswin` takes the reciprocal of a width relative to the window's half-length, where SciPy and tukey take the standard deviation in samples. `taylorwin` takes its sidelobe level as negative dB and doesn't normalise; at n = 1 it evaluates its formula, where `window.Taylor` returns `[1]`, as every window does. `hanning` is the Hann window without its zero endpoints.

## Conventions

- Offline and stateless
  - The whole signal is input as a slice and a new slice is returned
  - Inputs never modified
  - All functions are safe for concurrent use
  - For block processing, `FilterState` and `SOSFilterState` take a starting state and return the final one
- Errors and Panics
  - Error for an invalid input, this is a 'user' error
  - Panic for an internal error, this is a 'developer' error
  - Errors wrap an error sentinel so can be checked with `errors.Is`
- Stability checking
  - Filtering functions assume the filter is stable
  - Same as SciPy and MATLAB

## Roadmap

tukey is currently in initial development and the interface will signficantly change until v1.0.

The library is three packages, mirroring SciPy:

| Package        | Mirrors                | For                                                               |
| -------------- | ---------------------- | ----------------------------------------------------------------- |
| `tukey/dsp`    | `scipy.signal`         | Filter design and filtering, spectral estimation, peaks, analysis |
| `tukey/fft`    | `scipy.fft`            | FFTs                                                              |
| `tukey/window` | `scipy.signal.windows` | Window functions                                                  |

| Version | Status   | Theme                  | Adds                                                                                                                                                      |
| ------- | -------- | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| v0.1    | Complete | Zero-phase filtering   | `Butter`, `Filter`, `FiltFilt`                                                                                                                            |
| v0.2    | Complete | Filter representations | Filters as zeros/poles/gain with `(b, a)` and second-order-sections output, SOS filtering, filter state for block processing, frequency response          |
| v0.3    | Complete | FFT                    | `tukey/fft`: complex and real transforms of any length in O(n log n), normalisation and padding options, frequency and shift helpers                      |
| v0.4    | Complete | Windows                | `tukey/window`: every SciPy window except `dpss`, plus MATLAB's `hanning`, symmetric or periodic                                                          |
| v0.5    | Complete | Conditioning           | `Detrend`, `Interp1`, `FillMissing`, `RMS`, checked against SciPy and MATLAB                                                                              |
| v0.6    | Complete | Features               | `FindPeaks`, `PeakProminences`, `PeakWidths`, `Hilbert`, analytic and moving-RMS `Envelope`                                                               |
| v0.7    | Complete | Alignment              | `FindSignal` (template search), exact by default at FFT speed                                                                                             |
| v0.8    | Complete | Spectral analysis      | `STFT` (at FFT bins or chosen frequencies), `Spectrogram`, `Welch`, `Periodogram`                                                                         |
| v0.9    | Complete | Reductions             | `Sum`, `Mean`, `Diff`, `Grad`, `CumSum`, `Trapz`, `CumTrapz`, `Percentile`, with compensated sums                                                         |
| v0.10   | Complete | Filter design breadth  | `Cheby1`, `Cheby2`, `Ellip`, `Bessel`, order selection, `IIRNotch`, `IIRPeak`, `FIRWin`, `GroupDelay`, `Convolve`, `Correlate`, `MedFilt`, `SavGolFilter` |
| v0.11   | Planned  | MATLAB compatibility   | MATLAB-compatibility mode for `SOSFiltFilt`                                                                                                               |
| v0.12   | Planned  | API review             | The whole API reviewed for consistency before v1.0                                                                                                        |

Functions with no SciPy equivalent (`FillMissing`, the moving-RMS `Envelope`, `FindSignal`) are checked against MATLAB, and against MATLAB's documented definition built from SciPy and NumPy.

### Out of scope

- Streaming objects that hold state: tukey stays stateless, and block processing is supported by functions that take and return filter state
- Adaptive filters, wavelets, Parks–McClellan (`remez`) FIR design
- Audio I/O and file formats

## Development

```sh
just test       # go vet + go test (unit, and compare to SciPy and MATLAB reference vectors), on every architecture this machine can run
just tolerances # print measured error per function against SciPy and MATLAB, per architecture
just bench      # FFT benchmarks, six runs each, for comparing with benchstat
just vectors    # regenerate reference vectors: SciPy, then MATLAB from SciPy's inputs
```

Reference vectors are committed in `testdata/scipy/` and `testdata/matlab/`, so neither Python nor MATLAB is needed to run the tests. Their generators are in `reference/`: `reference/scipy` (uv, pinned SciPy) and `reference/matlab`, which reruns each SciPy vector's inputs through MATLAB so both are checked on identical inputs.

`just bectors-matlab` needs MATLAB with the Signal Processing Toolbox, on the `PATH` or set as `MATLAB=/path/to/matlab`.

### How tukey was built

tukey was developed with an LLM (Claude). Design decisions, review and most of the code was written directly by me, written from proposal snippets from the LLM.

All functions are directly tested against results from SciPy equivalents, and from MATLAB where it has one, with published tolerances that can be verified with `just tolerances`.

## Licence

MIT
