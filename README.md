# tukey

Offline Digital Signal Processing (DSP) for GO, matching SciPy and MATLAB functionality.

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
import "github.com/cxd309/tukey/fft"

spectrum, err := fft.ForwardReal(x)  // any length, O(n log n)
freqs := fft.RealFreq(len(x), 1/fs)  // each bin's frequency in Hz
```

```go
import "github.com/cxd309/tukey/window"

w, err := window.Hann(len(x), window.Periodic()) // periodic for spectral analysis
w, err = window.Kaiser(64, 8.6)                  // symmetric, for filter design
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

| Function       | SciPy equivalent                 | Worst measured | Tested to |
| -------------- | -------------------------------- | -------------- | --------- |
| `Butter`       | `butter`, as `(b, a)`            | 7.6e-15        | 1e-13     |
| `Butter` + SOS | `butter(output='sos')`           | 7.3e-15        | 1e-13     |
| `ZPK.SOS`      | `zpk2sos`                        | 1.1e-16        | 1e-14     |
| `Filter`       | `lfilter`                        | 6.4e-16        | 1e-14     |
| `FilterState`  | `lfilter` with `zi`              | 2.3e-16        | 1e-14     |
| `FilterZi`     | `lfilter_zi`                     | 3.3e-16        | 1e-14     |
| `FiltFilt`     | `filtfilt`                       | 8.8e-16        | 1e-14     |
| `SOSFilter`    | `sosfilt`, with and without `zi` | 4.1e-14        | 1e-12     |
| `SOSFilterZi`  | `sosfilt_zi`                     | 9.8e-17        | 1e-14     |
| `SOSFiltFilt`  | `sosfiltfilt`                    | 1.7e-13        | 1e-11     |
| `Freqz`        | `freqz`                          | 1.4e-15        | 1e-13     |
| `SOSFreqz`     | `freqz_sos`                      | 3.7e-13        | 1e-11     |
| `Detrend`      | `detrend`                        | 5.6e-16        | 1e-14     |
| `Interp1`      | `interp1d(bounds_error=False)`   | 1.3e-16        | 1e-14     |
| `RMS`          | `sqrt(mean(x**2))`, from NumPy   | 3.1e-16        | 1e-14     |

For scale, 1e-16 is ~ 1 ULP of a `float64`. The larger SOS figures all come from one deliberately difficult test filter (order 8, cutoff 0.01), whose poles sit so close to z = 1 that single rounding differences between Go and SciPy are amplified. Every other filter agrees within a few ULP.

Two very low cutoff filters, order 2 at 0.001 low-pass and 0.0008 high-pass, are checked separately. They're the original MATLAB application's own designs, and their poles sit even closer to z = 1, so rounding differences between implementations are magnified by roughly the square of 1/(1 − pole):

| Function      | Worst against SciPy | Tested to | Worst against MATLAB | Tested to |
| ------------- | ------------------- | --------- | -------------------- | --------- |
| `Filter`      | 9.5e-13             | 1e-11     | 6.2e-13              | 1e-11     |
| `FilterState` | 9.9e-14             | 1e-12     | 1.1e-13              | 1e-11     |
| `FiltFilt`    | 4.4e-12             | 1e-10     | 1.5e-11              | 1e-9      |
| `Freqz`       | 1.5e-12             | 1e-10     | 1.1e-11              | 1e-9      |

SciPy and MATLAB differ from each other on these filters by as much. At DC the 0.001 low-pass's gain is exactly 1, and `Freqz` gives 1 + 1.6e-12 where MATLAB gives 1 + 1.3e-11, so much of the difference is MATLAB's.

Functions MATLAB also has are checked against vectors generated by MATLAB R2025a (Signal Processing Toolbox 25.1), on the same inputs:

| Function                  | MATLAB equivalent               | Worst measured | Tested to |
| ------------------------- | ------------------------------- | -------------- | --------- |
| `Butter`                  | `[b, a] = butter`               | 1.2e-13        | 1e-11     |
| `Filter`                  | `filter`                        | 6.4e-16        | 1e-14     |
| `FilterState`             | `[y, zf] = filter(b, a, x, zi)` | 3.0e-16        | 1e-14     |
| `FiltFilt` + MATLABPadLen | `filtfilt(b, a, x)`             | 8.8e-16        | 1e-14     |
| `SOSFilter`               | `sosfilt`                       | 4.5e-14        | 1e-12     |
| `Freqz`                   | `freqz(b, a, n)`                | 1.9e-15        | 1e-13     |
| `SOSFreqz`                | `freqz(sos, n)`                 | 1.1e-12        | 1e-10     |
| `Detrend`                 | `detrend`                       | 1.2e-15        | 1e-13     |
| `Interp1`                 | `interp1`                       | 2.2e-16        | 1e-14     |
| `FillMissing`             | `fillmissing`                   | exact          | 1e-14     |
| `RMS`                     | `rms`                           | 3.1e-16        | 1e-14     |

The largest come from the same difficult filter, where MATLAB and SciPy differ from each other by as much. SciPy has no `fillmissing`, so `FillMissing` is checked against MATLAB's documented definition as well as MATLAB itself. See [Porting from MATLAB](#porting-from-matlab) for the functions not in this table.

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

| MATLAB                          | tukey                                                                      |
| ------------------------------- | -------------------------------------------------------------------------- |
| `[b, a] = butter(n, Wn)`        | `f, err := dsp.Butter(n, dsp.Lowpass(Wn))`, `b, a, err := f.BA()`          |
| `butter(n, Wn, 'high')`         | `dsp.Butter(n, dsp.Highpass(Wn))`                                          |
| `butter(n, [W1 W2])`            | `dsp.Butter(n, dsp.Bandpass(W1, W2))`                                      |
| `butter(n, [W1 W2], 'stop')`    | `dsp.Butter(n, dsp.Bandstop(W1, W2))`                                      |
| `[z, p, k] = butter(n, Wn)`     | `f, err := dsp.Butter(n, dsp.Lowpass(Wn))`: `f.Zeros`, `f.Poles`, `f.Gain` |
| `y = filter(b, a, x)`           | `y, err := dsp.Filter(b, a, x)`                                            |
| `y = filtfilt(b, a, x)`         | `y, err := dsp.FiltFilt(b, a, x, dsp.MATLABPadLen())`                      |
| `[h, w] = freqz(b, a, n)`       | `w, h, err := dsp.Freqz(b, a, n)` (note the order)                         |
| `y = sosfilt(sos, x)`           | `y, err := dsp.SOSFilter(sos, x)`                                          |
| `y = detrend(x)`                | `y, err := dsp.Detrend(x, dsp.DetrendLinear)`                              |
| `y = detrend(x, 'constant')`    | `y, err := dsp.Detrend(x, dsp.DetrendConstant)`                            |
| `yi = interp1(x, y, xi)`        | `yi, err := dsp.Interp1(x, y, xi)`                                         |
| `fillmissing(x, 'linear')`      | `y := dsp.FillMissing(x, dsp.FillLinear())`                                |
| `fillmissing(x, 'constant', v)` | `y := dsp.FillMissing(x, dsp.FillConstant(v))`                             |
| `r = rms(x)`                    | `r, err := dsp.RMS(x)`                                                     |

Frequencies are normalised to Nyquist in both: for a frequency `f`Hz at sample rate `fs`, use `f / (fs/2)`. `Freqz` returns frequencies in rad/sample, as MATLAB's `freqz` does.

`MATLABPadLen` uses MATLAB's `filtfilt` padding length, `3*(nfilt-1)`, instead of SciPy's `3*nfilt`. Only samples near the ends of the output are affected.

Where tukey and MATLAB differ, verified against MATLAB R2025a:

- **Second-order sections:** `f.SOS()` pairs and orders sections as SciPy's `zpk2sos` does, with the overall gain in the first section. MATLAB's `zp2sos` gives a first-order section one zero rather than two, and breaks ties in pole radius differently: the same filter, in different sections. `[sos, g] = zp2sos(z, p, k)` also returns the gain `g` separately
- **`SOSFiltFilt` with `MATLABPadLen`** does not match MATLAB's `filtfilt(sos, g, x)`, which filters forward then backward one section at a time, each padded by 3 × the total filter order. Expect differences near the ends of the output; a MATLAB-compatible mode is on the roadmap
- **`Freqz` with n = 0** returns empty slices, as SciPy does, where MATLAB's `freqz` errors
- **`RMS` of an empty signal** returns `ErrInvalidLength`, where MATLAB's `rms` gives NaN

`Interp1` returns NaN outside the sampled range, as MATLAB's `interp1` and SciPy's `interp1d(bounds_error=False)` do. NumPy's `np.interp` differs there: it repeats the end values. `Interp1` takes vectors, not the matrices MATLAB also accepts, and needs `x` strictly increasing. `FillMissing` likewise works on one vector, so fill a matrix column by column.

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

| Version | Status   | Theme                  | Adds                                                                                                                                             |
| ------- | -------- | ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| v0.1    | Complete | Zero-phase filtering   | `Butter`, `Filter`, `FiltFilt`                                                                                                                   |
| v0.2    | Complete | Filter representations | Filters as zeros/poles/gain with `(b, a)` and second-order-sections output, SOS filtering, filter state for block processing, frequency response |
| v0.3    | Complete | FFT                    | `tukey/fft`: complex and real transforms of any length in O(n log n), normalisation and padding options, frequency and shift helpers             |
| v0.4    | Complete | Windows                | `tukey/window`: every SciPy window except `dpss`, plus MATLAB's `hanning`, symmetric or periodic                                                 |
| v0.5    | Complete | Conditioning           | `Detrend`, `Interp1`, `FillMissing`, `RMS`, checked against SciPy and MATLAB                                                                     |
| v0.6    | Planned  | Features               | `FindPeaks`, a moving-RMS `Envelope`                                                                                                             |
| v0.7    | Planned  | Alignment              | `FindSignal` (template search)                                                                                                                   |
| tbd     | Concept  | Signal analysis        | Spectral estimation (`Welch`, `Spectrogram`), `Hilbert`, further reductions                                                                      |
| tbd     | Concept  | Filter design breadth  |                                                                                                                                                  |
| tbd     | Planned  | MATLAB compatibility   | MATLAB-compatibility mode for `SOSFiltFilt`                                                                                                      |

Two functions have no SciPy equivalent and will be checked against MATLAB instead: a moving-RMS envelope, and `FindSignal` (template search).

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
