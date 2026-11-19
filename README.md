<h1 align="center"><img src="images/logo.svg" alt="tukey" width="320"></h1>

[![Go Reference](https://pkg.go.dev/badge/github.com/cxd309/tukey.svg)](https://pkg.go.dev/github.com/cxd309/tukey)
[![Test](https://github.com/cxd309/tukey/actions/workflows/test.yaml/badge.svg)](https://github.com/cxd309/tukey/actions/workflows/test.yaml)
[![Lint](https://github.com/cxd309/tukey/actions/workflows/lint.yaml/badge.svg)](https://github.com/cxd309/tukey/actions/workflows/lint.yaml)
[![Coverage](https://codecov.io/gh/cxd309/tukey/graph/badge.svg)](https://codecov.io/gh/cxd309/tukey)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Digital Signal Processing (DSP) library for Go engineers, checked against SciPy and MATLAB.

Named after [John Tukey](https://en.wikipedia.org/wiki/John_Tukey), co-inventor of the FFT and considered by some as the father of Data Science.

## Packages

The library is three packages, mirroring SciPy, each usable on its own:

| Package        | Mirrors                | For                                                               |
| -------------- | ---------------------- | ----------------------------------------------------------------- |
| `tukey/fft`    | `scipy.fft`            | FFTs of any length                                                |
| `tukey/window` | `scipy.signal.windows` | Window functions                                                  |
| `tukey/dsp`    | `scipy.signal`         | Filter design and filtering, spectral estimation, peaks, analysis |

## Usage

```go
import "github.com/cxd309/tukey/fft"

spectrum, err := fft.ForwardReal(x)       // any length, O(n log n)
freqs, err := fft.RealFreq(len(x), 1/fs)  // each bin's frequency in Hz
```

```go
import "github.com/cxd309/tukey/window"

w, err := window.Hann(1024, window.Periodic()) // periodic, for spectral analysis
k, err := window.Kaiser(64, 8.6)               // symmetric, for filter design
```

```go
import "github.com/cxd309/tukey/dsp"

f, err := dsp.Butter(2, dsp.Lowpass(0.3))
b, a, err := f.BA()
y, err := dsp.FiltFilt(b, a, x)

// higher orders and narrow bands: second-order sections
f, err = dsp.Butter(8, dsp.Bandpass(0.30, 0.32))
sos, err := f.SOS()
y, err = dsp.SOSFiltFilt(sos, x)

// the lowest order that meets a specification, then the design
order, band, err := dsp.EllipOrd(dsp.Lowpass(0.2), dsp.Lowpass(0.25), 0.5, 60)
f, err = dsp.Ellip(order, 0.5, 60, band)

y, err = dsp.MedFilt(x, 5)          // spikes removed, steps kept
y, err = dsp.SavGolFilter(x, 21, 3) // noise smoothed, peaks kept

// the window sets the segment length
psd, err := dsp.Welch(x, fs, w)      // power spectral density, half-overlapped segments
tf, err := dsp.Spectrogram(x, fs, w) // power over time: tf.Values[segment][frequency]
```

## Why

In brief: I wanted to port an web application from matlab into my preferred language golang, but there were no existing DSP libraries that could support what I needed.

Therefore I decided to develop this package/set of tools and make it public source so that someone might not encounter the same challenge.

I did consider _"just using scipy"_ but I prefer the tooling and package management that comes with go.

## Documentation

- [Coming from MATLAB](docs/matlab.md)
- [Coming from SciPy](docs/scipy.md)
- [Accuracy](docs/accuracy.md) _measured errors against SciPy, MATLAB and exact answers_
- [Roadmap](docs/roadmap.md) _what's released, and what's planned or not_
- [Development](docs/development.md) _building, testing, reference vectors and conventions_
- [Contributing](CONTRIBUTING.md)

Each package's full documentation is on [pkg.go.dev](https://pkg.go.dev/github.com/cxd309/tukey).

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
- Checked on every change
  - Tests against SciPy and MATLAB reference vectors, on Linux, macOS and Windows
  - `go vet` and `golangci-lint`

## Status

tukey has been released as of v1.0.0 and the API should be considered as stable until v2 which is not currently planned.
See [roadmap](docs/roadmap.md) for planned and unplanned future work.

## Licence

MIT
