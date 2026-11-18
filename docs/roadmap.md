# Roadmap

What tukey has, and what's planned or not.

## Released

| Version | Theme                  | Added                                                                                                                                                     |
| ------- | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| v0.1    | Zero-phase filtering   | `Butter`, `Filter`, `FiltFilt`                                                                                                                            |
| v0.2    | Filter representations | zeros, poles and gain with (b, a) and second-order sections, SOS filtering, filter state for block processing, `Freqz`                                    |
| v0.3    | FFT                    | `tukey/fft`: complex and real transforms of any length in O(n log n), normalisation and length options, frequency and shift helpers                       |
| v0.4    | Windows                | `tukey/window`: every SciPy window except `dpss`, symmetric or periodic                                                                                   |
| v0.5    | Conditioning           | `Detrend`, `Interp1`, `FillMissing`, `RMS`                                                                                                                |
| v0.6    | Features               | `FindPeaks`, `PeakProminences`, `PeakWidths`, `Hilbert`, analytic and moving-RMS `Envelope`                                                               |
| v0.7    | Alignment              | `FindSignal`, a template search, exact by default at FFT speed                                                                                            |
| v0.8    | Spectral analysis      | `STFT`, at the FFT's bins or chosen frequencies, `Spectrogram`, `Welch`, `Periodogram`                                                                    |
| v0.9    | Reductions             | `Sum`, `Mean`, `Diff`, `Grad`, `CumSum`, `Trapz`, `CumTrapz`, `Percentile`, with compensated sums                                                         |
| v0.10   | Filter design breadth  | `Cheby1`, `Cheby2`, `Ellip`, `Bessel`, order selection, `IIRNotch`, `IIRPeak`, `FIRWin`, `GroupDelay`, `Convolve`, `Correlate`, `MedFilt`, `SavGolFilter` |
| v0.11   | MATLAB compatibility   | `BySection` and `MATLABPadLen`, matching MATLAB's `filtfilt(sos, g, x)`                                                                                   |
| v1.0    | Stable API             | code, documentation and API reviewed for consistency, intended to be stable from here                                                                     |

## Planned

Valuable, and planned as demand appears, in no particular order.

- Filter design:
  `FIRWin2`, a FIR design from any frequency response, and `KaiserOrd`, the order for a Kaiser-window design.
- Conversions:
  (b, a) to zeros, poles and gain, and to sections, the same as SciPy's `tf2zpk` and `tf2sos`.
- Resampling:
  `ResamplePoly`, `Decimate` and `UpFirDn`.
- Spectral analysis:
  `ISTFT`, the inverse of `STFT`, and `CSD`, `Coherence` and `LombScargle`.
- MATLAB options not yet supported:
  `findsignal`'s `'Normalization'`, `'Metric'`, `'TimeAlignment'`, `'MaxDistance'` and `'MaxNumSegments'`,
  `xcorr`'s normalisations,
  `findpeaks`' `'SortStr'`, `'NPeaks'` and `'WidthReference'`,
  and `envelope`'s `'peak'` method.
- Performance, to investigate:
  dedicated 7 and 11 point FFT butterflies, so lengths with those factors are as fast as 2, 3 and 5,
  after which `fft.FastLength` could allow them, same as SciPy's default.

## Not currently planned

- Streaming objects that hold state.
  tukey works on whole signals, and `FilterState` and `SOSFilterState` cover block processing.
- Adaptive filters, such as LMS and RLS.
- Wavelets, and FFTs of more than one dimension.
- Parks-McClellan FIR design, `remez`, unless it's asked for.
- Analog filter design output.
- The `dpss` window, which is an eigenvalue problem.
- Audio input and output, file formats, hardware and modulation.
- Checking a filter's stability by default, the way SciPy and MATLAB do by default.
