// Package dsp is offline signal processing for engineers:
// filter design, filtering, frequency response, conditioning, smoothing, peaks,
// alignment, spectral analysis and reductions.
// It works on whole signals as []float64, uses only the standard library,
// and is checked against SciPy and MATLAB, which the README's compatibility tables measure.
//
//	f, err := dsp.Butter(4, dsp.Lowpass(0.2))
//	sos, err := f.SOS()
//	y, err := dsp.SOSFiltFilt(sos, x)
//
// # Filter design
//
// Each design returns a ZPK, whose SOS method gives second-order sections and BA method (b, a).
// A Band is Lowpass, Highpass, Bandpass or Bandstop, with edges as fractions of Nyquist.
//
//	tukey                SciPy                MATLAB
//	Butter               butter               butter
//	Cheby1               cheby1               cheby1
//	Cheby2               cheby2               cheby2
//	Ellip                ellip                ellip
//	Bessel               bessel               besself, then bilinear
//	ButterOrd            buttord              buttord
//	Cheby1Ord            cheb1ord             cheb1ord
//	Cheby2Ord            cheb2ord             cheb2ord
//	EllipOrd             ellipord             ellipord
//	IIRNotch             iirnotch             iirnotch, DSP System Toolbox
//	IIRPeak              iirpeak              iirpeak, DSP System Toolbox
//	FIRWin               firwin               fir1
//	ZPK.BA               zpk2tf               zp2tf
//	ZPK.SOS              zpk2sos              zp2sos
//
// # Filtering and response
//
//	tukey                SciPy                MATLAB
//	Filter               lfilter              filter
//	FilterState          lfilter with zi      filter with zi
//	FilterZi             lfilter_zi
//	SOSFilter            sosfilt              sosfilt
//	SOSFilterState       sosfilt with zi
//	SOSFilterZi          sosfilt_zi
//	FiltFilt             filtfilt             filtfilt
//	SOSFiltFilt          sosfiltfilt          filtfilt(sos, g, x)
//	Freqz                freqz                freqz
//	SOSFreqz             freqz_sos            freqz(sos, n)
//	GroupDelay           group_delay          grpdelay
//	SOSGroupDelay                             grpdelay(sos, n)
//
// ConvertFrequency, and shortcuts such as HzToNyquist and RadPerSampleToHz,
// convert between Hz, rad/sample, cycles/sample and fractions of Nyquist.
//
// # Signals
//
//	tukey                SciPy                         MATLAB
//	Detrend              detrend                       detrend
//	Interp1              interp1d(bounds_error=False)  interp1
//	FillMissing                                        fillmissing
//	Convolve             convolve                      conv
//	Correlate            correlate, correlation_lags   xcorr
//	MedFilt              medfilt                       medfilt1
//	SavGolFilter         savgol_filter                 sgolayfilt
//	FindPeaks            find_peaks                    findpeaks
//	PeakProminences      peak_prominences              findpeaks fourth output
//	PeakWidths           peak_widths                   findpeaks third output
//	Hilbert              hilbert                       hilbert
//	Envelope             envelope                      envelope
//	FindSignal                                         findsignal
//
// # Spectral analysis
//
// Each takes the window as a slice, from tukey/window or anywhere, and its length is the segment length.
//
//	tukey                SciPy                MATLAB
//	STFT                 stft                 spectrogram S
//	Spectrogram          spectrogram          spectrogram P
//	Welch                welch                pwelch
//	Periodogram          periodogram          periodogram
//
// # Reductions and calculus
//
//	tukey                NumPy and SciPy               MATLAB
//	Sum                  numpy.sum                     sum
//	Mean                 numpy.mean                    mean
//	RMS                                                rms
//	Percentile           numpy.percentile              prctile, with PercentileHazen
//	Diff                 numpy.diff                    diff
//	Grad                 numpy.gradient                gradient
//	CumSum               numpy.cumsum                  cumsum
//	Trapz                trapezoid                     trapz
//	CumTrapz             cumulative_trapezoid          cumtrapz
//
// # Conventions
//
// Every function is offline and stateless:
// the whole signal goes in as a slice, and a new slice comes back.
// Inputs are never modified, and nothing is shared between calls, so all functions are safe for concurrent use.
// To filter in blocks, FilterState and SOSFilterState take a starting state and return the final one,
// to pass to the next block.
// An empty result is an empty slice, never nil.
//
// Invalid input returns an error wrapping one of the Err sentinels, so check for one with errors.Is.
// Filtering functions don't check that a filter is stable, and assume it is, same as SciPy and MATLAB.
//
// # Where tukey differs
//
// Where SciPy and MATLAB differ, tukey follows SciPy, and MATLAB's behaviour is an option.
// Where neither suits general use, tukey chooses its own default and says why.
//
// Frequencies:
//   - band edges, w0 and cutoffs are fractions of Nyquist, same as SciPy and MATLAB,
//     and HzToNyquist converts from Hz
//   - Freqz and GroupDelay return w in rad/sample, and w first, where MATLAB returns the response first
//   - the spectral functions work in Hz, given fs
//
// Filter design and filtering:
//   - (b, a) loses digits as the order rises or the band narrows,
//     so use SOS, SOSFilter and SOSFiltFilt above order 4 or so
//   - Bessel's band edge is the -3 dB point by default, same as Butter's, where SciPy's default is BesselNormPhase
//   - the order selections meet the passband edge, same as SciPy,
//     and MatchStopband meets the stopband's, same as MATLAB's buttord and cheb2ord
//   - IIRNotch and IIRPeak take the quality factor q, where MATLAB takes the bandwidth w0/q
//   - FIRWin takes the window itself, whose length is the number of taps, where fir1 takes the order, one less
//   - FiltFilt pads by SciPy's length, and MATLABPadLen uses MATLAB's
//   - SOSFiltFilt filters the whole cascade forward then backward, same as SciPy,
//     and BySection filters each section in turn, same as MATLAB
//
// Signals:
//   - indices are 0-based, in FindPeaks, PeakProminences, PeakWidths and FindSignal
//   - ConvolveSame starts one output earlier than MATLAB's conv 'same' for an even len(v),
//     and ConvolveValid with v longer than a isn't empty
//   - Correlate gives only the lags where a and v overlap, where xcorr pads the shorter with zeros
//   - MedFilt truncates its window at the ends by default, where SciPy and MATLAB pad with zeros, MedFiltZeroPad
//   - FindPeaks follows SciPy's rules, and MATLABPeaks MATLAB's,
//     whose MinPeakDistance d is MinPeakDistance(d+1)
//   - FindSignal's distance is exact, where MATLAB's rounds, and SearchFast rounds the same way as MATLAB
//
// Spectral analysis:
//   - there's no default window or segment length, as together they set the trade-off between time and frequency resolution
//   - there's no detrending by default, where SciPy's welch, spectrogram and periodogram subtract each segment's mean,
//     so pass SegmentDetrend(DetrendConstant()) to match
//   - TimeFrequency's Values are [segment][frequency], the transpose of SciPy's and MATLAB's
//   - STFT isn't divided by the window's sum, where SciPy's stft is
//   - TwoSided frequencies run from 0 up then negative, same as fft.Freq, where MATLAB's run up to fs
//
// Reductions:
//   - sums are compensated, so a long signal, or small values among large ones, keeps its low bits,
//     e.g. 1e8, 20,000 values of 0.001 and -1e8 sum to 20, where MATLAB's sum gives 20.00001
//   - Mean and RMS of an empty signal return ErrSignalTooShort, where MATLAB's mean and rms give NaN
//   - Percentile is linear by default, same as NumPy, and PercentileHazen is MATLAB's prctile,
//     but any NaN gives NaN, where prctile ignores it, so use FillMissing first
//
// # NaN and infinity
//
// A NaN or infinite sample affects only the outputs it contributes to, as plain arithmetic would.
// A sample is part of the signal, so it's never an error.
// What it contributes to depends on the function:
//   - a running sum, a filter with feedback, or FiltFilt's backward pass carries it on
//   - a transform, such as Hilbert, Envelope's analytic method or ConvolveFFT, mixes every sample into every output
//   - Detrend and Envelope's RMS method use the whole signal's trend or mean
//
// A sample and a filter coefficient of exactly 0 contribute nothing,
// so a zero coefficient never turns a NaN or infinite sample into NaN.
// Interpolating next to an infinite value gives that infinity, its limit.
// Each function's doc says where it differs, e.g. Percentile, where any NaN gives NaN.
//
// A NaN or infinite parameter, coefficient, state or window tap returns an error,
// unless the function's doc gives it a meaning, e.g. MinPeakHeight(math.Inf(-1)) is no limit.
package dsp
