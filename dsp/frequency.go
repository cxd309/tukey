package dsp

import (
	"fmt"
	"math"
)

// ConvertFrequency returns f, a frequency or a slice of them in the unit from,
// in the unit to.
// fs is the sample rate in Hz, used only when either unit is UnitHz.
// A slice gives a new slice, and a single value a single value:
//
//	wn, err := dsp.ConvertFrequency(50.0, dsp.UnitHz, dsp.UnitNyquist, 1000) // 0.1
//	f, err := dsp.ConvertFrequency(w, dsp.UnitRadPerSample, dsp.UnitHz, fs)
//
// A single value must be a float64, so write 50.0 rather than 50.
// Between UnitHz, UnitNyquist and UnitCyclesPerSample the result is correctly rounded.
// To or from UnitRadPerSample it's within an ULP or two, as π isn't exact.
//
// Returns:
//   - ErrInvalidParameter for a unit that isn't one of the four,
//     or for an fs that isn't positive and finite when either unit is UnitHz.
//
// Equivalent to:
//   - scipy.signal.freqz(b, a, fs=fs), which returns w in Hz
//   - MATLAB freqz(b, a, n, fs), which returns w in Hz
func ConvertFrequency[T float64 | []float64](f T, from, to FrequencyUnit, fs float64) (g T, err error) {
	fromNyquist, err := from.nyquist(fs)
	if err != nil {
		return
	}
	toNyquist, err := to.nyquist(fs)
	if err != nil {
		return
	}
	switch f := any(f).(type) {
	case float64:
		g = any(f * toNyquist / fromNyquist).(T)
	case []float64:
		converted := make([]float64, len(f))
		for i, v := range f {
			converted[i] = v * toNyquist / fromNyquist
		}
		g = any(converted).(T)
	}
	return
}

// FrequencyUnit is a unit of frequency, for ConvertFrequency.
type FrequencyUnit int

const (
	// UnitNyquist is a fraction of the Nyquist frequency, fs/2,
	// the unit of a Band's edges, IIRNotch's w0 and FIRWin's band.
	UnitNyquist FrequencyUnit = iota
	// UnitRadPerSample is radians per sample, where Nyquist is π,
	// the unit of Freqz's and GroupDelay's frequencies.
	UnitRadPerSample
	// UnitCyclesPerSample is cycles per sample, f/fs, where Nyquist is 0.5,
	// the unit of fft.Freq(n, 1).
	UnitCyclesPerSample
	// UnitHz is cycles per second, where Nyquist is fs/2,
	// the unit of the spectral functions' frequencies.
	UnitHz
)

// HzToNyquist returns f in Hz as a fraction of Nyquist, f/(fs/2),
// the unit of a Band's edges:
//
//	wn, err := dsp.HzToNyquist(50.0, 1000) // 0.1
//	f, err := dsp.Butter(4, dsp.Lowpass(wn))
//
// Returns:
//   - ErrInvalidParameter for an fs that isn't positive and finite.
//
// Equivalent to:
//   - the fs argument of scipy.signal.butter and the other designs
//   - MATLAB f/(fs/2)
func HzToNyquist[T float64 | []float64](f T, fs float64) (wn T, err error) {
	wn, err = ConvertFrequency(f, UnitHz, UnitNyquist, fs)
	return
}

// NyquistToHz returns wn, a fraction of Nyquist, in Hz, wn·fs/2.
//
// Returns:
//   - ErrInvalidParameter for an fs that isn't positive and finite.
//
// Equivalent to:
//   - MATLAB wn*fs/2
func NyquistToHz[T float64 | []float64](wn T, fs float64) (f T, err error) {
	f, err = ConvertFrequency(wn, UnitNyquist, UnitHz, fs)
	return
}

// HzToRadPerSample returns f in Hz in radians per sample, 2π·f/fs.
//
// Returns:
//   - ErrInvalidParameter for an fs that isn't positive and finite.
//
// Equivalent to:
//   - MATLAB 2*pi*f/fs
func HzToRadPerSample[T float64 | []float64](f T, fs float64) (w T, err error) {
	w, err = ConvertFrequency(f, UnitHz, UnitRadPerSample, fs)
	return
}

// RadPerSampleToHz returns w in radians per sample in Hz, w·fs/(2π),
// e.g. to label Freqz's response in Hz:
//
//	w, h, err := dsp.Freqz(b, a, 512)
//	f, err := dsp.RadPerSampleToHz(w, fs)
//
// Returns:
//   - ErrInvalidParameter for an fs that isn't positive and finite.
//
// Equivalent to:
//   - scipy.signal.freqz(b, a, fs=fs)'s w
//   - MATLAB freqz(b, a, n, fs)'s w
func RadPerSampleToHz[T float64 | []float64](w T, fs float64) (f T, err error) {
	f, err = ConvertFrequency(w, UnitRadPerSample, UnitHz, fs)
	return
}

// NyquistToRadPerSample returns wn, a fraction of Nyquist, in radians per sample, wn·π.
//
// Equivalent to:
//   - MATLAB wn*pi
func NyquistToRadPerSample[T float64 | []float64](wn T) (w T) {
	w, _ = ConvertFrequency(wn, UnitNyquist, UnitRadPerSample, 0) // can't fail without UnitHz
	return
}

// RadPerSampleToNyquist returns w in radians per sample as a fraction of Nyquist, w/π,
// e.g. to compare Freqz's response with a Band's edges.
//
// Equivalent to:
//   - scipy.signal.freqz(b, a, fs=2)'s w
//   - MATLAB w/pi
func RadPerSampleToNyquist[T float64 | []float64](w T) (wn T) {
	wn, _ = ConvertFrequency(w, UnitRadPerSample, UnitNyquist, 0) // can't fail without UnitHz
	return
}

// nyquist returns the Nyquist frequency in unit u, at sample rate fs
func (u FrequencyUnit) nyquist(fs float64) (nyquist float64, err error) {
	switch u {
	case UnitNyquist:
		nyquist = 1
	case UnitRadPerSample:
		nyquist = math.Pi
	case UnitCyclesPerSample:
		nyquist = 0.5
	case UnitHz:
		if !(fs > 0) || math.IsInf(fs, 1) {
			return 0, fmt.Errorf("%w: fs %g, must be positive and finite", ErrInvalidParameter, fs)
		}
		nyquist = fs / 2
	default:
		return 0, fmt.Errorf("%w: frequency unit %d", ErrInvalidParameter, u)
	}
	return
}
