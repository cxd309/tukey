package dsp

import (
	"fmt"
	"math"
	"math/cmplx"
)

// butterworthPrototype returns the analog Butterworth lowpass prototype of the given order:
// -3dB at 1 rad/s, no finite zeros, unity gain,
// and poles spaced evenly around the left half of the unit circle
// equivalent to scipy.signal.buttap
func butterworthPrototype(order int) (proto zpk) {
	proto.poles = make([]complex128, order)
	for k := range order {
		m := float64(-order + 1 + 2*k) // -N+1, -N+3, ..., N-1
		theta := math.Pi * m / (2 * float64(order))
		proto.poles[k] = -cmplx.Exp(complex(0, theta))
	}
	proto.gain = 1
	return
}

// Butter designs a digital lowpass or highpass Butterworth filter and returns
// its transfer-function coefficients b (numerator) and a (denominator)
//
// order must be >= 1. wn is the -3dB cutoff normalised to Nyquist, so 0 < wn < 1;
// for a cutoff fc in Hz at sample rate fs, wn = fc / (fs/2)
//
// bt must be LowPass or HighPass; use ButterBand for BandPass and BandStop
// equivalent to scipy.signal.butter(order, wn, btype)
func Butter(order int, wn float64, bt BandType) (b, a []float64, err error) {
	if order < 1 {
		return nil, nil, fmt.Errorf("dsp: order must be >= 1, got %d", order)
	}
	if wn <= 0 || wn >= 1 {
		return nil, nil, fmt.Errorf("dsp: wn must be in (0,1), got %v", wn)
	}

	proto := butterworthPrototype(order)
	warped := prewarp(wn)

	var analog zpk
	switch bt {
	case LowPass:
		analog = proto.toLowpass(warped)
	case HighPass:
		analog = proto.toHighpass(warped)
	case BandPass, BandStop:
		return nil, nil, fmt.Errorf("dsp: Butter does not support %v; use ButterBand", bt)
	default:
		return nil, nil, fmt.Errorf("dsp: unknown band type %v", bt)
	}

	b, a = analog.bilinear().transferFunction()
	return
}

// ButterBand designs a digital bandpass or bandstop Butterworth filter and returns
// its transfer-function coefficients b (numerator) and a (denominator)
//
// order must be >= 1; the resulting filter has order 2*order, since each
// prototype pole becomes a pair. low and high are the -3dB band edges normalised
// to Nyquist, so 0 < low < high < 1
//
// bt must be BandPass or BandStop; use Butter for LowPass and HighPass
// equivalent to scipy.signal.butter(order, [low, high], btype)
func ButterBand(order int, low, high float64, bt BandType) (b, a []float64, err error) {
	if order < 1 {
		return nil, nil, fmt.Errorf("dsp: order must be >= 1, got %d", order)
	}
	if low <= 0 || high >= 1 || low >= high {
		return nil, nil, fmt.Errorf("dsp: require 0 < low < high < 1, got low=%v high=%v", low, high)
	}

	proto := butterworthPrototype(order)
	warpedLow := prewarp(low)
	warpedHigh := prewarp(high)
	bw := warpedHigh - warpedLow
	wo := math.Sqrt(warpedLow * warpedHigh)

	var analog zpk
	switch bt {
	case BandPass:
		analog = proto.toBandpass(wo, bw)
	case BandStop:
		analog = proto.toBandstop(wo, bw)
	case LowPass, HighPass:
		return nil, nil, fmt.Errorf("dsp: ButterBand does not support %v; use Butter", bt)
	default:
		return nil, nil, fmt.Errorf("dsp: unknown band type %v", bt)
	}

	b, a = analog.bilinear().transferFunction()
	return
}
