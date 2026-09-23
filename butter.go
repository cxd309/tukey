package dsp

import (
	"fmt"
	"math"
	"math/cmplx"
)

// fs is the fixed "sampling frequency" SciPy and MATLAB both use internally
// for normalised digital filter design where wn=1.0 means Nyquist
// It has no meaning outside this internal convention
// exists only so that the bilinear transform's frequency warping has a concrete fs to warp against
const (
	fs  = 2.0
	fs2 = 2 * fs // 4.0, used directly by bilinear's (fs2±s)/(fs2∓s) form
)

type BandType int

const (
	LowPass BandType = iota
	HighPass
	BandPass
	BandStop
)

// buttap return the poles of the analog Butterworth prototype lowpass filter of the given order
// normalised to a -3dB cutoff at 1rad/s
// has no zerps and unity gain
func buttap(order int) (poles []complex128) {
	poles = make([]complex128, order)
	for k := 0; k < order; k++ {
		m := float64(-order + 1 + 2*k) // -N+1, -N+3, ...., N-1
		theta := math.Pi * m / (2 * float64(order))
		poles[k] = -cmplx.Exp(complex(0, theta))
	}
	return
}

// prewarp maps a digital cutoff wn (normalised to Nyquist, 1.0 == Nyquist)
// to the analog frequency the bilinear transform should target
// using fixed fs=2 convention shared by SciPy and MATLAB
func prewarp(wn float64) (warped float64) {
	warped = fs2 * math.Tan(math.Pi*wn/fs)
	return
}

// lp21p scales the prototype's poles from a 1 rad/s cutoff to wo
// the prototype has no zeros so the "degree" (numPoles-numZeros) is len(poles)
func lp2lp(poles []complex128, wo float64) (newPoles []complex128, gain float64) {
	newPoles = make([]complex128, len(poles))
	for i, p := range poles {
		newPoles[i] = p * complex(wo, 0)
	}
	gain = math.Pow(wo, float64(len(poles)))
	return
}

// lp2hp converts the prototype to a highpass filter with cutoff wo
// every pole reflects through wo, and a zero appears at the origin for each pole
// (since the prototype start with none)
func lp2hp(poles []complex128, wo float64) (newZeros, newPoles []complex128, gain float64) {
	n := len(poles)
	newZeros = make([]complex128, n) // zero value complex128 is 0+0i - n zeros at the origin
	newPoles = make([]complex128, n)
	prodNegP := complex(1, 0)
	for i, p := range poles {
		newPoles[i] = complex(wo, 0) / p
		prodNegP *= -p
	}
	gain = real(complex(1, 0) / prodNegP)
	return
}

// lp2bp converts the prototype to a bandpass filter centred on a wo with bandwidth bw
// eah prototype pole becomes two (via the quandratic formula)
// and "degree" zeros appear at the origin
func lp2bp(poles []complex128, bw, wo float64) (newZeros, newPoles []complex128, gain float64) {
	degree := len(poles)
	newZeros = make([]complex128, degree) // zero value 0+0i - degree zeros at the origin
	newPoles = make([]complex128, 2*degree)
	for i, p := range poles {
		pLp := complex(0.5*bw, 0) * p
		disc := cmplx.Sqrt(pLp*pLp - complex(wo*wo, 0))
		newPoles[2*i] = pLp + disc
		newPoles[2*i+1] = pLp - disc
	}
	gain = math.Pow(bw, float64(degree))
	return
}

func lp2bs(poles []complex128, bw, wo float64) (newZeros, newPoles []complex128, gain float64) {
	degree := len(poles)
	newZeros = make([]complex128, 2*degree)
	for i := 0; i < degree; i++ {
		newZeros[i] = complex(0, wo)
		newZeros[degree+i] = complex(0, -wo)
	}
	newPoles = make([]complex128, 2*degree)
	prodNegP := complex(1, 0)
	for i, p := range poles {
		pLp := complex(0.5*bw, 0) / p
		disc := cmplx.Sqrt(pLp*pLp - complex(wo*wo, 0))
		newPoles[2*i] = pLp + disc
		newPoles[2*i+1] = pLp - disc
		prodNegP *= -p
	}
	gain = real(complex(1, 0) / prodNegP)
	return
}

// bilinear applies the bilinear transform (fs=2 convention)
// to move an alalog zero/pole/gain filter into the digital domain
// any "zeros at infinity" implied by having fewer zeros than poles map to z=-1
func bilinear(zeros, poles []complex128, gain float64) (digitalZeros, digitalPoles []complex128, digitalGain float64) {
	digitalZeros = make([]complex128, len(poles))
	prodNumber := complex(1, 0)
	for i, z := range zeros {
		digitalZeros[i] = (complex(fs2, 0) + z) / (complex(fs2, 0) - z)
		prodNumber *= complex(fs2, 0) - z
	}
	for i := len(zeros); i < len(poles); i++ {
		digitalZeros[i] = complex(-1, 0)
	}
	digitalPoles = make([]complex128, len(poles))
	prodDenom := complex(1, 0)
	for i, p := range poles {
		digitalPoles[i] = (complex(fs2, 0) + p) / (complex(fs2, 0) - p)
		prodDenom *= complex(fs2, 0) - p
	}

	digitalGain = gain * real(prodNumber/prodDenom)
	return
}

// zpk2tf expands a digital zero/pole/gain filter into
// transfer-function coefficients b (numerator) and a (denominator)
func zpk2tf(zeros, poles []complex128, gain float64) (b, a []float64) {
	bComplex := polyFromRoots(zeros)
	for i := range bComplex {
		bComplex[i] *= complex(gain, 0)
	}
	b = realCoeffs(bComplex)
	a = realCoeffs(polyFromRoots(poles))
	return
}

func Butter(order int, wn float64, bt BandType) (b, a []float64, err error) {
	if order < 1 {
		return nil, nil, fmt.Errorf("dsp: order must be >= 1, got %d", order)
	}
	if wn <= 0 || wn >= 1 {
		return nil, nil, fmt.Errorf("dsp: wn must be in (0,1), got %v", wn)
	}

	proto := buttap(order)
	warped := prewarp(wn)

	var zeros, poles []complex128
	var gain float64
	switch bt {
	case LowPass:
		poles, gain = lp2lp(proto, warped)
	case HighPass:
		zeros, poles, gain = lp2hp(proto, warped)
	default:
		return nil, nil, fmt.Errorf("dsp: Butter does not support band type %v; use ButterBand", bt)
	}

	dz, dp, dk := bilinear(zeros, poles, gain)
	b, a = zpk2tf(dz, dp, dk)
	return
}

func ButterBand(order int, low, high float64, bt BandType) (b, a []float64, err error) {
	if order < 1 {
		return nil, nil, fmt.Errorf("dsp: order must be >= 1, got %d", order)
	}
	if low <= 0 || high >= 1 || low >= high {
		return nil, nil, fmt.Errorf("dsp: require 0 < low < high < 1, got low=%v high=%v", low, high)
	}

	proto := buttap(order)
	warpedLow := prewarp(low)
	warpedHigh := prewarp(high)
	bw := warpedHigh - warpedLow
	wo := math.Sqrt(warpedLow * warpedHigh)

	var zeros, poles []complex128
	var gain float64
	switch bt {
	case BandPass:
		zeros, poles, gain = lp2bp(proto, bw, wo)
	case BandStop:
		zeros, poles, gain = lp2bs(proto, bw, wo)
	default:
		return nil, nil, fmt.Errorf("dsp: ButterBand does not support band type %v; use Butter", bt)
	}

	dz, dp, dk := bilinear(zeros, poles, gain)
	b, a = zpk2tf(dz, dp, dk)

	return
}
