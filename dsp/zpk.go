package dsp

import (
	"math"
	"math/cmplx"
)

// fs is the fixed "sampling frequency" SciPy and MATLAB use internally for
// normalised digital filter design, where wn=1.0 means Nyquist
// it has no meaning outside this convention; it only gives the bilinear
// transform's frequency warping a concrete fs to warp against
const (
	fs  = 2.0
	fs2 = 2 * fs // 4.0, used directly in bilinear's (fs2+r)/(fs2-r)
)

// prewarp maps a digital cutoff wn (normalised to Nyquist, 1.0 == Nyquist)
// to the analog frequency the bilinear transform should target
// using fixed fs=2 convention shared by SciPy and MATLAB
func prewarp(wn float64) (warped float64) {
	warped = fs2 * math.Tan(math.Pi*wn/fs)
	return
}

// zpk is a filter described by its zeros, poles and gain
// in either the s-plane (analog) or the z-plane (digital)
type zpk struct {
	zeros []complex128
	poles []complex128
	gain  float64
}

// relativeDegree is how many more poles than zeros the filter has
// each surplus pole is balanced by a zero "at infinity",
// which frequency transforms and bilinear have to place somewhere finite
func (f zpk) relativeDegree() (degree int) {
	degree = len(f.poles) - len(f.zeros)
	return
}

// toLowpass moves a 1 rad/s lowpass prototype to a cutoff of wo rad/s
// by substituting s -> s/wo, which scales every zero and pole by wo
// equivalent to scipy.signal.lp2lp_zpk
func (f zpk) toLowpass(wo float64) (lowpass zpk) {
	scale := func(r complex128) complex128 { return complex(wo, 0) * r }
	lowpass.zeros = mapRoots(f.zeros, scale)
	lowpass.poles = mapRoots(f.poles, scale)
	lowpass.gain = f.gain * math.Pow(wo, float64(f.relativeDegree()))
	return
}

// toHighpass turns a 1 rad/s lowpass prototype into a highpass with cutoff wo rad/s
// by substituting s -> wo/s, which inverts every zero and pole through wo
// each zero at infinity moves to the origin (DC)
// equivalent to scipy.signal.lp2hp_zpk
func (f zpk) toHighpass(wo float64) (highpass zpk) {
	invert := func(r complex128) complex128 { return complex(wo, 0) / r }

	highpass.zeros = appendRepeated(mapRoots(f.zeros, invert), 0, f.relativeDegree())
	highpass.poles = mapRoots(f.poles, invert)
	highpass.gain = f.gain * real(prod(mapRoots(f.zeros, negate))/prod(mapRoots(f.poles, negate)))
	return
}

// toBandpass turns a 1 rad/s lowpass prototype into a bandpass centred on wo rad/s
// with bandwidth bw rad/s by substituting s -> (s² + wo²)/(s·bw)
// that's quadratic in s, so every zero and pole becomes a pair
// each zero at infinity becomes one zero at the origin and one that stays at infinity
// equivalent to scipy.signal.lp2bp_zpk
func (f zpk) toBandpass(wo, bw float64) (bandpass zpk) {
	scale := func(r complex128) complex128 { return complex(bw/2, 0) * r }
	degree := f.relativeDegree()

	bandpass.zeros = appendRepeated(splitRoots(mapRoots(f.zeros, scale), wo), 0, degree)
	bandpass.poles = splitRoots(mapRoots(f.poles, scale), wo)
	bandpass.gain = f.gain * math.Pow(bw, float64(degree))
	return
}

// toBandstop turns a 1 rad/s lowpass prototype into a bandstop centred on wo rad/s
// with bandwidth bw rad/s by substituting s -> (s·bw)/(s² + wo²)
// every zero and pole becomes a pair, and each zero at infinity becomes
// a pair of zeros at ±j·wo, which is the notch
// equivalent to scipy.signal.lp2bs_zpk
func (f zpk) toBandstop(wo, bw float64) (bandstop zpk) {
	invert := func(r complex128) complex128 { return complex(bw/2, 0) / r }
	degree := f.relativeDegree()

	bandstop.zeros = splitRoots(mapRoots(f.zeros, invert), wo)
	bandstop.zeros = appendRepeated(bandstop.zeros, complex(0, wo), degree)
	bandstop.zeros = appendRepeated(bandstop.zeros, complex(0, -wo), degree)
	bandstop.poles = splitRoots(mapRoots(f.poles, invert), wo)
	bandstop.gain = f.gain * real(prod(mapRoots(f.zeros, negate))/prod(mapRoots(f.poles, negate)))
	return
}

// bilinear maps an analog filter to the digital domain via s = fs2*(z-1)/(z+1)
// each root r becomes (fs2+r)/(fs2-r), and each zero at infinity lands at z=-1 (Nyquist)
// equivalent to scipy.signal.bilinear_zpk with fs=2
func (f zpk) bilinear() (digital zpk) {
	toZ := func(r complex128) complex128 { return (fs2 + r) / (fs2 - r) }
	distance := func(r complex128) complex128 { return fs2 - r }

	digital.zeros = appendRepeated(mapRoots(f.zeros, toZ), -1, f.relativeDegree())
	digital.poles = mapRoots(f.poles, toZ)
	digital.gain = f.gain * real(prod(mapRoots(f.zeros, distance))/prod(mapRoots(f.poles, distance)))
	return
}

// transferFunction expands the filter into transfer-function coefficients
// b (numerator/feedforward) and a (denominator/feedback), highest power first
// equivalent to scipy.signal.zpk2tf
func (f zpk) transferFunction() (b, a []float64) {
	bComplex := polyFromRoots(f.zeros)
	for i := range bComplex {
		bComplex[i] *= complex(f.gain, 0)
	}
	b = realCoeffs(bComplex)
	a = realCoeffs(polyFromRoots(f.poles))
	return
}

// negate returns -r, used for the prod(-z)/prod(-p) gain corrections
func negate(r complex128) (negated complex128) {
	negated = -r
	return
}

// splitRoots replaces each root r with the pair r ± sqrt(r²-wo²),
// the two solutions of x² - 2rx + wo² = 0 that the bandpass and bandstop
// substitutions reduce to for each root
// all the "+" roots come first, then all the "-" roots, matching SciPy ordering
func splitRoots(roots []complex128, wo float64) (split []complex128) {
	n := len(roots)
	split = make([]complex128, 2*n)
	for i, r := range roots {
		disc := cmplx.Sqrt(r*r - complex(wo*wo, 0))
		split[i] = r + disc
		split[n+i] = r - disc
	}
	return
}
