package dsp

import (
	"math"
	"math/cmplx"
)

// analogZPK is a filter described by its zeros, poles and gain
// in either the s-plane (analog) or the z-plane (digital)
type analogZPK ZPK

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

// relativeDegree is how many more poles than zeros the filter has
// each surplus pole is balanced by a zero "at infinity",
// which frequency transforms and bilinear have to place somewhere finite
func (f analogZPK) relativeDegree() (degree int) {
	degree = len(f.Poles) - len(f.Zeros)
	return
}

// toLowpass moves a 1 rad/s lowpass prototype to a cutoff of wo rad/s
// by substituting s -> s/wo, which scales every zero and pole by wo
// equivalent to scipy.signal.lp2lp_zpk
func (f analogZPK) toLowpass(wo float64) (lowpass analogZPK) {
	scale := func(r complex128) complex128 { return complex(wo, 0) * r }
	lowpass.Zeros = mapRoots(f.Zeros, scale)
	lowpass.Poles = mapRoots(f.Poles, scale)
	lowpass.Gain = f.Gain * math.Pow(wo, float64(f.relativeDegree()))
	return
}

// toHighpass turns a 1 rad/s lowpass prototype into a highpass with cutoff wo rad/s
// by substituting s -> wo/s, which inverts every zero and pole through wo
// each zero at infinity moves to the origin (DC)
// equivalent to scipy.signal.lp2hp_zpk
func (f analogZPK) toHighpass(wo float64) (highpass analogZPK) {
	invert := func(r complex128) complex128 { return complex(wo, 0) / r }

	highpass.Zeros = appendRepeated(mapRoots(f.Zeros, invert), 0, f.relativeDegree())
	highpass.Poles = mapRoots(f.Poles, invert)
	highpass.Gain = f.Gain * real(prod(mapRoots(f.Zeros, negate))/prod(mapRoots(f.Poles, negate)))
	return
}

// toBandpass turns a 1 rad/s lowpass prototype into a bandpass centred on wo rad/s
// with bandwidth bw rad/s by substituting s -> (s² + wo²)/(s·bw)
// that's quadratic in s, so every zero and pole becomes a pair
// each zero at infinity becomes one zero at the origin and one that stays at infinity
// equivalent to scipy.signal.lp2bp_zpk
func (f analogZPK) toBandpass(wo, bw float64) (bandpass analogZPK) {
	scale := func(r complex128) complex128 { return complex(bw/2, 0) * r }
	degree := f.relativeDegree()

	bandpass.Zeros = appendRepeated(splitRoots(mapRoots(f.Zeros, scale), wo), 0, degree)
	bandpass.Poles = splitRoots(mapRoots(f.Poles, scale), wo)
	bandpass.Gain = f.Gain * math.Pow(bw, float64(degree))
	return
}

// toBandstop turns a 1 rad/s lowpass prototype into a bandstop centred on wo rad/s
// with bandwidth bw rad/s by substituting s -> (s·bw)/(s² + wo²)
// every zero and pole becomes a pair, and each zero at infinity becomes
// a pair of zeros at ±j·wo, which is the notch
// equivalent to scipy.signal.lp2bs_zpk
func (f analogZPK) toBandstop(wo, bw float64) (bandstop analogZPK) {
	invert := func(r complex128) complex128 { return complex(bw/2, 0) / r }
	degree := f.relativeDegree()

	bandstop.Zeros = splitRoots(mapRoots(f.Zeros, invert), wo)
	bandstop.Zeros = appendRepeated(bandstop.Zeros, complex(0, wo), degree)
	bandstop.Zeros = appendRepeated(bandstop.Zeros, complex(0, -wo), degree)
	bandstop.Poles = splitRoots(mapRoots(f.Poles, invert), wo)
	bandstop.Gain = f.Gain * real(prod(mapRoots(f.Zeros, negate))/prod(mapRoots(f.Poles, negate)))
	return
}

// bilinear maps an analog filter to the digital domain via s = fs2*(z-1)/(z+1)
// each root r becomes (fs2+r)/(fs2-r), and each zero at infinity lands at z=-1 (Nyquist)
// equivalent to scipy.signal.bilinear_zpk with fs=2
func (f analogZPK) bilinear() (digital ZPK) {
	toZ := func(r complex128) complex128 { return (fs2 + r) / (fs2 - r) }
	distance := func(r complex128) complex128 { return fs2 - r }

	digital.Zeros = appendRepeated(mapRoots(f.Zeros, toZ), -1, f.relativeDegree())
	digital.Poles = mapRoots(f.Poles, toZ)
	digital.Gain = f.Gain * real(prod(mapRoots(f.Zeros, distance))/prod(mapRoots(f.Poles, distance)))
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
