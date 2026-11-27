package waveform

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/special"
)

// GaussianPulse is a Guassian-modulated sinusoid, sampled at each time in t
type GaussianPulse struct {
	// InPhase is the pulse, the envelope times cos(2π·fc·t)
	InPhase []float64
	// Quadrature is the envelope times sin(2π·fc·t), 90° behind InPhase
	Quadrature []float64
	// Envelope is exp(-a·t²), peaking at 1 at t= 0
	Envelope []float64
}

// GaussPulse returns a Guassian-modulated sinusoid centred on t= 0 at each sample time in t,
// with the centre frequency fc, its in-phase part, its quadrature part, and its envelope.
// The envelope is exp(-a·t²), where a sets the bandwidth,
// so the spectrum falls to BandwidthLevel, -6dB by default, at fc·(1± bw/ 2),
// with Bandwidth bw, 0.5 by default.
//
//	pulse, err:= waveform.GaussPulse(t, 50e3, waveform.Bandwidth(0.6)) // 50kHz ultrasound pulse
//
// The envelope is computed in double-double, and the carrier's phase reduced to its fraction of a cycle first,
// so each sample is within 1ULP of exact.
// A NaN time gives NaN at that sample, and an infinite time gives 0, the limit of the pulse.
//
// Units:
//   - t seconds
//   - fc hertz, above 0
//   - outputs are dimensionless, from -1 to 1
//
// Returns:
//   - ErrInvalidParameter for:
//     an fc that isn't above 0 and finite,
//     a Bandwidth that isn't above 0 and finite, or,
//     a BandwidthLevel that isn't below 0 and finite.
//   - ErrInvalidOption for an option not in Bandwidth or BandwidthLevel
//
// Equivalent to:
//   - scipy.signal.gausspulse(t, fc, bw, bwr, retquad= True, retenv= True),
//     which takes fc= 0 too, and gives NaN at an infinite time
//   - MATLAB [yi, yq, ye]= gausepuls(t, fc, bw, bwr)
func GaussPulse(t []float64, fc float64, opts ...Option) (pulse GaussianPulse, err error) {
	cfg, err := newConfig("GaussPulse", opts, optionSet{optionBandwidth: true, optionBandwidthLevel: true})
	if err != nil {
		return GaussianPulse{}, err
	}
	shape, err := newGaussShape("GaussPulse", fc, cfg)
	if err != nil {
		return GaussianPulse{}, err
	}
	a := shape.rate()
	pulse = GaussianPulse{
		InPhase:    make([]float64, len(t)),
		Quadrature: make([]float64, len(t)),
		Envelope:   make([]float64, len(t)),
	}
	for i, ti := range t {
		pulse.InPhase[i], pulse.Quadrature[i], pulse.Envelope[i] = gaussPulseAt(ti, fc, a)
	}
	return
}

// GaussPulseCutoff returns the time after t= 0 where GaussPulse's envelope falls to CutoffLevel,
// -60dB by default, for the same fc, Bandwidth and BandwidthLevel.
// The pulse is symmetric, so it's below that level outside -tc to tc.
//
//	tc, err:= waveform.GaussPulseCutoff(50e3, waveform.Bandwidth(0.6), waveform.CutoffLevel(-40))
//
// It's tc= ln(10)· √(tpr· bwr)/ (10π· fc· bw), with tpr the CutoffLevel and bwr the BandwidthLevel,
// computed in double-double, so it's within an ULP of exact.
//
// Units:
//   - fc hertz, above 0
//   - tc seconds
//
// Returns:
//   - ErrInvalidParameter for:
//     an fc that isn't above 0 and finite,
//     a Bandwidth that isn't above 0 and finite,
//     a BandwidthLevel that isn't below 0 and finite, or,
//     a CutoffLevel that isn't below 0 and finite.
//   - ErrInvalidOption for an option not in Bandwidth, BandwidthLevel or CutoffLevel
//
// Equivalent to:
//   - scipy.signal.gausspulse("cutoff", fc, bw, bwr, tpr)
//   - MATLAB gauspuls("cutoff", fc, bw, bwr, tpr)
func GaussPulseCutoff(fc float64, opts ...Option) (tc float64, err error) {
	accepted := optionSet{optionBandwidth: true, optionBandwidthLevel: true, optionCutoffLevel: true}
	cfg, err := newConfig("GaussPulseCutoff", opts, accepted)
	if err != nil {
		return 0, err
	}
	shape, err := newGaussShape("GaussPulseCutoff", fc, cfg)
	if err != nil {
		return 0, err
	}
	cutoffLevel := -60.0
	if cfg.given[optionCutoffLevel] {
		cutoffLevel = cfg.cutoffLevel
	}
	if !(cutoffLevel < 0) || math.IsInf(cutoffLevel, -1) {
		return 0, fmt.Errorf("%w: GaussPulseCutoff's CutoffLevel is %g dB, must be below 0 and finite",
			ErrInvalidParameter, cutoffLevel)
	}
	// √(tpr·bwr) in double-double: the float64 root, then one Newton step
	product := special.ProductOf(cutoffLevel, shape.level)
	root := math.Sqrt(product.Hi)
	exactRoot := special.Exact(root).Add(product.Sub(special.ProductOf(root, root)).Div(special.Exact(2 * root)))

	width := ddPi.Mul(special.ProductOf(fc, shape.bandwidth)).Mul(special.Exact(10))
	tc = ddLn10.Mul(exactRoot).Div(width).Float64()
	return
}

// π and ln(10) as double-doubles, from mpmath.
var (
	ddPi   = special.DoubleDouble{Hi: 3.141592653589793, Lo: 1.2246467991473532e-16}
	ddLn10 = special.DoubleDouble{Hi: 2.302585092994046, Lo: -2.1707562233822494e-16}
)

// gaussShape is a Gaussian pulse's centre frequency and bandwidth
type gaussShape struct {
	fc        float64 // hertz
	bandwidth float64 // fraction of fc
	level     float64 // dB, where the spectrum is bandwidth wide
}

// newGaussShape check fc and the bandwidth options, filling in their defaults
func newGaussShape(waveform string, fc float64, cfg config) (shape gaussShape, err error) {
	shape = gaussShape{fc: fc, bandwidth: 0.5, level: -6}
	if cfg.given[optionBandwidth] {
		shape.bandwidth = cfg.bandwidth
	}
	if cfg.given[optionBandwidthLevel] {
		shape.level = cfg.bandwidthLevel
	}
	switch {
	case !(fc > 0) || math.IsInf(fc, 1):
		err = fmt.Errorf("%w: %s's fc is %g, must be greater than 0 and finite", ErrInvalidParameter, waveform, fc)
	case !(shape.bandwidth > 0) || math.IsInf(shape.bandwidth, 1):
		err = fmt.Errorf("%w: %s's Bandwidth is %g, must be greater than 0 and finite", ErrInvalidParameter, waveform, shape.bandwidth)
	case !(shape.level < 0) || math.IsInf(shape.level, -1):
		err = fmt.Errorf("%w: %s's BandwidthLevel is %g dB, must be below 0 and finite", ErrInvalidParameter, waveform, shape.level)
	}
	return
}

// rate is a, where the envelope is exp(-a·t²): a= 5· (π· fc· bw)²/ (-bwr· ln(10)).
// The spectrum of exp(-a· t²) is a Guassian falling as exp(-π²· f²/ a),
// which is bwr dB, 10^ (bwr/ 20), at f= fc· bw/ 2
func (shape gaussShape) rate() (a special.DoubleDouble) {
	w := ddPi.Mul(special.ProductOf(shape.fc, shape.bandwidth))
	a = w.Mul(w).Mul(special.Exact(5)).Div(ddLn10.Mul(special.Exact(-shape.level)))
	return
}

// gaussPulseAt is the pulse at t, for a centre frequency fc and an envelope exp(-a·t²)
func gaussPulseAt(t, fc float64, a special.DoubleDouble) (inPhase, quadrature, envelope float64) {
	switch {
	case math.IsNaN(t):
		return math.NaN(), math.NaN(), math.NaN()
	case math.IsInf(t, 0):
		return 0, 0, 0
	}
	e := special.Exp(a.Mul(special.ProductOf(t, t)).Neg())
	sin, cos := sinCosCycles(special.ProductOf(fc, t))
	inPhase = e.Mul(special.Exact(cos)).Float64()
	quadrature = e.Mul(special.Exact(sin)).Float64()
	envelope = e.Float64()
	return
}
