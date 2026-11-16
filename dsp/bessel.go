package dsp

import (
	"github.com/cxd309/tukey/internal/poly"
)

// Bessel designs a digital Bessel filter for band and returns ZPK.
// Its analog prototype has a maximally flat group delay,
// so a pulse in the passband keeps its shape,
// after the bilinear transform the delay is only close to flat,
// more so the lower the band edge.
// Its transition is the gentlest of the IIR designs for its order.
// By default a band edge is the -3 dB point, same as Butter,
// BesselNormPhase and BesselNormDelay move it.
//
//	f, err := dsp.Bessel(4, dsp.Lowpass(0.3))
//
// Its poles are found exactly, which takes about 2 ms at order 8 and 30 ms at
// order 25, growing as the cube of the order.
//
// Returns:
//   - ErrInvalidOrder if order < 1
//   - ErrInvalidBand for a band Butter rejects
//
// Equivalent to:
//   - scipy.signal.bessel(order, wn, btype, norm='mag', output='zpk'), whose
//     default is norm='phase', BesselNormPhase()
//   - matlab besself(order, 4·tan(π·wn/2), ftype), through
//     bilinear(z, p, k, 2), with BesselNormPhase(), up to
//     order 25. MATLAB has no digital Bessel design.
func Bessel(order int, band Band, opts ...BesselOption) (f ZPK, err error) {
	var cfg besselConfig
	for _, opt := range opts {
		opt(&cfg)
	}
	return designDigital(order, band, func() analogZPK { return besselPrototype(order, cfg.norm) })
}

// BesselOption configures Bessel
type BesselOption func(*besselConfig)

// BesselNormMagnitude puts the band edge at the -3 dB point, same as Butter.
// The default, where SciPy's is BesselNormPhase.
func BesselNormMagnitude() BesselOption {
	return func(c *besselConfig) { c.norm = besselNormMagnitude }
}

// BesselNormPhase scales the prototype so it falls as 1/ω^N at high frequency,
// same as a Butterworth filter of the same edge does,
// and its phase is near -N·π/4 at the edge.
// SciPy's default, and MATLAB's besself.
func BesselNormPhase() BesselOption {
	return func(c *besselConfig) { c.norm = besselNormPhase }
}

// BesselNormDelay gives the prototype a group delay of 1 s at DC,
// so the digital filter's delay at DC is about 1/(π·edge) samples,
// the edge normalised to Nyquist.
func BesselNormDelay() BesselOption {
	return func(c *besselConfig) { c.norm = besselNormDelay }
}

// besselNorm is which property of the analog prototype is fixed at 1 rad/s
// before it's moved onto the band
type besselNorm int

const (
	besselNormMagnitude besselNorm = iota
	besselNormPhase
	besselNormDelay
)

// besselConfig is Bessel's settings, built from BesselOptions
type besselConfig struct {
	norm besselNorm
}

// besselPrototype returns the analog Bessel lowpass prototype:
// poles at the roots of the reverse Bessel polynomial,
// scaled for norm, no finite zeros, and 0 dB at DC.
//
// Equivalent to scipy.signal.besselap(order, norm), but every pole exact,
// rounded once, where SciPy's are within about 1e-15, or 1e-13 for norm="mag".
func besselPrototype(order int, norm besselNorm) (proto analogZPK) {
	proto.Poles = poly.ReverseBesselRoots(order, besselScales[norm])
	proto.Gain = real(poly.Prod(poly.MapRoots(proto.Poles, negate)))
	return
}

// besselScales is each besselNorm's scaling of the reverse Bessel polynomial's roots
var besselScales = [...]poly.ReverseBesselScale{
	besselNormMagnitude: poly.UnitHalfPower,
	besselNormPhase:     poly.UnitProduct,
	besselNormDelay:     poly.UnitDelay,
}
