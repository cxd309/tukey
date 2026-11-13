package dsp

import (
	"fmt"

	"github.com/cxd309/tukey/internal/poly"
)

// BesselNorm is where a Bessel filter's band edge is:
// which property of the analog prototype is fixed at 1 rad/s before it's moved onto the band
type BesselNorm int

const (
	// BesselNormMagnitude puts the band edge at the -3 dB point, same as Butter.
	// The default, where SciPy's is BesselNormPhase.
	BesselNormMagnitude BesselNorm = iota
	// BesselNormPhase scales the prototype so it falls as 1/ω^N at high frequency,
	// same as a Butterworth filter of the same edge does,
	// and its phase is near −N·π/4 at the edge.
	// scipy default, and matlab besself.
	BesselNormPhase
	// BesselNormDelay gives the prototype a group delay of 1 s at DC,
	// so the digital filter's delay at DC is about 1/(π·edge) samples,
	// the edge normalised to Nyquist.
	BesselNormDelay
)

// besselScales is each BesselNorm's scaling of the reverse Bessel polynomial's roots
var besselScales = [...]poly.ReverseBesselScale{
	BesselNormMagnitude: poly.UnitHalfPower,
	BesselNormPhase:     poly.UnitProduct,
	BesselNormDelay:     poly.UnitDelay,
}

// BesselOption configures Bessel
type BesselOption func(*besselConfig)

// besselConfig is Bessel's settings, built from BesselOptions
type besselConfig struct {
	norm BesselNorm
}

// BesselNormalised sets where the band edge is. BesselNormMagnitude by default.
func BesselNormalised(norm BesselNorm) BesselOption {
	return func(c *besselConfig) { c.norm = norm }
}

// besselPrototype returns the analog Bessel lowpass prototype:
// poles at the roots of the reverse Bessel polynomial,
// scaled for norm, no finite zeros, and 0 dB at DC.
//
// Equivalent to scipy.signal.besselap(order, norm), but every pole exact,
// rounded once, where SciPy's are within about 1e-15, or 1e-13 for norm="mag".
func besselPrototype(order int, norm BesselNorm) (proto analogZPK) {
	proto.Poles = poly.ReverseBesselRoots(order, besselScales[norm])
	proto.Gain = real(poly.Prod(poly.MapRoots(proto.Poles, negate)))
	return
}

// Bessel designs a digital Bessel filter for band and returns ZPK.
// Its analog prototype has a maximally flat group delay,
// so a pulse in the passband keeps its shape,
// after the bilinear transform the delay is only close to flat,
// more so the lower the band edge.
// Its transition is the gentlest of the IIR designs for its order.
// By default a band edge is the -3 dB point, same as Butter,
// BesselNormalised moves it.
//
//	f, err := dsp.Bessel(4, dsp.Lowpass(0.3))
//
// Its poles are found exactly, which takes about 2 ms at order 8 and 30 ms at
// order 25, growing as the cube of the order.
//
// Returns:
//   - ErrInvalidOrder if order < 1
//   - ErrInvalidBand for a band Butter rejects
//   - ErrInvalidMode for an unknown BesselNorm
//
// Equivalent to:
//   - scipy.signal.bessel(order, wn, btype, norm='mag', output='zpk'), whose
//     default is norm='phase', BesselNormalised(BesselNormPhase)
//   - matlab besself(order, 4·tan(π·wn/2), ftype), through
//     bilinear(z, p, k, 2), with BesselNormalised(BesselNormPhase), up to
//     order 25. MATLAB has no digital Bessel design.
func Bessel(order int, band Band, opts ...BesselOption) (f ZPK, err error) {
	var cfg besselConfig
	for _, opt := range opts {
		opt(&cfg)
	}
	if cfg.norm < BesselNormMagnitude || cfg.norm > BesselNormDelay {
		return f, fmt.Errorf("%w: BesselNorm %d", ErrInvalidMode, cfg.norm)
	}
	return designDigital(order, band, func() analogZPK { return besselPrototype(order, cfg.norm) })
}
