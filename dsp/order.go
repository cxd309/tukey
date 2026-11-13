package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/special"
)

// MatchBand is the band edge an order selection meets exactly.
// The smallest order meeting a specification usually has margin to spare,
// and the band it returns decides where the margin goes.
type MatchBand int

const (
	// MatchPassband meets the passband edge exactly, at -rippleDB,
	// and the stopband is attenuated more than asked from its edge.
	// The default, same as scipy buttord, cheb1ord, cheb2ord and ellipord
	MatchPassband MatchBand = iota
	// MatchStopband meets the stopband edge exactly, at -attenuationDB,
	// and the passband reaches further than its edge.
	// Same as matlab buttord and cheb2ord
	MatchStopband
)

// OrderOption configures ButterOrd, Cheby1Ord, Cheby2Ord and EllipOrd
type OrderOption func(*orderConfig)

// orderConfig is the order selection settings, built from OrderOptions
type orderConfig struct {
	match MatchBand
}

// Matching sets the band edge met exactly. Default MatchPassband
func Matching(band MatchBand) OrderOption {
	return func(c *orderConfig) { c.match = band }
}

// orderDesign is the design an order is selected for
type orderDesign int

const (
	orderButter orderDesign = iota
	orderCheby1
	orderCheby2
	orderEllip
)

// maxSelectedOrder bounds a slected order,
// which must fit an int on every platform.
// Far below it, edges this close are better met by an FIR design.
const maxSelectedOrder = math.MaxInt32

// lowpassSpec is a specification as the analog lowpass prototype sees it:
// the passband edge at 1 rad/s and the stopband edge at nar,
// with gpass = 10^(rippleDB/10)-1 and gstop= 10^(attenuationDB/10)- 1,
// so the gains there are 1/(1+ gpass) and 1/(1+ gstop).
type lowpassSpec struct {
	nat, gpass, gstop float64
}

// rawOrder is the real-valued order that meets the specification exactly, the order is its ceiling
func (spec lowpassSpec) rawOrder(design orderDesign) (raw float64) {
	switch design {
	case orderButter:
		raw = math.Log(spec.gstop/spec.gpass) / (2 * math.Log(spec.nat))
	case orderCheby1, orderCheby2:
		raw = math.Acosh(math.Sqrt(spec.gstop/spec.gpass)) / math.Acosh(spec.nat)
	default:
		// the degree equation, N= K(m)·K′(m1) / (K′(m)·K(m1)),
		// with m= 1/nat² and m1 = gpass/gstop,
		// and 1- m= (nat- 1)(nat+ 1)/ nat², which keeps its digits as nat nears 1
		m := 1 / (spec.nat * spec.nat)
		mc := (spec.nat - 1) * (spec.nat + 1) * m
		m1 := spec.gpass / spec.gstop
		raw = special.EllipKM1(mc) * special.EllipKM1(m1) / (special.EllipKM1(m) * special.EllipKM1(1-m1))
	}
	return
}

// prototypeWn is rho, where the order's prototype puts the band edge,
// in the specification's units, so match's edge is met exactly
func (spec lowpassSpec) prototypeWn(design orderDesign, match MatchBand, order int) (rho float64) {
	n := float64(order)
	switch design {
	case orderButter: // Wn is the -3dB point
		if match == MatchPassband {
			rho = math.Exp(-math.Log(spec.gpass) / (2 * n))
		} else {
			rho = spec.nat * math.Exp(-math.Log(spec.gstop)/(2*n))
		}
	case orderCheby1: // Wn is the passband edge
		rho = 1
		if match == MatchStopband {
			rho = spec.nat / spec.selectivity(n)
		}
	case orderCheby2: // Wn is the stopband edge
		rho = spec.nat
		if match == MatchPassband {
			rho = spec.selectivity(n)
		}
	default: // Wn is the passband edge, and the stopband edge is Wn/√m
		rho = 1
		if match == MatchStopband {
			rho = spec.nat * math.Sqrt(ellipticParameter(order, spec.gpass/spec.gstop))
		}
	}
	return
}

// selectivity is how far past the passbnd edge a Chebyshev design of order n
// reaches -attenuationDB, in units of that edge
func (spec lowpassSpec) selectivity(n float64) (ratio float64) {
	ratio = math.Cosh(math.Acosh(math.Sqrt(spec.gstop/spec.gpass)) / n)
	return
}

// ellipticParameter is m, the parameter of an elliptic design of
// order N for the specification's parameter m1, whose stopband edge is 1/√m:
// from the degree equation by the nome, q=q1^(1/N), same as ellipticPrototype
func ellipticParameter(order int, m1 float64) (m float64) {
	k1, k1Prime := special.EllipKM1(1-m1), special.EllipKM1(m1)
	m = special.ParameterFromNome(math.Exp(-math.Pi * k1Prime / (float64(order) * k1)))
	return
}

// checkOrderBands reports whether pass and stop describe one specification:
// the same type of band, each valid, with the stopband edges outside the
// passband's for a lowpass, highpass or bandpass, and inside them for a bandstop
func checkOrderBands(pass, stop Band) (err error) {
	if err = pass.validate(); err != nil {
		return err
	}
	if err = stop.validate(); err != nil {
		return err
	}
	if pass.btype != stop.btype {
		return fmt.Errorf("%w: passband %v and stopband %v are different types of band", ErrInvalidBand, pass, stop)
	}
	p, s := pass.edges, stop.edges
	var ordered bool
	switch pass.btype {
	case lowpassType:
		ordered = p[0] < s[0]
	case highpassType:
		ordered = s[0] < p[0]
	case bandpassType:
		ordered = s[0] < p[0] && p[1] < s[1]
	case bandstopType:
		ordered = p[0] < s[0] && s[1] < p[1]
	}
	if !ordered {
		return fmt.Errorf("%w: passband %v and stopband %v: the transition bands must lie between their edges", ErrInvalidBand, pass, stop)
	}
	return
}

// checkRippleAttenuation reports whether rippleDB and attenuationDB are positive and finite,
// with attenuationDB the greater
func checkRippleAttenuation(rippleDB, attenuationDB float64) (err error) {
	if err = checkDB("rippleDB", rippleDB); err != nil {
		return err
	}
	if err = checkDB("attenuationDB", attenuationDB); err != nil {
		return err
	}
	if attenuationDB <= rippleDB {
		return fmt.Errorf("%w: attenuationDB must be greater than rippleDB, got %g and %g", ErrInvalidParameter, attenuationDB, rippleDB)
	}
	return
}

// warpEdges is tan(π·w/2) of each edge:
// prewarp, without its scale, which cancels in every ratio here
func warpEdges(edges []float64) (warped []float64) {
	warped = make([]float64, len(edges))
	for i, w := range edges {
		warped[i] = math.Tan(math.Pi * w / 2)
	}
	return
}

// symmetricBandstop moves a bandstop's passband edges p inward until p1·p2 = s1·s2,
// so the transformed stopband edges are equal and both bind:
// the lowest order there is.
// The edge on the side with margin moves.
func symmetricBandstop(p, s []float64) (symmetric []float64) {
	product := s[0] * s[1]
	if p[0]*p[1] < product {
		return []float64{product / p[1], p[1]}
	}
	return []float64{p[0], product / p[0]}
}

// stopbandRatio is nat:
// the stopband edge of the equivalent lowpass,
// whose passband edge is 1, from the tighter of two stopband edges
func stopbandRatio(btype bandType, p, s []float64) (nat float64) {
	switch btype {
	case lowpassType:
		return s[0] / p[0]
	case highpassType:
		return p[0] / s[0]
	}
	width, centreSq := p[1]-p[0], p[0]*p[1]
	nat = math.Inf(1)
	for _, x := range s {
		ratio := (x*x - centreSq) / (x * width)
		if btype == bandstopType {
			ratio = 1 / ratio
		}
		nat = min(nat, math.Abs(ratio))
	}
	return
}

// bandEdges is where the equivalent lowpass is at rho, back on the band:
// the inverse of the lowpass-to-band transform.
// A band's two edges multiply to p1·p2, so the lower comes from the upper rather than by a difference.
func bandEdges(btype bandType, p []float64, rho float64) (edges []float64) {
	switch btype {
	case lowpassType:
		return []float64{rho * p[0]}
	case highpassType:
		return []float64{p[0] / rho}
	}
	width, centreSq := p[1]-p[0], p[0]*p[1]
	var upper float64
	if btype == bandpassType {
		upper = math.Hypot(rho*width/2, math.Sqrt(centreSq)) + rho*width/2
	} else {
		upper = (math.Hypot(width, 2*rho*math.Sqrt(centreSq)) + width) / (2 * rho)
	}
	return []float64{centreSq / upper, upper}
}

// selectOrder is the order selection for design, shared by ButterOrd,
// Cheby1Ord, Cheby2Ord and EllipOrd
func selectOrder(design orderDesign, pass, stop Band, rippleDB, attenuationDB float64, opts []OrderOption) (order int, wn Band, err error) {
	var cfg orderConfig
	for _, opt := range opts {
		opt(&cfg)
	}
	if cfg.match < MatchPassband || cfg.match > MatchStopband {
		return 0, wn, fmt.Errorf("%w: MatchBand %d", ErrInvalidMode, cfg.match)
	}
	if err = checkRippleAttenuation(rippleDB, attenuationDB); err != nil {
		return 0, wn, err
	}
	if err = checkOrderBands(pass, stop); err != nil {
		return 0, wn, err
	}

	p, s := warpEdges(pass.edges), warpEdges(stop.edges)
	if pass.btype == bandstopType {
		p = symmetricBandstop(p, s)
	}
	spec := lowpassSpec{
		nat:   stopbandRatio(pass.btype, p, s),
		gpass: math.Expm1(rippleDB * math.Ln10 / 10),
		gstop: math.Expm1(attenuationDB * math.Ln10 / 10),
	}
	raw := spec.rawOrder(design)
	if !(raw <= maxSelectedOrder) {
		return 0, wn, fmt.Errorf("%w: passband %v and stopband %v are too close to meet", ErrInvalidBand, pass, stop)
	}
	order = int(math.Ceil(raw))

	rho := spec.prototypeWn(design, cfg.match, order)
	edges := bandEdges(pass.btype, p, rho)
	for i, x := range edges {
		edges[i] = 2 * math.Atan(x) / math.Pi
	}
	wn = Band{btype: pass.btype, edges: edges}
	return
}

// ButterOrd returns the lowest order of Butter filter that meets a specification,
// and the band to design it for:
// at most rippleDB of loss up to the passband's edges, and at least attenuationDB from the stopband's.
// pass and stop are the same type of band, with stop's edges in the stopband:
//
//	order, band, err := dsp.ButterOrd(dsp.Lowpass(0.2), dsp.Lowpass(0.3), 1, 40)
//	f, err := dsp.Butter(order, band)
//
// The order usually has margin to spare.
// By default the passband edge is met exactly and the stopband gets the margin, Matching(MatchStopband) swaps them.
// A bandpass's tighter stopband edge is the one met.
// A bandstop's passband edges are first moved inward until the stopband's edges are equally hard to meet,
// the lowest order there is.
//
// Returns:
//   - ErrInvalidBand for a band Butter rejects, if pass and stop are different
//     types of band, or if stop's edges aren't beyond pass's
//   - ErrInvalidParameter if rippleDB or attenuationDB isn't positive and finite,
//     or attenuationDB isn't greater than rippleDB
//   - ErrInvalidMode for an unknown MatchBand
//
// Equivalent to:
//   - scipy.signal.buttord(wp, ws, rippleDB, attenuationDB), whose bandstop
//     edges come from a numerical minimiser, to about 1e-5
//   - matlab buttord(wp, ws, rippleDB, attenuationDB) with Matching(MatchStopband),
//     with the same minimiser for a bandstop
func ButterOrd(pass, stop Band, rippleDB, attenuationDB float64, opts ...OrderOption) (order int, wn Band, err error) {
	return selectOrder(orderButter, pass, stop, rippleDB, attenuationDB, opts)
}

// Cheby1Ord returns the lowest order of Cheby1 filter that meets a specification,
// and the band to design it for, as ButterOrd:
//
//	order, band, err := dsp.Cheby1Ord(dsp.Lowpass(0.2), dsp.Lowpass(0.3), 1, 40)
//	f, err := dsp.Cheby1(order, 1, band)
//
// Returns the errors ButterOrd does.
//
// Equivalent to:
//   - scipy.signal.cheb1ord(wp, ws, rippleDB, attenuationDB), whose bandstop
//     edges come from a numerical minimiser, to about 1e-5
//   - matlab cheb1ord(wp, ws, rippleDB, attenuationDB), except for a bandstop,
//     whose order MATLAB selects for the moved passband edges but returns
//     the edges unmoved, which can fall short of attenuationDB
func Cheby1Ord(pass, stop Band, rippleDB, attenuationDB float64, opts ...OrderOption) (order int, wn Band, err error) {
	return selectOrder(orderCheby1, pass, stop, rippleDB, attenuationDB, opts)
}

// Cheby2Ord returns the lowest order of Cheby2 filter that meets a specification,
// and the band to design it for, same as ButterOrd.
// Cheby2's band edges are its stopband's,
// so by default those returned lie inside the stopband, where -attenuationDB is first reached:
//
//	order, band, err := dsp.Cheby2Ord(dsp.Lowpass(0.2), dsp.Lowpass(0.3), 1, 40)
//	f, err := dsp.Cheby2(order, 40, band)
//
// Returns the errors ButterOrd does.
//
// Equivalent to:
//   - scipy.signal.cheb2ord(wp, ws, rippleDB, attenuationDB), whose bandstop
//     edges come from a numerical minimiser, to about 1e-5
//   - matlab cheb2ord(wp, ws, rippleDB, attenuationDB) with Matching(MatchStopband),
//     except for a bandpass, where MATLAB returns the stopband edges, both
//     met exactly, and tukey only the tighter
func Cheby2Ord(pass, stop Band, rippleDB, attenuationDB float64, opts ...OrderOption) (order int, wn Band, err error) {
	return selectOrder(orderCheby2, pass, stop, rippleDB, attenuationDB, opts)
}

// EllipOrd returns the lowest order of Ellip filter that meets a specification,
// and the band to design it for, as ButterOrd:
//
//	order, band, err := dsp.EllipOrd(dsp.Lowpass(0.2), dsp.Lowpass(0.3), 1, 40)
//	f, err := dsp.Ellip(order, 1, 40, band)
//
// Returns the errors ButterOrd does.
//
// Equivalent to:
//   - scipy.signal.ellipord(wp, ws, rippleDB, attenuationDB), whose bandstop
//     edges come from a numerical minimiser, to about 1e-5
//   - matlab ellipord(wp, ws, rippleDB, attenuationDB), except for a bandstop,
//     where MATLAB returns the passband edges unmoved
func EllipOrd(pass, stop Band, rippleDB, attenuationDB float64, opts ...OrderOption) (order int, wn Band, err error) {
	return selectOrder(orderEllip, pass, stop, rippleDB, attenuationDB, opts)
}
