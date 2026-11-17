package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/special"
)

// ButterOrd returns the lowest order of Butter filter that meets a specification, and the band to design it for:
// at most rippleDB of loss up to the passband's edges, and at least attenuationDB from the stopband's.
// pass and stop are the same type of band, with stop's edges in the stopband.
// The order usually has margin to spare.
// By default the passband edge is met exactly and the stopband gets the margin, and MatchStopband swaps them.
// A bandpass's tighter stopband edge is the one met.
// A bandstop's passband edges are first moved inward until the stopband's edges are equally hard to meet,
// which gives the lowest order there is.
//
//	order, wn, err := dsp.ButterOrd(dsp.Lowpass(0.2), dsp.Lowpass(0.3), 1, 40)
//	f, err := dsp.Butter(order, wn)
//
// Units:
//   - pass's and stop's edges, and wn's, are fractions of Nyquist
//   - rippleDB and attenuationDB are in dB, positive
//   - order is the prototype's order
//
// Returns:
//   - ErrInvalidBand for the zero Band, an edge outside (0, 1), or a low edge not below the high,
//     for pass and stop of different types, for stop's edges not beyond pass's,
//     or for edges too close together to meet.
//   - ErrInvalidParameter for a rippleDB or attenuationDB that isn't positive and finite,
//     or an attenuationDB that isn't greater than rippleDB.
//
// Equivalent to:
//   - scipy.signal.buttord(wp, ws, rippleDB, attenuationDB), whose bandstop edges come from a numerical minimiser, to about 1e-5
//   - MATLAB buttord(wp, ws, rippleDB, attenuationDB), with MatchStopband, and the same minimiser for a bandstop
func ButterOrd(pass, stop Band, rippleDB, attenuationDB float64, opts ...OrderOption) (order int, wn Band, err error) {
	return selectOrder(orderButter, pass, stop, rippleDB, attenuationDB, opts)
}

// Cheby1Ord returns the lowest order of Cheby1 filter that meets a specification, and the band to design it for:
// at most rippleDB of loss up to the passband's edges, and at least attenuationDB from the stopband's.
// pass and stop are the same type of band, with stop's edges in the stopband.
// By default the passband edge is met exactly and the stopband gets the margin, and MatchStopband swaps them.
// A bandpass's tighter stopband edge is the one met,
// and a bandstop's passband edges are first moved inward until its stopband's edges are equally hard to meet.
//
//	order, wn, err := dsp.Cheby1Ord(dsp.Lowpass(0.2), dsp.Lowpass(0.3), 1, 40)
//	f, err := dsp.Cheby1(order, 1, wn)
//
// Units:
//   - pass's and stop's edges, and wn's, are fractions of Nyquist
//   - rippleDB and attenuationDB are in dB, positive
//   - order is the prototype's order
//
// Returns:
//   - ErrInvalidBand for the zero Band, an edge outside (0, 1), or a low edge not below the high,
//     for pass and stop of different types, for stop's edges not beyond pass's,
//     or for edges too close together to meet.
//   - ErrInvalidParameter for a rippleDB or attenuationDB that isn't positive and finite,
//     or an attenuationDB that isn't greater than rippleDB.
//
// Equivalent to:
//   - scipy.signal.cheb1ord(wp, ws, rippleDB, attenuationDB), whose bandstop edges come from a numerical minimiser, to about 1e-5
//   - MATLAB cheb1ord(wp, ws, rippleDB, attenuationDB), except for a bandstop,
//     whose order MATLAB selects for the moved passband edges but whose edges it returns unmoved,
//     which can fall short of attenuationDB
func Cheby1Ord(pass, stop Band, rippleDB, attenuationDB float64, opts ...OrderOption) (order int, wn Band, err error) {
	return selectOrder(orderCheby1, pass, stop, rippleDB, attenuationDB, opts)
}

// Cheby2Ord returns the lowest order of Cheby2 filter that meets a specification, and the band to design it for:
// at most rippleDB of loss up to the passband's edges, and at least attenuationDB from the stopband's.
// pass and stop are the same type of band, with stop's edges in the stopband.
// Cheby2's band edges are its stopband's, so by default those returned lie inside the stopband,
// where -attenuationDB is first reached, and MatchStopband puts them at stop's edges.
// A bandpass's tighter stopband edge is the one met,
// and a bandstop's passband edges are first moved inward until its stopband's edges are equally hard to meet.
//
//	order, wn, err := dsp.Cheby2Ord(dsp.Lowpass(0.2), dsp.Lowpass(0.3), 1, 40)
//	f, err := dsp.Cheby2(order, 40, wn)
//
// Units:
//   - pass's and stop's edges, and wn's, are fractions of Nyquist
//   - rippleDB and attenuationDB are in dB, positive
//   - order is the prototype's order
//
// Returns:
//   - ErrInvalidBand for the zero Band, an edge outside (0, 1), or a low edge not below the high,
//     for pass and stop of different types, for stop's edges not beyond pass's,
//     or for edges too close together to meet.
//   - ErrInvalidParameter for a rippleDB or attenuationDB that isn't positive and finite,
//     or an attenuationDB that isn't greater than rippleDB.
//
// Equivalent to:
//   - scipy.signal.cheb2ord(wp, ws, rippleDB, attenuationDB), whose bandstop edges come from a numerical minimiser, to about 1e-5
//   - MATLAB cheb2ord(wp, ws, rippleDB, attenuationDB), with MatchStopband, except for a bandpass,
//     where MATLAB meets both stopband edges exactly, and Cheby2Ord only the tighter
func Cheby2Ord(pass, stop Band, rippleDB, attenuationDB float64, opts ...OrderOption) (order int, wn Band, err error) {
	return selectOrder(orderCheby2, pass, stop, rippleDB, attenuationDB, opts)
}

// EllipOrd returns the lowest order of Ellip filter that meets a specification, and the band to design it for:
// at most rippleDB of loss up to the passband's edges, and at least attenuationDB from the stopband's.
// pass and stop are the same type of band, with stop's edges in the stopband.
// By default the passband edge is met exactly and the stopband gets the margin, and MatchStopband swaps them.
// A bandpass's tighter stopband edge is the one met,
// and a bandstop's passband edges are first moved inward until its stopband's edges are equally hard to meet.
//
//	order, wn, err := dsp.EllipOrd(dsp.Lowpass(0.2), dsp.Lowpass(0.3), 1, 40)
//	f, err := dsp.Ellip(order, 1, 40, wn)
//
// Units:
//   - pass's and stop's edges, and wn's, are fractions of Nyquist
//   - rippleDB and attenuationDB are in dB, positive
//   - order is the prototype's order
//
// Returns:
//   - ErrInvalidBand for the zero Band, an edge outside (0, 1), or a low edge not below the high,
//     for pass and stop of different types, for stop's edges not beyond pass's,
//     or for edges too close together to meet.
//   - ErrInvalidParameter for a rippleDB or attenuationDB that isn't positive and finite,
//     or an attenuationDB that isn't greater than rippleDB.
//
// Equivalent to:
//   - scipy.signal.ellipord(wp, ws, rippleDB, attenuationDB), whose bandstop edges come from a numerical minimiser, to about 1e-5
//   - MATLAB ellipord(wp, ws, rippleDB, attenuationDB), except for a bandstop, whose passband edges MATLAB returns unmoved
func EllipOrd(pass, stop Band, rippleDB, attenuationDB float64, opts ...OrderOption) (order int, wn Band, err error) {
	return selectOrder(orderEllip, pass, stop, rippleDB, attenuationDB, opts)
}

// OrderOption configures ButterOrd, Cheby1Ord, Cheby2Ord and EllipOrd: MatchPassband or MatchStopband.
type OrderOption func(*orderConfig)

// MatchPassband meets the passband edge exactly, at -rippleDB,
// so the stopband is attenuated more than asked from its edge.
// The lowest order meeting a specification usually has margin to spare,
// and the band returned decides where it goes.
// It's the default.
//
// Equivalent to:
//   - scipy.signal.buttord, cheb1ord, cheb2ord and ellipord, which have no other
//   - MATLAB cheb1ord and ellipord
func MatchPassband() OrderOption {
	return func(c *orderConfig) { c.match = matchPassband }
}

// MatchStopband meets the stopband edge exactly, at -attenuationDB,
// so the passband reaches further than its edge.
//
// Equivalent to:
//   - SciPy has none
//   - MATLAB buttord and cheb2ord
func MatchStopband() OrderOption {
	return func(c *orderConfig) { c.match = matchStopband }
}

// matchBand is the band edge an order selection meets exactly.
type matchBand int

const (
	matchPassband matchBand = iota
	matchStopband
)

// orderConfig is an order selection's settings, built from OrderOptions.
type orderConfig struct {
	match matchBand
}

// selectOrder is the order selection for design, shared by ButterOrd, Cheby1Ord, Cheby2Ord and EllipOrd.
func selectOrder(design orderDesign, pass, stop Band, rippleDB, attenuationDB float64, opts []OrderOption) (order int, wn Band, err error) {
	var cfg orderConfig
	for _, opt := range opts {
		opt(&cfg)
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

// orderDesign is the design an order is selected for.
type orderDesign int

const (
	orderButter orderDesign = iota
	orderCheby1
	orderCheby2
	orderEllip
)

// checkRippleAttenuation returns ErrInvalidParameter unless rippleDB and attenuationDB are positive and finite,
// with attenuationDB the greater.
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

// checkOrderBands returns ErrInvalidBand unless pass and stop describe one specification:
// the same type of band, each valid, with the stopband edges outside the passband's for a lowpass, highpass or bandpass,
// and inside them for a bandstop.
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

// warpEdges is tan(π·w/2) of each edge:
// prewarp, without its scale, which cancels in every ratio here.
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

// lowpassSpec is a specification as the analog lowpass prototype sees it:
// the passband edge at 1 rad/s and the stopband edge at nat,
// with gpass = 10^(rippleDB/10) - 1 and gstop = 10^(attenuationDB/10) - 1,
// so the power gains there are 1/(1 + gpass) and 1/(1 + gstop).
type lowpassSpec struct {
	nat, gpass, gstop float64
}

// rawOrder is the real-valued order that meets the specification exactly.
// The order is its ceiling.
func (spec lowpassSpec) rawOrder(design orderDesign) (raw float64) {
	switch design {
	case orderButter:
		raw = math.Log(spec.gstop/spec.gpass) / (2 * math.Log(spec.nat))
	case orderCheby1, orderCheby2:
		raw = math.Acosh(math.Sqrt(spec.gstop/spec.gpass)) / math.Acosh(spec.nat)
	default:
		// the degree equation, N = K(m)·K′(m1)/(K′(m)·K(m1)),
		// with m = 1/nat² and m1 = gpass/gstop,
		// and 1 - m = (nat - 1)(nat + 1)/nat², which keeps its digits as nat nears 1
		m := 1 / (spec.nat * spec.nat)
		mc := (spec.nat - 1) * (spec.nat + 1) * m
		m1 := spec.gpass / spec.gstop
		raw = special.EllipKM1(mc) * special.EllipKM1(m1) / (special.EllipKM1(m) * special.EllipKM1(1-m1))
	}
	return
}

// prototypeWn is rho, where the order's prototype puts the band edge, in the specification's units,
// so match's edge is met exactly.
func (spec lowpassSpec) prototypeWn(design orderDesign, match matchBand, order int) (rho float64) {
	n := float64(order)
	switch design {
	case orderButter: // Wn is the -3 dB point
		if match == matchPassband {
			rho = math.Exp(-math.Log(spec.gpass) / (2 * n))
		} else {
			rho = spec.nat * math.Exp(-math.Log(spec.gstop)/(2*n))
		}
	case orderCheby1: // Wn is the passband edge
		rho = 1
		if match == matchStopband {
			rho = spec.nat / spec.selectivity(n)
		}
	case orderCheby2: // Wn is the stopband edge
		rho = spec.nat
		if match == matchPassband {
			rho = spec.selectivity(n)
		}
	default: // Wn is the passband edge, and the stopband edge is Wn/√m
		rho = 1
		if match == matchStopband {
			rho = spec.nat * math.Sqrt(ellipticParameter(order, spec.gpass/spec.gstop))
		}
	}
	return
}

// selectivity is how far past the passband edge a Chebyshev design of order n reaches -attenuationDB,
// in units of that edge.
func (spec lowpassSpec) selectivity(n float64) (ratio float64) {
	ratio = math.Cosh(math.Acosh(math.Sqrt(spec.gstop/spec.gpass)) / n)
	return
}

// maxSelectedOrder bounds a selected order, which must fit an int on every platform.
// Far below it, edges this close are better met by an FIR design.
const maxSelectedOrder = math.MaxInt32

// ellipticParameter is m, the parameter of an elliptic design of order N for the specification's parameter m1,
// whose stopband edge is 1/√m.
// It's from the degree equation by the nome, q = q1^(1/N), same as ellipticPrototype.
func ellipticParameter(order int, m1 float64) (m float64) {
	k1, k1Prime := special.EllipKM1(1-m1), special.EllipKM1(m1)
	m = special.ParameterFromNome(math.Exp(-math.Pi * k1Prime / (float64(order) * k1)))
	return
}

// stopbandRatio is nat, the stopband edge of the equivalent lowpass whose passband edge is 1,
// from the tighter of two stopband edges.
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
