package dsp

import (
	"fmt"
	"math"
	"math/big"
	"testing"

	"github.com/cxd309/tukey/internal/special"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// orderSelections are the four order selections, each with the design it
// selects for, given the specification's dB
var orderSelections = []struct {
	name        string
	design      orderDesign
	selectOrder func(pass, stop Band, rippleDB, attenuationDB float64, opts ...OrderOption) (int, Band, error)
	filter      func(order int, rippleDB, attenuationDB float64, band Band) (ZPK, error)
}{
	{"ButterOrd", orderButter, ButterOrd, func(order int, _, _ float64, band Band) (ZPK, error) { return Butter(order, band) }},
	{"Cheby1Ord", orderCheby1, Cheby1Ord, func(order int, rippleDB, _ float64, band Band) (ZPK, error) { return Cheby1(order, rippleDB, band) }},
	{"Cheby2Ord", orderCheby2, Cheby2Ord, func(order int, _, attenuationDB float64, band Band) (ZPK, error) {
		return Cheby2(order, attenuationDB, band)
	}},
	{"EllipOrd", orderEllip, EllipOrd, Ellip},
}

// orderSpecs are specifications of each band type. The asymmetric bandpass
// has one stopband edge much tighter than the other, and the asymmetric
// bandstop moves its upper passband edge from 0.8 to about 0.4, where
// MATLAB's cheb1ord leaves it and misses attenuationDB. 3 dB against 6 dB
// selects order 1 for every design.
var orderSpecs = []struct {
	name                  string
	pass, stop            Band
	rippleDB, attenuation float64
}{
	{"lowpass", Lowpass(0.2), Lowpass(0.3), 1, 40},
	{"lowpass sharp", Lowpass(0.3), Lowpass(0.32), 0.1, 60},
	{"lowpass order 1", Lowpass(0.1), Lowpass(0.5), 3, 6},
	{"highpass", Highpass(0.3), Highpass(0.2), 1, 40},
	{"bandpass", Bandpass(0.2, 0.5), Bandpass(0.1, 0.6), 1, 40},
	{"bandpass asymmetric", Bandpass(0.2, 0.3), Bandpass(0.1, 0.5), 0.5, 50},
	{"bandstop", Bandstop(0.1, 0.6), Bandstop(0.2, 0.5), 1, 40},
	{"bandstop asymmetric", Bandstop(0.1, 0.8), Bandstop(0.2, 0.3), 1, 40},
}

// edgeGains is f's gain in dB at each edge, the lowest first
func edgeGains(f ZPK, edges []float64) (lowest, highest float64) {
	lowest, highest = math.Inf(1), math.Inf(-1)
	for _, w := range edges {
		gain := gainDB(f, w)
		lowest, highest = min(lowest, gain), max(highest, gain)
	}
	return
}

// The selected order, designed for the selected band, meets the specification:
// at most rippleDB of loss at every passband edge and at least attenuationDB at
// every stopband edge, with the matched band's worst edge met exactly.
func TestOrderMeetsSpecification(t *testing.T) {
	for _, sel := range orderSelections {
		for _, spec := range orderSpecs {
			for _, match := range []MatchBand{MatchPassband, MatchStopband} {
				name := fmt.Sprintf("%s, %s, match %d", sel.name, spec.name, match)
				order, band, err := sel.selectOrder(spec.pass, spec.stop, spec.rippleDB, spec.attenuation, Matching(match))
				require.NoError(t, err, name)
				f, err := sel.filter(order, spec.rippleDB, spec.attenuation, band)
				require.NoError(t, err, name)

				passLowest, _ := edgeGains(f, spec.pass.edges)
				_, stopHighest := edgeGains(f, spec.stop.edges)
				if match == MatchPassband {
					assert.InDelta(t, -spec.rippleDB, passLowest, gainSlackDB, name)
					assert.Less(t, stopHighest, -spec.attenuation, name)
				} else {
					assert.InDelta(t, -spec.attenuation, stopHighest, gainSlackDB, name)
					assert.Greater(t, passLowest, -spec.rippleDB, name)
				}
			}
		}
	}
}

// One order lower can't meet the specification: designed to meet the
// passband exactly, it falls short of attenuationDB at a stopband edge.
func TestOrderIsTheLowest(t *testing.T) {
	for _, sel := range orderSelections {
		for _, spec := range orderSpecs {
			name := fmt.Sprintf("%s, %s", sel.name, spec.name)
			order, _, err := sel.selectOrder(spec.pass, spec.stop, spec.rippleDB, spec.attenuation)
			require.NoError(t, err, name)
			if order == 1 {
				continue
			}

			// the passband-matched band for order − 1, as selectOrder finds it
			p, s := warpEdges(spec.pass.edges), warpEdges(spec.stop.edges)
			if spec.pass.btype == bandstopType {
				p = symmetricBandstop(p, s)
			}
			lowpass := lowpassSpec{
				nat:   stopbandRatio(spec.pass.btype, p, s),
				gpass: math.Expm1(spec.rippleDB * math.Ln10 / 10),
				gstop: math.Expm1(spec.attenuation * math.Ln10 / 10),
			}
			rho := lowpass.prototypeWn(sel.design, MatchPassband, order-1)
			edges := bandEdges(spec.pass.btype, p, rho)
			for i, x := range edges {
				edges[i] = 2 * math.Atan(x) / math.Pi
			}

			f, err := sel.filter(order-1, spec.rippleDB, spec.attenuation, Band{btype: spec.pass.btype, edges: edges})
			require.NoError(t, err, name)
			_, stopHighest := edgeGains(f, spec.stop.edges)
			assert.Greater(t, stopHighest, -spec.attenuation, name)
		}
	}
}

// Without an option, the passband edge is met, as Matching(MatchPassband).
func TestOrderMatchesPassbandByDefault(t *testing.T) {
	for _, sel := range orderSelections {
		for _, spec := range orderSpecs {
			order, band, err := sel.selectOrder(spec.pass, spec.stop, spec.rippleDB, spec.attenuation)
			require.NoError(t, err)
			wantOrder, wantBand, err := sel.selectOrder(spec.pass, spec.stop, spec.rippleDB, spec.attenuation, Matching(MatchPassband))
			require.NoError(t, err)
			assert.Equal(t, wantOrder, order, "%s, %s", sel.name, spec.name)
			assert.Equal(t, wantBand, band, "%s, %s", sel.name, spec.name)
		}
	}
}

// Where a design's band edge is the matched edge itself, the selected band is
// that edge: Cheby1 and Ellip's passband edges, and Cheby2's stopband edges.
// A bandpass's looser stopband edge and a bandstop's moved passband edge are
// the exceptions, so only lowpass and highpass are compared.
func TestOrderBandIsTheMatchedEdge(t *testing.T) {
	cases := []struct {
		sel   int
		match MatchBand
		stop  bool
	}{
		{1, MatchPassband, false}, // Cheby1Ord
		{3, MatchPassband, false}, // EllipOrd
		{2, MatchStopband, true},  // Cheby2Ord
	}
	for _, c := range cases {
		sel := orderSelections[c.sel]
		for _, spec := range orderSpecs[:4] {
			_, band, err := sel.selectOrder(spec.pass, spec.stop, spec.rippleDB, spec.attenuation, Matching(c.match))
			require.NoError(t, err)
			want := spec.pass.edges
			if c.stop {
				want = spec.stop.edges
			}
			assert.InEpsilonSlice(t, want, band.Edges(), 1e-15, "%s, %s", sel.name, spec.name)
		}
	}
}

// symmetricBandstop moves one passband edge inward until p1·p2 = s1·s2,
// leaving the other where it was.
func TestSymmetricBandstop(t *testing.T) {
	cases := []struct{ p, s []float64 }{
		{[]float64{0.1, 5}, []float64{0.5, 2}},    // p1·p2 < s1·s2: p1 moves up
		{[]float64{0.1, 20}, []float64{0.5, 2}},   // p1·p2 > s1·s2: p2 moves down
		{[]float64{0.25, 4}, []float64{0.5, 2}},   // already symmetric
		{[]float64{1e-3, 1e3}, []float64{1, 1.1}}, // far apart
	}
	for _, c := range cases {
		got := symmetricBandstop(c.p, c.s)
		assert.InEpsilon(t, c.s[0]*c.s[1], got[0]*got[1], 1e-15, "%v, %v", c.p, c.s)
		assert.True(t, got[0] == c.p[0] || got[1] == c.p[1], "one edge stays: %v from %v", got, c.p)
		assert.GreaterOrEqual(t, got[0], c.p[0], "inward: %v from %v", got, c.p)
		assert.LessOrEqual(t, got[1], c.p[1], "inward: %v from %v", got, c.p)
	}
}

// bandEdges inverts the lowpass-to-band transform: the equivalent lowpass is
// at rho at both edges, and they multiply to p1·p2, even for a rho so large
// or small that the lower edge would cancel if found by a difference.
func TestBandEdgesInvertTheTransform(t *testing.T) {
	p := []float64{0.5, 2}
	width, centreSq := p[1]-p[0], p[0]*p[1]
	for _, btype := range []bandType{bandpassType, bandstopType} {
		for _, rho := range []float64{1e-6, 0.5, 1, 3, 1e6} {
			edges := bandEdges(btype, p, rho)
			assert.InEpsilon(t, centreSq, edges[0]*edges[1], 1e-15, "%d, rho %g", btype, rho)
			for _, x := range edges {
				lowpass := math.Abs((x*x - centreSq) / (x * width))
				if btype == bandstopType {
					lowpass = 1 / lowpass
				}
				assert.InEpsilon(t, rho, lowpass, 1e-9, "%d, rho %g, edge %g", btype, rho, x)
			}
		}
	}
}

// The elliptic raw order keeps its digits as nat nears 1: against the same
// formula with 1 − 1/nat² found exactly, it agrees to 1e-14. 1 − 1/nat²
// computed directly loses up to 9e-9 of itself near nat − 1 = 1e-9, about
// 2e-10 of the order, which this would catch.
func TestEllipRawOrderNearOne(t *testing.T) {
	gpass, gstop := math.Expm1(0.1*math.Ln10/10), math.Expm1(40*math.Ln10/10)
	m1 := gpass / gstop
	for i := range 400 {
		nat := 1 + math.Pow(10, -12+8*float64(i)/399) // nat − 1 from 1e-12 to 1e-4

		// 1 − 1/nat² = (nat² − 1)/nat², exactly from nat as passed, then rounded
		exactNat := new(big.Float).SetPrec(200).SetFloat64(nat)
		natSq := new(big.Float).SetPrec(200).Mul(exactNat, exactNat)
		excess := new(big.Float).SetPrec(200).Sub(natSq, big.NewFloat(1))
		mc, _ := excess.Quo(excess, natSq).Float64()

		want := special.EllipKM1(mc) * special.EllipKM1(m1) / (special.EllipKM1(1-mc) * special.EllipKM1(1-m1))
		got := lowpassSpec{nat: nat, gpass: gpass, gstop: gstop}.rawOrder(orderEllip)
		assert.InEpsilon(t, want, got, 1e-14, "nat − 1 = %g", nat-1)
	}
}

// Edges returns a copy, so changing it leaves the band as it was.
func TestBandEdgesIsACopy(t *testing.T) {
	band := Bandpass(0.2, 0.5)
	edges := band.Edges()
	assert.Equal(t, []float64{0.2, 0.5}, edges)
	edges[0] = 0.9
	assert.Equal(t, []float64{0.2, 0.5}, band.Edges())
}

func TestOrderRejectsInvalidInput(t *testing.T) {
	for _, sel := range orderSelections {
		check := func(pass, stop Band, rippleDB, attenuationDB float64, opts []OrderOption, want error, what string) {
			t.Helper()
			_, _, err := sel.selectOrder(pass, stop, rippleDB, attenuationDB, opts...)
			assert.ErrorIs(t, err, want, "%s: %s", sel.name, what)
		}
		lowPass, lowStop := Lowpass(0.2), Lowpass(0.3)

		for _, dB := range []float64{0, -1, nan, math.Inf(1), math.Inf(-1)} {
			check(lowPass, lowStop, dB, 40, nil, ErrInvalidParameter, fmt.Sprintf("rippleDB %g", dB))
			check(lowPass, lowStop, 1, dB, nil, ErrInvalidParameter, fmt.Sprintf("attenuationDB %g", dB))
		}
		check(lowPass, lowStop, 1, 1, nil, ErrInvalidParameter, "attenuationDB equal to rippleDB")
		check(lowPass, lowStop, 3, 1, nil, ErrInvalidParameter, "attenuationDB below rippleDB")

		bands := []struct {
			pass, stop Band
			what       string
		}{
			{Band{}, lowStop, "zero passband"},
			{lowPass, Band{}, "zero stopband"},
			{Lowpass(1.2), lowStop, "passband outside (0, 1)"},
			{lowPass, Lowpass(0), "stopband outside (0, 1)"},
			{lowPass, Highpass(0.3), "different band types"},
			{lowPass, Lowpass(0.1), "lowpass stopband below passband"},
			{lowPass, Lowpass(0.2), "lowpass stopband at passband"},
			{Highpass(0.2), Highpass(0.3), "highpass stopband above passband"},
			{Bandpass(0.2, 0.5), Bandpass(0.3, 0.6), "bandpass stopband inside passband"},
			{Bandpass(0.2, 0.5), Bandpass(0.1, 0.4), "bandpass stopband inside passband, upper"},
			{Bandstop(0.2, 0.5), Bandstop(0.1, 0.4), "bandstop stopband outside passband"},
			{Bandstop(0.2, 0.5), Bandstop(0.3, 0.6), "bandstop stopband outside passband, upper"},
		}
		for _, b := range bands {
			check(b.pass, b.stop, 1, 40, nil, ErrInvalidBand, b.what)
		}

		for _, match := range []MatchBand{-1, 2} {
			check(lowPass, lowStop, 1, 40, []OrderOption{Matching(match)}, ErrInvalidMode, fmt.Sprintf("MatchBand %d", match))
		}
	}
}

// Edges one ULP apart put nat − 1 near 1e-16. A Butterworth order grows as
// 1/(nat − 1), past what an int holds; the others grow as 1/√(nat − 1) or
// ln(1/(nat − 1)), and their true orders are returned.
func TestOrderEdgesTooClose(t *testing.T) {
	pass, stop := Lowpass(0.3), Lowpass(math.Nextafter(0.3, 1))
	_, _, err := ButterOrd(pass, stop, 1, 40)
	assert.ErrorIs(t, err, ErrInvalidBand)

	for _, sel := range orderSelections[1:] {
		order, _, err := sel.selectOrder(pass, stop, 1, 40)
		require.NoError(t, err, sel.name)
		assert.Greater(t, order, 1, sel.name)
	}
}
