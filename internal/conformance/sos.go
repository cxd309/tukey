package conformance

import (
	"math"
	"slices"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

type zpk2sosVector struct {
	reference.Meta
	Params struct {
		Z reference.Complex `json:"z"`
		P reference.Complex `json:"p"`
		K float64           `json:"k"`
	} `json:"params"`
	Output struct {
		Sos [][6]float64 `json:"sos"`
	} `json:"output"`
}

// zpk2sosSuite checks ZPK.SOS against scipy.signal.zpk2sos
// on SciPy's own zeros, poles and gain,
// isolating the pairing from any difference in designs
// measured: section order and pairing match exactly on all 29 vectors, worst
// coefficient error 1.1e-16 (under 1 ULP); set to 1e-14 to match Filter
var zpk2sosSuite = newSuite("ZPK.SOS", "zpk2sos",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	func(v zpk2sosVector) (outputs []Output, err error) {
		zeros, err := v.Params.Z.Values()
		if err != nil {
			return nil, err
		}
		poles, err := v.Params.P.Values()
		if err != nil {
			return nil, err
		}
		sos, err := dsp.ZPK{Zeros: zeros, Poles: poles, Gain: v.Params.K}.SOS()
		if err != nil {
			return nil, err
		}
		outputs, err = sosOutputs(sos, v.Output.Sos)
		return
	})

// butterSOSSuite checks dsp.Butter(...).SOS() against
// scipy.signal.butter(..., output='sos') end to end,
// using the sos output recorded in the Butter vectors
//
// sections are compared as a set, with the overall gain taken out and compared
// separately: when poles tie in distance from the unit circle, rounding decides
// both the section order and which section carries the gain, and either way is
// the same filter. bandstop_order2_wn0.01-0.99 is such a case: its edges are
// symmetric about half Nyquist, so its pole pairs are p and -p. Order and gain
// placement are checked exactly by zpk2sosSuite, on identical inputs.
//
// measured: worst 7.3e-15 (~33 ULP), on the same vector as Butter's own worst
// case; the conversion adds essentially nothing (see zpk2sosSuite), so SOS
// inherits the design's pole accuracy; set to match Butter
var butterSOSSuite = newSuite("Butter (SOS)", "butter",
	reference.Tolerance{Rel: 1e-13, Scaled: 1e-13},
	func(v butterVector) (outputs []Output, err error) {
		wn, err := parseWn(v.Params.Wn)
		if err != nil {
			return nil, err
		}
		band, err := bandFromVector(v.Params.BType, wn)
		if err != nil {
			return nil, err
		}
		f, err := dsp.Butter(v.Params.Order, band)
		if err != nil {
			return nil, err
		}
		sos, err := f.SOS()
		if err != nil {
			return nil, err
		}
		gotSections, gotGain := normalizedSections(sos)
		wantSections, wantGain := normalizedSections(v.Output.Sos)
		outputs, err = sosOutputs(matchSections(gotSections, wantSections), wantSections)
		if err != nil {
			return nil, err
		}
		outputs = append(outputs, Output{Name: "gain", Got: []float64{gotGain}, Want: []float64{wantGain}})
		return
	})

// normalizedSections takes the gain out of each section's b, dividing it by its
// leading non-zero coefficient; the product of those is the overall gain, so
// sections can be compared regardless of which one carries it
func normalizedSections(sos [][6]float64) (sections [][6]float64, gain float64) {
	gain = 1
	sections = slices.Clone(sos) // arrays are values, so this copies every coefficient
	for i := range sections {
		lead := leadingCoefficient(sections[i][:3])
		gain *= lead
		for j := range 3 {
			sections[i][j] /= lead
		}
	}
	return
}

// leadingCoefficient is b's first non-zero coefficient, or 1 if there is none
func leadingCoefficient(b []float64) (lead float64) {
	for _, c := range b {
		if c != 0 {
			return c
		}
	}
	return 1
}

// matchSections reorders got so each section lines up with the closest section
// in want; see butterSOSSuite for why order can legitimately differ
func matchSections(got dsp.SOS, want [][6]float64) (matched dsp.SOS) {
	used := make([]bool, len(got))
	for _, w := range want {
		best, bestDist := -1, math.Inf(1)
		for i, g := range got {
			if used[i] {
				continue
			}
			if d := sectionDistance(g, w); d < bestDist {
				best, bestDist = i, d
			}
		}
		if best < 0 {
			break // more sections wanted than got; sosOutputs reports the count mismatch
		}
		used[best] = true
		matched = append(matched, got[best])
	}
	return
}

// sectionDistance is the largest coefficient difference between two sections
func sectionDistance(g, w [6]float64) (d float64) {
	for i := range g {
		d = math.Max(d, math.Abs(g[i]-w[i]))
	}
	return
}
