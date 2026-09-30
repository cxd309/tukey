package dsp

import "github.com/cxd309/tukey/internal/floats"

// cascade is a filter as sections applied in turn, each section's output feeding the next
// a single (b, a) filter is a cascade of one, and an SOS is a cascade of second-order sections
type cascade []digitalFilter

// apply runs x through each section in turn, starting each from its state in initial
// every section at rest when initial is nil
// returns the output and each section's final state
// equivalent to scipy.signal.sosfilter, or lfilter for a single section
func (c cascade) apply(x []float64, initial [][]float64) (y []float64, final [][]float64) {
	y = x
	final = make([][]float64, len(c))
	for i, f := range c {
		var zi []float64
		if initial != nil {
			zi = initial[i]
		}
		y, final[i] = f.apply(y, zi)
	}
	return
}

// steadyState returns each section's state after a long unit step into the cascade
// each section's own steady state, scaled by the DC gain of the sections before it
// equivalent to scipy.signal.sosfilt_zi, or lfilter_zi for a single section
func (c cascade) steadyState() (zi [][]float64, err error) {
	zi = make([][]float64, len(c))
	scale := 1.0
	for i, f := range c {
		sectionZi, err := f.steadyState()
		if err != nil {
			return nil, err
		}
		zi[i] = floats.Scaled(sectionZi, scale)
		sumB, sumA := f.sums()
		scale *= sumB / sumA
	}
	return
}

// scaledState returns a copy of every section's state multiplied by s
func scaledState(zi [][]float64, s float64) (scaled [][]float64) {
	scaled = make([][]float64, len(zi))
	for i := range zi {
		scaled[i] = floats.Scaled(zi[i], s)
	}
	return
}
