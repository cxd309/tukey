package dsp

import (
	"math"
	"slices"
)

// Fill is how FillMissing replaces missing (NaN) samples:
// FillLinear (the zero Fill) or FillConstant
type Fill struct {
	constant bool
	value    float64
}

// FillLinear fills each gap by linear interpolation between the known samples
// either side, and each end by extending the line through its nearest two
func FillLinear() (fill Fill) {
	return
}

// FillConstant fills every missing sample with value
func FillConstant(value float64) (fill Fill) {
	fill = Fill{constant: true, value: value}
	return
}

// FillMissing returns x with each NaN replaced as fill says
// with FillLinear, fewer than two known samples leaves x unchanged,
// as there's no line to draw
//
// equivalent to MATLAB's fillmissing(x, "linear") and
// fillmissing(x, "constant", value).
// No SciPy equivalent
func FillMissing(x []float64, fill Fill) (y []float64) {
	if fill.constant {
		y = fillConstant(x, fill.value)
		return
	}
	y = fillLinear(x)
	return
}

// fillConstant replaces each NaN in x with value
func fillConstant(x []float64, value float64) (y []float64) {
	y = slices.Clone(x)
	for i, v := range y {
		if math.IsNaN(v) {
			y[i] = value
		}
	}
	return
}

// fillLinear replaces each NaN in x by the line through the known samples
// either side of it, or, before the first or after the last known sample,
// through the nearest two
func fillLinear(x []float64) (y []float64) {
	y = slices.Clone(x)
	known := knownIndices(x)
	if len(known) < 2 {
		return
	}
	k := 0 // known[k] and known[k+1] bracket the sample being filled
	for i, v := range x {
		if !math.IsNaN(v) {
			continue
		}
		for k+2 < len(known) && known[k+1] < i {
			k++
		}
		lo, hi := known[k], known[k+1]
		y[i] = line(float64(lo), x[lo], float64(hi), x[hi], float64(i))
	}
	return
}

// knownIndices is the index of every sample of x that isn't NaN
func knownIndices(x []float64) (known []int) {
	for i, v := range x {
		if !math.IsNaN(v) {
			known = append(known, i)
		}
	}
	return
}
