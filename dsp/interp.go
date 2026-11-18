package dsp

import (
	"fmt"
	"math"
	"sort"

	"github.com/cxd309/tukey/internal/floats"
)

// Interp1 returns y, sampled at the strictly increasing points x, linearly interpolated at each point of xi.
// A point of xi outside [x[0], x[len(x)-1]], or a NaN one, gives NaN.
// Next to an infinite y the line is that infinity, its limit.
//
//	yi, err := dsp.Interp1([]float64{0, 1, 2}, []float64{0, 10, 0}, []float64{0.5, 1.5}) // [5 5]
//
// Units:
//   - x and xi share a unit
//   - yi is in y's unit
//
// Returns:
//   - ErrInvalidLength for an x and y of different lengths, or shorter than 2.
//   - ErrInvalidSamplePoints for an x that isn't strictly increasing, including one with a NaN.
//
// Equivalent to:
//   - scipy.interpolate.interp1d(x, y, bounds_error=False)(xi), where numpy.interp repeats the end values outside x's range instead
//   - MATLAB interp1(x, y, xi)
func Interp1(x, y, xi []float64) (yi []float64, err error) {
	if len(x) != len(y) || len(x) < 2 {
		return nil, fmt.Errorf("%w: x and y need the same length, at least 2, got %d and %d",
			ErrInvalidLength, len(x), len(y))
	}
	for i := 1; i < len(x); i++ {
		if !(x[i] > x[i-1]) { // false for NaN too
			return nil, fmt.Errorf("%w: x must be strictly increasing, x[%d] = %g and x[%d] = %g",
				ErrInvalidSamplePoints, i-1, x[i-1], i, x[i])
		}
	}
	yi = make([]float64, len(xi))
	for j, q := range xi {
		yi[j] = interpolate(x, y, q)
	}
	return
}

// FillMissing returns x with each NaN replaced as fill says.
// With FillLinear, fewer than two known samples leave x unchanged, as there's no line to draw.
//
//	y := dsp.FillMissing([]float64{1, math.NaN(), 3}, dsp.FillLinear()) // [1 2 3]
//
// Units: same as x
//
// Equivalent to:
//   - SciPy has none
//   - MATLAB fillmissing(x, 'linear'), or fillmissing(x, 'constant', value)
func FillMissing(x []float64, fill Fill) (y []float64) {
	if fill.constant {
		y = fillConstant(x, fill.value)
		return
	}
	y = fillLinear(x)
	return
}

// Fill is how FillMissing replaces missing samples, which are NaN:
// FillLinear, the zero Fill, or FillConstant.
type Fill struct {
	constant bool
	value    float64
}

// FillLinear fills each gap by linear interpolation between the known samples either side,
// and each end by extending the line through its nearest two.
// It's the zero Fill.
//
// Equivalent to:
//   - SciPy has none
//   - MATLAB fillmissing(x, 'linear')
func FillLinear() (fill Fill) {
	return
}

// FillConstant fills every missing sample with value.
// FillConstant(math.NaN()) leaves x as it is.
//
// Units: value is in x's unit
//
// Equivalent to:
//   - SciPy has none
//   - MATLAB fillmissing(x, 'constant', value)
func FillConstant(value float64) (fill Fill) {
	fill = Fill{constant: true, value: value}
	return
}

// interpolate is y at q:
// y itself on a sample point, the line through the samples either side between them,
// and NaN outside x's range or for a NaN q.
func interpolate(x, y []float64, q float64) (v float64) {
	i := sort.SearchFloat64s(x, q) // the first x[i] >= q, len(x) for NaN
	switch {
	case i < len(x) && x[i] == q:
		return y[i]
	case i == 0 || i == len(x):
		return math.NaN()
	}
	v = line(x[i-1], y[i-1], x[i], y[i], q)
	return
}

// line is the value at q of the straight line through (x0, y0) and (x1, y1),
// or the limit next to an infinite y.
func line(x0, y0, x1, y1, q float64) (v float64) {
	if v, ok := floats.InfiniteEnd(y0, y1, (q-x0)/(x1-x0)); ok {
		return v
	}
	slope := (y1 - y0) / (x1 - x0)
	v = slope*(q-x0) + y0
	return
}

// fillConstant replaces each NaN in x with value.
func fillConstant(x []float64, value float64) (y []float64) {
	y = make([]float64, len(x))
	copy(y, x)
	for i, v := range y {
		if math.IsNaN(v) {
			y[i] = value
		}
	}
	return
}

// fillLinear replaces each NaN in x by the line through the known samples either side of it,
// or, before the first or after the last known sample, through the nearest two.
func fillLinear(x []float64) (y []float64) {
	y = make([]float64, len(x))
	copy(y, x)
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

// knownIndices is the index of every sample of x that isn't NaN.
func knownIndices(x []float64) (known []int) {
	for i, v := range x {
		if !math.IsNaN(v) {
			known = append(known, i)
		}
	}
	return
}
