package dsp

import (
	"fmt"
	"math"
	"sort"
)

// Interp1 returns y, sampled at strictly increasing points x, linearly,
// interpolated at each point of xi.
// A point outside [x[0], x[len(x)-1]], or NaN, gives NaN
//
// returns ErrInvalidLength unless x and y have the same length, of at least 2.
// returns ErrInvalidSamplePoints unless x is strictly increasing (so no NaN)
//
// equivalent to MATLAB interp1(x, y, xi).
// equivalent to scipy.interpolate.interp1d(x, y, bounds_error=False)(xi).
// numpy.interp differs outside x's range, where it repeats the end values
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

// interpolate is y at q: y itself on a sample point, the line through the
// samples either side between them, and NaN outside x's range or for NaN q
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

// line is the value at q of the straight line through (x0, y0) and (x1, y1)
func line(x0, y0, x1, y1, q float64) (v float64) {
	slope := (y1 - y0) / (x1 - x0)
	v = slope*(q-x0) + y0
	return
}
