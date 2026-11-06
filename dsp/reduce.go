package dsp

import (
	"fmt"
	"math"

	"github.com/cxd309/tukey/internal/floats"
)

// RMS returns the root mean square of x, sqrt(mean(x^2))
//
// returns ErrInvalidLength if x is empty (MATLAB's rms gives NaN)
//
// equivalent to MATLAB's rms(x).
// equivalent to numpy.sqrt(numpy.mean(x**2)), no SciPy equivalent
func RMS(x []float64) (rms float64, err error) {
	if len(x) == 0 {
		return 0, fmt.Errorf("%w: RMS of an empty signal", ErrInvalidLength)
	}
	var sum floats.Sum
	for _, v := range x {
		sum.Add(v * v)
	}
	rms = math.Sqrt(sum.Total() / float64(len(x)))
	return
}

// mean is Sum(x) divided by its length, for internal use where x is known not to be empty.
// NaN for an empty x
func mean(x []float64) (m float64) {
	m = Sum(x) / float64(len(x))
	return
}

// Sum returns the total of x, 0 for an empty x.
// It's summed with compensation, to maintain low bit accuracy.
//
// Equivalent to numpy.sum(x), which sums pairwise, and MATLAB sum(x), which sums plainly.
// On such signals Sum is closer to the exact total than either:
// one value of 1e8 among 20,000 of 0.001, then cancelled, sum to 20, where MATLAB gives 20.00001
func Sum(x []float64) (sum float64) {
	var s floats.Sum
	for _, v := range x {
		s.Add(v)
	}
	sum = s.Total()
	return
}

// Mean returns the arithmetic mean of x: its Sum divided by its length.
//
// Returns ErrInvalidLength for an empty x, which has no mean (MATLAB's mean gives NaN)
//
// Equivalent to numpy.mean(x) and MATLAB mean(x), summed with compensation like Sum
func Mean(x []float64) (m float64, err error) {
	if len(x) == 0 {
		return 0, fmt.Errorf("%w: Mean of an empty signal", ErrInvalidLength)
	}
	m = mean(x)
	return
}
