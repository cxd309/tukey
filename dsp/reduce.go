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

// mean is the arithmetic mean of x. NaN for an empty x
func mean(x []float64) (m float64) {
	var sum floats.Sum
	for _, v := range x {
		sum.Add(v)
	}
	m = sum.Total() / float64(len(x))
	return
}
