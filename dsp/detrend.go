package dsp

import "fmt"

// DetrendMode selects what Detrend removes
type DetrendMode int

const (
	// DetrendLinear removes the least-squares straight line,
	// the default, same as SciPy and MATLAB
	DetrendLinear DetrendMode = iota
	// DetrendConstant removes the mean
	DetrendConstant
)

func (m DetrendMode) String() (name string) {
	switch m {
	case DetrendLinear:
		name = "DetrendLinear"
	case DetrendConstant:
		name = "DetrendConstant"
	default:
		name = fmt.Sprintf("DetrendMode(%d)", int(m))
	}
	return
}

// Detrend returns x with its trend removed:
// the least-squares line for DetrendLinear, the mean for DetrendConstant,
// empty x returns empty, and a single sample returns [0]
//
// returns ErrInvalidMode for an unknown mode
//
// equivalent to scipy.signal.detrend(x, type) and MATLAB's detrend(x, type)
func Detrend(x []float64, mode DetrendMode) (y []float64, err error) {
	switch mode {
	case DetrendLinear:
		y = detrendLinear(x)
	case DetrendConstant:
		y = detrendConstant(x)
	default:
		err = fmt.Errorf("%w: unknown detrend mode %v", ErrInvalidMode, mode)
	}
	return
}

// detrendConstant subtracts the mean
func detrendConstant(x []float64) (y []float64) {
	m := mean(x)
	y = make([]float64, len(x))
	for i, v := range x {
		y[i] = v - m
	}
	return
}

// detrendLinear subtracts the least-squares line, fitted about the centre
// (n-1)/2 so the intercept is just the mean and the sums stay small
func detrendLinear(x []float64) (y []float64) {
	centre := float64(len(x)-1) / 2
	m := mean(x)
	var sxy, sxx float64
	for i, v := range x {
		t := float64(i) - centre
		sxy += t * (v - m)
		sxx += t * t
	}
	slope := 0.0
	if sxx > 0 { // a single sample has no slope
		slope = sxy / sxx
	}
	y = make([]float64, len(x))
	for i, v := range x {
		y[i] = v - m - slope*(float64(i)-centre)
	}
	return
}
