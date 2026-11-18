package dsp

import (
	"github.com/cxd309/tukey/internal/floats"
)

// Detrend returns x with its trend removed:
// the least-squares line for DetrendLinear, or the mean for DetrendConstant.
// An empty x gives an empty result, and a single sample gives [0].
//
//	y := dsp.Detrend(x, dsp.DetrendLinear())
//
// Units: same as x
//
// Equivalent to:
//   - scipy.signal.detrend(x, type='linear'), or type='constant'
//   - MATLAB detrend(x), or detrend(x, 0)
func Detrend(x []float64, mode DetrendMode) (y []float64) {
	if mode.constant {
		y = detrendConstant(x)
	} else {
		y = detrendLinear(x)
	}
	return
}

// DetrendMode is what Detrend removes, built by DetrendLinear or DetrendConstant.
// Its zero value is DetrendLinear.
type DetrendMode struct {
	constant bool
}

// DetrendLinear removes the least-squares straight line.
// It's DetrendMode's zero value.
//
// Equivalent to:
//   - scipy.signal.detrend's type='linear', its default
//   - MATLAB detrend(x), or detrend(x, 1)
func DetrendLinear() (mode DetrendMode) {
	return
}

// DetrendConstant removes the mean.
//
// Equivalent to:
//   - scipy.signal.detrend's type='constant'
//   - MATLAB detrend(x, 0)
func DetrendConstant() (mode DetrendMode) {
	mode.constant = true
	return
}

// detrendConstant subtracts the mean.
func detrendConstant(x []float64) (y []float64) {
	m := mean(x)
	y = make([]float64, len(x))
	for i, v := range x {
		y[i] = v - m
	}
	return
}

// detrendLinear subtracts the least-squares line, fitted about the centre (n-1)/2,
// so the intercept is just the mean and the sums stay small.
func detrendLinear(x []float64) (y []float64) {
	centre := float64(len(x)-1) / 2
	m := mean(x)
	var sxy, sxx floats.Sum
	for i, v := range x {
		t := float64(i) - centre
		sxy.Add(t * (v - m))
		sxx.Add(t * t)
	}
	slope := 0.0
	if sxx.Total() > 0 { // a single sample has no slope
		slope = sxy.Total() / sxx.Total()
	}
	y = make([]float64, len(x))
	for i, v := range x {
		y[i] = v - m - slope*(float64(i)-centre)
	}
	return
}
