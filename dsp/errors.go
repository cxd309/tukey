package dsp

import "errors"

// Errors returned for invalid input
// each is wrapped with the specific detail
// so can be checked with errors.Is(err, dsp.ErrSignalTooShort)
var (
	ErrInvalidOrder        = errors.New("dsp: invalid filter order")
	ErrInvalidBand         = errors.New("dsp: invalid band")
	ErrInvalidCoefficients = errors.New("dsp: invalid filter coefficients")
	ErrInvalidPadding      = errors.New("dsp: invalid padding")
	ErrSignalTooShort      = errors.New("dsp: signal too short")
	ErrNoSteadyState       = errors.New("dsp: filter has no steady state")
)
