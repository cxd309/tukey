package dsp

import "errors"

// Errors returned for invalid input.
// Each is wrapped with the detail of what was wrong,
// so check for one with errors.Is, e.g. errors.Is(err, dsp.ErrSignalTooShort).
// The first four mean the same in every tukey package.
var (
	// ErrInvalidLength is a length or count that's wrong in itself,
	// or lengths that should match and don't
	ErrInvalidLength = errors.New("dsp: invalid length")
	// ErrSignalTooShort is an input with too few samples for what's asked, including none
	ErrSignalTooShort = errors.New("dsp: signal too short")
	// ErrInvalidParameter is a value outside its domain
	ErrInvalidParameter = errors.New("dsp: invalid parameter")
	// ErrInvalidOption is an option or enum value that's unknown or doesn't apply
	ErrInvalidOption = errors.New("dsp: invalid option")

	// ErrInvalidBand is a Band that's zero, or whose edges aren't in (0, 1) and in order
	ErrInvalidBand = errors.New("dsp: invalid band")
	// ErrInvalidOrder is a filter order below 1
	ErrInvalidOrder = errors.New("dsp: invalid filter order")
	// ErrInvalidCoefficients is an empty b or a, or an a[0] of 0
	ErrInvalidCoefficients = errors.New("dsp: invalid filter coefficients")
	// ErrInvalidZPK is complex zeros or poles that aren't in conjugate pairs
	ErrInvalidZPK = errors.New("dsp: invalid zeros, poles and gain")
	// ErrInvalidSOS is no sections, or a section whose a0 isn't 1
	ErrInvalidSOS = errors.New("dsp: invalid second-order sections")
	// ErrInvalidState is a filter state of the wrong size
	ErrInvalidState = errors.New("dsp: invalid filter state")
	// ErrNoSteadyState is a filter with a pole at z = 1, which has no steady state
	ErrNoSteadyState = errors.New("dsp: filter has no steady state")
	// ErrInvalidSamplePoints is sample points that aren't strictly increasing
	ErrInvalidSamplePoints = errors.New("dsp: invalid sample points")
)
