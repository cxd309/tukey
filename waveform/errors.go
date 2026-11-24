package waveform

import "errors"

// Errors returned for invalid input.
// Each is wrapped with the detail of what was wrong, so check for one with errors.Is.
// They mean the same in every tukey package.
var (
	// ErrInvalidParameter is a value outside its domain.
	ErrInvalidParameter = errors.New("waveform: invalid parameter")
	// ErrInvalidOption is an option that's unknown or doesn't apply.
	ErrInvalidOption = errors.New("waveform: invalid option")
)
