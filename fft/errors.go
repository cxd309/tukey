package fft

import "errors"

// Errors returned for invalid input.
// Each is wrapped with the detail of what was wrong, so check for one with errors.Is.
// They mean the same in every tukey package.
var (
	// ErrInvalidLength is a length that's wrong in itself.
	ErrInvalidLength = errors.New("fft: invalid length")
	// ErrSignalTooShort is an input with too few samples, including none.
	ErrSignalTooShort = errors.New("fft: signal too short")
	// ErrInvalidParameter is a value outside its domain.
	ErrInvalidParameter = errors.New("fft: invalid parameter")
	// ErrInvalidOption is an option that's unknown or doesn't apply.
	ErrInvalidOption = errors.New("fft: invalid option")
)
