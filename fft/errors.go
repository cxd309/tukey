package fft

import "errors"

// Errors returned for invalid input
// each is wrapped with the specifics of what was wrong
// check with errors.Is
var (
	ErrInvalidLength = errors.New("fft: invalid length")
	ErrInvalidOption = errors.New("fft: invalid option")
)
