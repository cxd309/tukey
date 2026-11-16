package window

import "errors"

// Errors returned for invalid input.
// Each is wrapped with the detail of what was wrong, so check for one with errors.Is.
// They mean the same in every tukey package.
var (
	// ErrInvalidLength is a window length that's wrong in itself
	ErrInvalidLength = errors.New("window: invalid length")
	// ErrInvalidParameter is a value outside its domain
	ErrInvalidParameter = errors.New("window: invalid parameter")
	// ErrInvalidOption is an option that doesn't apply to the window
	ErrInvalidOption = errors.New("window: invalid option")
)
