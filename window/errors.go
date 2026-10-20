package window

import "errors"

// Errors returned for invalid input
// each is wrapped with the specifics of what was wrong
// check with errors.Is
var (
	ErrInvalidLength    = errors.New("window: invalid length")
	ErrInvalidParameter = errors.New("window: invalid parameter")
)
