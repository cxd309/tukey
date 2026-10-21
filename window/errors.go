package window

import (
	"errors"
	"fmt"
)

// Errors returned for invalid input
// each is wrapped with the specifics of what was wrong
// check with errors.Is
var (
	ErrInvalidLength    = errors.New("window: invalid length")
	ErrInvalidParameter = errors.New("window: invalid parameter")
)

// positive returns ErrInvalidParameter unless v > 0, which also rules our NaN
func positive(window, name string, v float64) (err error) {
	if !(v > 0) {
		err = fmt.Errorf("%w: %s's %s must be >0, got %g", ErrInvalidParameter, window, name, v)
	}
	return
}
