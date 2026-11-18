package poly

import "errors"

// ErrUnpairedConjugate is returned when a complex root has no matching conjugate,
// so a polynomial built from the roots can't have real coefficients.
// Callers in dsp wrap it with their own sentinel, as poly is internal.
var ErrUnpairedConjugate = errors.New("complex roots must come in conjugate pairs")
