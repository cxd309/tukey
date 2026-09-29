package dsp

import "fmt"

// Padding selects how FiltFilt extends the signal past each end before filtering
// so that each pass's startup transient falls on the extension rather than the data
type Padding int

const (
	// PaddingOdd reflect the signal through its end sample
	// so both value and slope continue smoothly
	// default option like SciPy and MATLAB
	PaddingOdd Padding = iota
	// PaddingEven mirror the signal about its end sample
	PaddingEven
	// PaddingConstant repeats the end sample
	PaddingConstant
	// PaddingNone filters the signal as given
	PaddingNone
)

func (p Padding) String() (name string) {
	switch p {
	case PaddingOdd:
		name = "PaddingOdd"
	case PaddingEven:
		name = "PaddingEven"
	case PaddingConstant:
		name = "PaddingConstant"
	case PaddingNone:
		name = "PaddingNone"
	default:
		name = fmt.Sprintf("Padding(%d)", int(p))
	}
	return
}

// extend returns a new slice holding x with n extension samples before and after it
// it panics unless 0 <= n < len(x), and n == 0 for PaddingNone
// FiltFilt validates caller input first, so a panic here means a bug inside dsp
// equvalent to scipy.signal's odd_ext, even_ext and const_ext
func (p Padding) extend(x []float64, n int) (extended []float64) {
	if !(n >= 0 && n < len(x)) {
		panic(fmt.Sprintf("dsp: extend needs 0 <= n < len(x), got n=%d, len(x)=%d", n, len(x)))
	}
	if p == PaddingNone && n != 0 {
		panic(fmt.Sprintf("dsp: extend with PaddingNone needs n = 0, got n=%d", n))
	}
	first, last := x[0], x[len(x)-1]
	extended = make([]float64, 0, len(x)+2*n)
	for i := n; i >= 1; i-- {
		extended = append(extended, p.reflect(first, x[i]))
	}
	extended = append(extended, x...)
	for i := 1; i <= n; i++ {
		extended = append(extended, p.reflect(last, x[len(x)-1-i]))
	}
	return
}

// reflect is the extensions sample that mirrors sample v about the edge sample edge
func (p Padding) reflect(edge, v float64) (sample float64) {
	switch p {
	case PaddingOdd:
		sample = 2*edge - v
	case PaddingEven:
		sample = v
	case PaddingConstant:
		sample = edge
	default:
		panic(fmt.Sprintf("dsp: unknown padding %v", p))
	}
	return
}
