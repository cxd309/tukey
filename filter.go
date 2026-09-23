package dsp

import "fmt"

// BandType selects which frequencies a filter passes
type BandType int

const (
	LowPass  BandType = iota // passes below the cutoff
	HighPass                 // passes above the cutoff
	BandPass                 // passes between the two band edges
	BandStop                 // rejects between the two band edges
)

// String returns the band type's name, matching SciPy's btype strings
func (bt BandType) String() (name string) {
	switch bt {
	case LowPass:
		name = "lowpass"
	case HighPass:
		name = "highpass"
	case BandPass:
		name = "bandpass"
	case BandStop:
		name = "bandstop"
	default:
		name = fmt.Sprintf("BandType(%d)", int(bt))
	}
	return
}
