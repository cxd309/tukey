package dsp

import (
	"fmt"
	"math"
)

// Band is the frequency response a filter design targets: which frequencies
// it passes and where its edges are
// Build one with Lowpass, Highpass, Bandpass or Bandstop.
//
// Frequencies are normalised to Nyquist, so they lie in (0, 1);
// for a frequency f in Hz at sample rate fs, use f / (fs/2)
type Band struct {
	btype bandType  // SciPy's name for the same parameter
	edges []float64 // the cutoff for lowpass/highpass; low, high for bandpass/bandstop
}

// bandType starts at 1 so the zero Band is detectably invalid
type bandType int

const (
	lowpassType bandType = iota + 1
	highpassType
	bandpassType
	bandstopType
)

// Lowpass passes frequencies below cutoff
func Lowpass(cutoff float64) (band Band) {
	band = Band{btype: lowpassType, edges: []float64{cutoff}}
	return
}

// Highpass passes frequencies above cutoff
func Highpass(cutoff float64) (band Band) {
	band = Band{btype: highpassType, edges: []float64{cutoff}}
	return
}

// Bandpass passes frequencies between low and high
func Bandpass(low, high float64) (band Band) {
	band = Band{btype: bandpassType, edges: []float64{low, high}}
	return
}

// Bandstop rejects frequencies between low and high
func Bandstop(low, high float64) (band Band) {
	band = Band{btype: bandstopType, edges: []float64{low, high}}
	return
}

// String formats the band the way it was constructed, e.g. "Bandpass(0.2, 0.5)"
func (band Band) String() (s string) {
	switch band.btype {
	case lowpassType:
		s = fmt.Sprintf("Lowpass(%v)", band.edges[0])
	case highpassType:
		s = fmt.Sprintf("Highpass(%v)", band.edges[0])
	case bandpassType:
		s = fmt.Sprintf("Bandpass(%v, %v)", band.edges[0], band.edges[1])
	case bandstopType:
		s = fmt.Sprintf("Bandstop(%v, %v)", band.edges[0], band.edges[1])
	default:
		s = "Band{}"
	}
	return
}

// validate reports whether the band can be designed:
// 1. it must come from a constructor
// 2. its edges must be ascending within (0, 1)
func (band Band) validate() (err error) {
	if band.btype == 0 {
		return fmt.Errorf("%w: zero Band; build one with Lowpass, Highpass, Bandpass or Bandstop", ErrInvalidBand)
	}
	for i, edge := range band.edges {
		// written as !(in range) rather than (out of range) so NaN is rejected
		if !(edge > 0 && edge < 1) {
			return fmt.Errorf("%w: %v: frequencies must be in (0, 1), normalised to Nyquist", ErrInvalidBand, band)
		}
		if i > 0 && edge <= band.edges[i-1] {
			return fmt.Errorf("%w: %v: low edge must be below high edge", ErrInvalidBand, band)
		}
	}
	return
}

// transform prewarps the band edges and moves the lowpass prototype onto them,
// giving the analog filter that bilinear will map onto this band
// the band must already have passed validate
func (band Band) transform(proto zpk) (analog zpk) {
	switch band.btype {
	case lowpassType:
		analog = proto.toLowpass(prewarp(band.edges[0]))
	case highpassType:
		analog = proto.toHighpass(prewarp(band.edges[0]))
	case bandpassType, bandstopType:
		low, high := prewarp(band.edges[0]), prewarp(band.edges[1])
		wo := math.Sqrt(low * high) // geometric centre
		bw := high - low
		if band.btype == bandpassType {
			analog = proto.toBandpass(wo, bw)
		} else {
			analog = proto.toBandstop(wo, bw)
		}
	default:
		panic(fmt.Sprintf("dsp: transform called on unvalidated band %v", band))
	}
	return
}
