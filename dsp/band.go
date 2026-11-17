package dsp

import (
	"fmt"
)

// Band is the frequency response a filter design targets:
// which frequencies it passes, and where its edges are.
// Build one with Lowpass, Highpass, Bandpass or Bandstop.
// Its zero value is no band, which every design rejects with ErrInvalidBand.
//
// Its edges are fractions of Nyquist, so they lie in (0, 1).
// HzToNyquist converts from Hz:
//
//	wn, err := dsp.HzToNyquist(100.0, fs)
//	f, err := dsp.Butter(4, dsp.Lowpass(wn))
type Band struct {
	btype bandType  // SciPy's name for the same parameter
	edges []float64 // the cutoff for a lowpass or highpass, and low, high for a bandpass or bandstop
}

// Lowpass is the band below cutoff.
// A design checks cutoff, so Lowpass itself can't fail.
//
// Units: cutoff is a fraction of Nyquist
//
// Equivalent to:
//   - SciPy's btype='lowpass', with Wn = cutoff
//   - MATLAB's default ftype, with Wn = cutoff
func Lowpass(cutoff float64) (band Band) {
	band = Band{btype: lowpassType, edges: []float64{cutoff}}
	return
}

// Highpass is the band above cutoff.
// A design checks cutoff, so Highpass itself can't fail.
//
// Units: cutoff is a fraction of Nyquist
//
// Equivalent to:
//   - SciPy's btype='highpass', with Wn = cutoff
//   - MATLAB's 'high', with Wn = cutoff
func Highpass(cutoff float64) (band Band) {
	band = Band{btype: highpassType, edges: []float64{cutoff}}
	return
}

// Bandpass is the band between low and high.
// A design checks the edges, so Bandpass itself can't fail.
//
// Units: low and high are fractions of Nyquist
//
// Equivalent to:
//   - SciPy's btype='bandpass', with Wn = [low, high]
//   - MATLAB's 'bandpass', with Wn = [low high]
func Bandpass(low, high float64) (band Band) {
	band = Band{btype: bandpassType, edges: []float64{low, high}}
	return
}

// Bandstop is everything but the band between low and high.
// A design checks the edges, so Bandstop itself can't fail.
//
// Units: low and high are fractions of Nyquist
//
// Equivalent to:
//   - SciPy's btype='bandstop', with Wn = [low, high]
//   - MATLAB's 'stop', with Wn = [low high]
func Bandstop(low, high float64) (band Band) {
	band = Band{btype: bandstopType, edges: []float64{low, high}}
	return
}

// String formats the band the way it was built, e.g. "Bandpass(0.2, 0.5)", or "Band{}" for the zero Band.
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

// Edges returns a copy of the band's edges:
// the cutoff for a lowpass or highpass, and low and high for a bandpass or bandstop.
// The zero Band has none.
//
// Units: the edges are fractions of Nyquist
func (band Band) Edges() (edges []float64) {
	edges = make([]float64, len(band.edges))
	copy(edges, band.edges)
	return
}

// bandType starts at 1, so the zero Band is detectably invalid.
type bandType int

const (
	lowpassType bandType = iota + 1
	highpassType
	bandpassType
	bandstopType
)

// validate returns ErrInvalidBand unless the band can be designed:
// it must come from a constructor, and its edges must ascend within (0, 1).
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
