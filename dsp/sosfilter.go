package dsp

import "fmt"

// toCascade validates the sections as SciPy's sosfilt does
// at least on section, and a0 == 1 in every one
// return sections as a cascade
func (s SOS) toCascade() (c cascade, err error) {
	if len(s) == 0 {
		return nil, fmt.Errorf("%w: no sections", ErrInvalidSOS)
	}
	c = make(cascade, len(s))
	for i, section := range s {
		if section[3] != 1 {
			return nil, fmt.Errorf("%w: section %d has a0 = %v, every section's a0 must be 1", ErrInvalidSOS, i, section[3])
		}
		if c[i], err = newDigitalFilter(section[:3], section[3:]); err != nil {
			return nil, err
		}
	}
	return
}

// ntaps is SciPy's tap count for sosfiltfilt's default padding:
// 2*sections + 1 minus 1 for each first-order section (b2 and a2 both zero)
func (s SOS) ntaps() (n int) {
	zeroB2, zeroA2 := 0, 0
	for _, section := range s {
		if section[2] == 0 {
			zeroB2++
		}
		if section[5] == 0 {
			zeroA2++
		}
	}
	n = 2*len(s) + 1 - min(zeroB2, zeroA2)
	return
}

// toCascadeState converts per-section states from SOS's fixed pairs to cascade's slices
func toCascadeState(zi [][2]float64) (state [][]float64) {
	state = make([][]float64, len(zi))
	for i := range zi {
		state[i] = []float64{zi[i][0], zi[i][1]}
	}
	return
}

// fromCascadeState converts per-section states from cascade's slices to SOS's fixed pairs
func fromCascadeState(state [][]float64) (zi [][2]float64) {
	zi = make([][2]float64, len(state))
	for i := range state {
		zi[i] = [2]float64{state[i][0], state[i][1]}
	}
	return
}

// SOSFilter applies the filter sos to x causally, section by section, starting from rest
// stays accurate at high orders and narrow bands, where Filter's (b, a) form loses precision
//
// returns ErrInvalidSOS if sos has no sections or any section's a0 isn't 1
//
// equivalent to scipy.signal.sosfilt(sos, x) and MATLAB's sosfilt(sos, x)
func SOSFilter(sos SOS, x []float64) (y []float64, err error) {
	c, err := sos.toCascade()
	if err != nil {
		return nil, err
	}
	y, _ = c.apply(x, nil)
	return
}

// SOSFilterState is SOSFilter starting from state zi, one pair per section, and also returns the final state zf
// passing each block's zf as the next block's zi filters a long signal in blocks exactly as filtering it in one go would
//
// returns ErrInvalidSOS as SOSFilter does, and ErrInvalidState if zi doesn't have one pair per section
//
// equivalent to scipy.signal.sosfilt(sos, x, zi=zi)
func SOSFilterState(sos SOS, x []float64, zi [][2]float64) (y []float64, zf [][2]float64, err error) {
	c, err := sos.toCascade()
	if err != nil {
		return nil, nil, err
	}
	if len(zi) != len(sos) {
		return nil, nil, fmt.Errorf("%w: zi has %d pairs but sos has %d sections", ErrInvalidState, len(zi), len(sos))
	}
	y, final := c.apply(x, toCascadeState(zi))
	zf = fromCascadeState(final)
	return
}

// SOSFilterZi returns the state each section settles into after a long unit step input
// scaled by a signal's first sample, it starts SOSFilterState without a startup transient
//
// returns ErrInvalidSOS as SOSFilter does, and ErrNoSteadyState if a section has a pole at DC
//
// equivalent to scipy.signal.sosfilt_zi(sos)
func SOSFilterZi(sos SOS) (zi [][2]float64, err error) {
	c, err := sos.toCascade()
	if err != nil {
		return nil, err
	}
	state, err := c.steadyState()
	if err != nil {
		return nil, err
	}
	zi = fromCascadeState(state)
	return
}

// SOSFiltFilt is FiltFilt for a filter given as second-order sections:
// it filters x forward then backward for zero phase, with the same edge handling options
// the default padding length is SciPy's 3*ntaps, where ntaps is 2*len(sos) + 1,
// minus one for each first-order section
//
// with MATLABPadLen the padding length is 3*(ntaps-1), same as MATLAB's filtfilt for (b, a)
// this hasn't been verified against MATLAB's filtfilt with an SOS matrix
//
// returns ErrInvalidSOS as SOSFilter does, and ErrInvalidPadding,
// ErrSignalTooShort and ErrNoSteadyState as FiltFilt does
//
// equivalent to scipy.signal.sosfiltfilt(sos, x)
func SOSFiltFilt(sos SOS, x []float64, opts ...FiltOption) (y []float64, err error) {
	c, err := sos.toCascade()
	if err != nil {
		return nil, err
	}
	y, err = c.zeroPhase(x, sos.ntaps(), opts)
	return
}
