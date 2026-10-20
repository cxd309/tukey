package window

import (
	"math"
	"math/cmplx"

	"github.com/cxd309/tukey/fft"
)

// Chebwin returns the n-point Dolph–Chebyshev window:
// for its main-lobe width it has the lowest possible sidelobes,
// all equal, at attenuation dB below the main lobe.
// Normalised so its peak is 1
//
// SciPy warns that below about 45 dB it isn't suitable for spectral analysis
//
// returns ErrInvalidLength if n < 0, and ErrInvalidParameter unless attenuation > 0
//
// equivalent to scipy.signal.windows.chebwin(n, attenuation),
// equivalent to MATLAB's chebwin(n, attenuation)
func Chebwin(n int, attenuation float64, opts ...Option) (w []float64, err error) {
	cfg, err := newConfig("Chebwin", opts, optionSet{optionPeriodic: true})
	if err != nil {
		return nil, err
	}
	if err = positive("Chebwin", "attenuation", attenuation); err != nil {
		return nil, err
	}
	if w, short, err := shortWindow(n); short {
		return w, err
	}
	m := n
	if cfg.periodic {
		m = n + 1
	}
	w = chebwin(m, attenuation)[:n]
	return
}

// chebwin is the symmetric Chebyshev window of length m >= 2:
// the DFT of the Chebyshev polynomial T_{m-1}
// sampled around a circle scaled so its sidelobes sit at attenuation dB,
// mirrored about its centre and normalised to peak at 1
func chebwin(m int, attenuation float64) (w []float64) {
	order := m - 1
	beta := math.Cosh(math.Acosh(math.Pow(10, attenuation/20)) / float64(order))
	p := make([]complex128, m)
	for k := range p {
		x := beta * math.Cos(math.Pi*float64(k)/float64(m))
		p[k] = complex(chebyshevT(order, x), 0)
		if m%2 == 0 {
			// half-sample shift, so the even-length window is centred between samples
			p[k] *= cmplx.Exp(complex(0, math.Pi*float64(k)/float64(m)))
		}
	}
	spectrum, err := fft.Forward(p)
	if err != nil {
		panic(err) // p has m >= 2 points and no options: Forward can't fail
	}

	// the window is the spectrum's real part, mirrored:
	// odd m keeps bin 0 at its centre,
	// even m drops it and repeats bins 1...m/2 on both sides
	w = make([]float64, m)
	centre := (m - 1) / 2
	if m%2 == 1 {
		for j := range centre + 1 {
			w[centre-j], w[centre+j] = real(spectrum[j]), real(spectrum[j])
		}
	} else {
		for j := 1; j <= m/2; j++ {
			w[centre+1-j], w[centre+j] = real(spectrum[j]), real(spectrum[j])
		}
	}

	peak := w[0]
	for _, v := range w {
		peak = math.Max(peak, v)
	}
	for i := range w {
		w[i] /= peak
	}
	return
}

// chebyshevT is the Chebyshev polynomial of the first kind T_order(x),
// from its trigonometric form inside [-1, 1] and hyperbolic form outside,
// same as SciPy (its power series loses accuracy).
// T_order(-x) = (-1)^order T_order(x)
func chebyshevT(order int, x float64) (t float64) {
	n := float64(order)
	switch {
	case x > 1:
		t = math.Cosh(n * math.Acosh(x))
	case x < -1:
		t = math.Cosh(n * math.Acosh(-x))
		if order%2 == 1 {
			t = -t
		}
	default:
		t = math.Cos(n * math.Acos(x))
	}
	return
}
