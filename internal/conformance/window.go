package conformance

import (
	"github.com/cxd309/tukey/internal/reference"
	"github.com/cxd309/tukey/window"
)

type windowVector struct {
	reference.Meta
	Params windowParams `json:"params"`
	Output struct {
		W []float64 `json:"w"`
	} `json:"output"`
}

// windowParams is every window's parameters, each window reads its own
type windowParams struct {
	N      int       `json:"n"`
	Sym    bool      `json:"sym"`
	Alpha  float64   `json:"alpha"`  // tukey, general_hamming
	Beta   float64   `json:"beta"`   // kaiser, kaiser_bessel_derived
	A      []float64 `json:"a"`      // general_cosine
	Std    float64   `json:"std"`    // gaussian
	P      float64   `json:"p"`      // general_gaussian
	Sig    float64   `json:"sig"`    // general_gaussian
	Tau    float64   `json:"tau"`    // exponential
	Center *float64  `json:"center"` // exponential. null or absent: the middle
	Nbar   int       `json:"nbar"`   // taylor
	Sll    float64   `json:"sll"`    // taylor, positive dB as SciPy
	Norm   bool      `json:"norm"`   // taylor
	At     float64   `json:"at"`     // chebwin
}

// windowCall calls one window function with a vector's parameters and options
type windowCall func(p windowParams, opts []window.Option) (w []float64, err error)

// windowCheck is one window function, checked against testdata/<source>/window/<category>
// from each source it has; matlab is nil where MATLAB has no equivalent
type windowCheck struct {
	name, category string
	call           windowCall
	scipy          reference.Tolerance
	matlab         *reference.Tolerance
}

// windowChecks is every window, grouped by family, each with its measured tolerances
var windowChecks = []windowCheck{
	// cosine-sum windows
	// measured: within 3 ULP of SciPy, most at n = 64 symmetric, where SciPy's
	// own asymmetry (up to 3.5 ULP) meets tukey's exact mirroring; within 4 ULP
	// of MATLAB, whose windows are symmetric; set ~10x above
	{name: "window.Hann", category: "hann", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Hann(p.N, opts...) }},
	{name: "window.Hamming", category: "hamming", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Hamming(p.N, opts...) }},
	{name: "window.Blackman", category: "blackman", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Blackman(p.N, opts...) }},
	{name: "window.BlackmanHarris", category: "blackmanharris", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			return window.BlackmanHarris(p.N, opts...)
		}},
	{name: "window.Nuttall", category: "nuttall", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Nuttall(p.N, opts...) }},
	{name: "window.FlatTop", category: "flattop", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.FlatTop(p.N, opts...) }},
	// GeneralCosine measured: worst 8 ULP (1.7e-15) on HFT90D, whose large,
	// alternating coefficients nearly cancel where the window is small; set ~10x above
	// MATLAB has no general_cosine or general_hamming
	{name: "window.GeneralCosine", category: "general_cosine", scipy: tolerance(1e-13),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			return window.GeneralCosine(p.N, p.A, opts...)
		}},
	{name: "window.GeneralHamming", category: "general_hamming", scipy: tolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			return window.GeneralHamming(p.N, p.Alpha, opts...)
		}},

	// piecewise and sine windows
	// measured: within 3 ULP of SciPy and MATLAB, worst Tukey at alpha 0.25;
	// Boxcar and Triang are exact against both: one integer division per sample
	// (Triang) or none (Boxcar) leaves nothing to round differently, so anything
	// but exact is a bug; set ~10x above the rest
	{name: "window.Boxcar", category: "boxcar", scipy: tolerance(0), matlab: matlabTolerance(0),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Boxcar(p.N, opts...) }},
	{name: "window.Triang", category: "triang", scipy: tolerance(0), matlab: matlabTolerance(0),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Triang(p.N, opts...) }},
	{name: "window.Bartlett", category: "bartlett", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Bartlett(p.N, opts...) }},
	{name: "window.Parzen", category: "parzen", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Parzen(p.N, opts...) }},
	{name: "window.Bohman", category: "bohman", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Bohman(p.N, opts...) }},
	{name: "window.Barthann", category: "barthann", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Barthann(p.N, opts...) }},
	{name: "window.Tukey", category: "tukey", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			return window.Tukey(p.N, p.Alpha, opts...)
		}},
	{name: "window.Hanning", category: "hanning", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Hanning(p.N, opts...) }},
	// MATLAB has no cosine or lanczos
	{name: "window.Cosine", category: "cosine", scipy: tolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Cosine(p.N, opts...) }},
	{name: "window.Lanczos", category: "lanczos", scipy: tolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) { return window.Lanczos(p.N, opts...) }},

	// exponential family
	// measured: within 1 ULP of SciPy and MATLAB (gausswin, alpha converted from std);
	// exp of a single well-conditioned argument leaves little to differ; set ~10x above
	{name: "window.Gaussian", category: "gaussian", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			return window.Gaussian(p.N, p.Std, opts...)
		}},
	// MATLAB has no general_gaussian or exponential
	{name: "window.GeneralGaussian", category: "general_gaussian", scipy: tolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			return window.GeneralGaussian(p.N, p.P, p.Sig, opts...)
		}},
	{name: "window.Exponential", category: "exponential", scipy: tolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			if p.Center != nil {
				opts = append(opts, window.Center(*p.Center))
			}
			return window.Exponential(p.N, p.Tau, opts...)
		}},

	// Bessel windows
	// Kaiser measured: worst 8 ULP against SciPy (beta 14) and 9 against MATLAB
	// (beta 8.6): I0 grows like exp, so its relative error is about its argument's
	// absolute error, which scales with beta; set ~10x above
	{name: "window.Kaiser", category: "kaiser", scipy: tolerance(1e-13), matlab: matlabTolerance(1e-13),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			return window.Kaiser(p.N, p.Beta, opts...)
		}},
	// KaiserBesselDerived measured: exact against SciPy, as tukey's I0 (Cephes, as
	// SciPy's) matches it to the bit at these arguments; that's coincidence, not
	// construction, so it takes Kaiser's tolerance, whose error it inherits
	// MATLAB has no kaiser_bessel_derived; its vectors are symmetric and even n only
	{name: "window.KaiserBesselDerived", category: "kaiser_bessel_derived", scipy: tolerance(1e-13),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			return window.KaiserBesselDerived(p.N, p.Beta, opts...)
		}},

	// Taylor and Chebyshev
	// Taylor measured: within 2 ULP of SciPy and MATLAB (its raw cases); set ~10x above
	{name: "window.Taylor", category: "taylor", scipy: tolerance(1e-14), matlab: matlabTolerance(1e-14),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			if !p.Norm {
				opts = append(opts, window.Unnormalised())
			}
			return window.Taylor(p.N, p.Nbar, p.Sll, opts...)
		}},
	// Chebwin measured: worst 1.5e-12 against SciPy (50 dB, n = 255 periodic) and
	// 3.7e-13 against MATLAB, as far apart as SciPy and MATLAB are (3.2e-13):
	// near x = ±1, T_{n-1}(x) magnifies an ULP of difference in x by about
	// (n-1)/sqrt(1-x^2), and Go, NumPy and MATLAB round x differently; the FFT
	// then spreads that evenly over the window; set ~10x above
	{name: "window.Chebwin", category: "chebwin", scipy: tolerance(1e-11), matlab: matlabTolerance(1e-11),
		call: func(p windowParams, opts []window.Option) ([]float64, error) {
			return window.Chebwin(p.N, p.At, opts...)
		}},
}

// windowSuites is a suite for each window in windowChecks with vectors from source,
// in windowChecks' order
func windowSuites(source reference.Source) (suites []Suite) {
	for _, c := range windowChecks {
		tol, ok := c.tolerance(source)
		if !ok {
			continue
		}
		suites = append(suites, newWindowSuite(c.name, source, c.category, tol, c.call))
	}
	return
}

// tolerance is c's tolerance against source, and whether source has vectors for it
func (c windowCheck) tolerance(source reference.Source) (tol reference.Tolerance, ok bool) {
	switch {
	case source == reference.SourceSciPy:
		return c.scipy, true
	case source == reference.SourceMATLAB && c.matlab != nil:
		return *c.matlab, true
	}
	return
}

// newWindowSuite builds a suite checking call against testdata/<source>/window/<category>
func newWindowSuite(name string, source reference.Source, category string, tol reference.Tolerance, call windowCall) (s Suite) {
	s = newSuite(name, source, "window/"+category, tol, func(v windowVector) (outputs []Output, err error) {
		var opts []window.Option
		if !v.Params.Sym {
			opts = append(opts, window.Periodic())
		}
		w, err := call(v.Params, opts)
		if err != nil {
			return nil, err
		}
		// Scale 1: windows peak around 1, and some are all zeros (Hann at n = 2),
		// which relative to their own scale would have to match exactly
		outputs = []Output{{Name: "w", Got: w, Want: v.Output.W, Scale: 1}}
		return
	})
	return
}

// tolerance is a Tolerance bounding both relative and scaled error by bound;
// 0 is exact
func tolerance(bound float64) (tol reference.Tolerance) {
	tol = reference.Tolerance{Rel: bound, Scaled: bound}
	return
}

// matlabTolerance is tolerance(bound), for windowCheck's optional matlab field
func matlabTolerance(bound float64) (tol *reference.Tolerance) {
	t := tolerance(bound)
	tol = &t
	return
}
