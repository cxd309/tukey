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
	A      []float64 `json:"a"`      // general_cosine
	Std    float64   `json:"std"`    // gaussian
	P      float64   `json:"p"`      // general_gaussian
	Sig    float64   `json:"sig"`    // general_gaussian
	Tau    float64   `json:"tau"`    // exponential
	Center *float64  `json:"center"` // exponential; null or absent: the middle
}

// windowCall calls one window function with a vector's parameters and options
type windowCall func(p windowParams, opts []window.Option) (w []float64, err error)

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

// cosine-sum windows
// measured: within 3 ULP of SciPy, most at n = 64 symmetric,
// where SciPy's own asymmetry (up to 3.5 ULP) meets tukey's exact mirroring
// within 4 ULP of MATLAB, whose windows are symmetric, set ~10x above

var (
	hannSciPySuite            = newWindowSuite("window.Hann", reference.SourceSciPy, "hann", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, hannCall)
	hannMATLABSuite           = newWindowSuite("window.Hann", reference.SourceMATLAB, "hann", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, hannCall)
	hammingSciPySuite         = newWindowSuite("window.Hamming", reference.SourceSciPy, "hamming", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, hammingCall)
	hammingMATLABSuite        = newWindowSuite("window.Hamming", reference.SourceMATLAB, "hamming", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, hammingCall)
	blackmanSciPySuite        = newWindowSuite("window.Blackman", reference.SourceSciPy, "blackman", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, blackmanCall)
	blackmanMATLABSuite       = newWindowSuite("window.Blackman", reference.SourceMATLAB, "blackman", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, blackmanCall)
	blackmanHarrisSciPySuite  = newWindowSuite("window.BlackmanHarris", reference.SourceSciPy, "blackmanharris", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, blackmanHarrisCall)
	blackmanHarrisMATLABSuite = newWindowSuite("window.BlackmanHarris", reference.SourceMATLAB, "blackmanharris", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, blackmanHarrisCall)
	nuttallSciPySuite         = newWindowSuite("window.Nuttall", reference.SourceSciPy, "nuttall", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, nuttallCall)
	nuttallMATLABSuite        = newWindowSuite("window.Nuttall", reference.SourceMATLAB, "nuttall", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, nuttallCall)
	flatTopSciPySuite         = newWindowSuite("window.FlatTop", reference.SourceSciPy, "flattop", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, flatTopCall)
	flatTopMATLABSuite        = newWindowSuite("window.FlatTop", reference.SourceMATLAB, "flattop", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, flatTopCall)
	// MATLAB has no general_cosine or general_hamming
	// GeneralCosine measured: worst 8 ULP (1.7e-15) on HFT90D, whose large,
	// alternating coefficients nearly cancel where the window is small, set ~10x above
	generalCosineSciPySuite  = newWindowSuite("window.GeneralCosine", reference.SourceSciPy, "general_cosine", reference.Tolerance{Rel: 1e-13, Scaled: 1e-13}, generalCosineCall)
	generalHammingSciPySuite = newWindowSuite("window.GeneralHamming", reference.SourceSciPy, "general_hamming", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, generalHammingCall)
	// piecewise and sine windows
	// measured: within 3 ULP of SciPy and MATLAB, worst Tukey at alpha 0.25;
	// Boxcar and Triang are exact against both: one integer division per sample
	// (Triang) or none (Boxcar) leaves nothing to round differently, so anything
	// but exact is a bug; set ~10x above the rest
	boxcarSciPySuite    = newWindowSuite("window.Boxcar", reference.SourceSciPy, "boxcar", reference.Tolerance{}, boxcarCall)
	boxcarMATLABSuite   = newWindowSuite("window.Boxcar", reference.SourceMATLAB, "boxcar", reference.Tolerance{}, boxcarCall)
	triangSciPySuite    = newWindowSuite("window.Triang", reference.SourceSciPy, "triang", reference.Tolerance{}, triangCall)
	triangMATLABSuite   = newWindowSuite("window.Triang", reference.SourceMATLAB, "triang", reference.Tolerance{}, triangCall)
	bartlettSciPySuite  = newWindowSuite("window.Bartlett", reference.SourceSciPy, "bartlett", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, bartlettCall)
	bartlettMATLABSuite = newWindowSuite("window.Bartlett", reference.SourceMATLAB, "bartlett", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, bartlettCall)
	parzenSciPySuite    = newWindowSuite("window.Parzen", reference.SourceSciPy, "parzen", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, parzenCall)
	parzenMATLABSuite   = newWindowSuite("window.Parzen", reference.SourceMATLAB, "parzen", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, parzenCall)
	bohmanSciPySuite    = newWindowSuite("window.Bohman", reference.SourceSciPy, "bohman", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, bohmanCall)
	bohmanMATLABSuite   = newWindowSuite("window.Bohman", reference.SourceMATLAB, "bohman", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, bohmanCall)
	barthannSciPySuite  = newWindowSuite("window.Barthann", reference.SourceSciPy, "barthann", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, barthannCall)
	barthannMATLABSuite = newWindowSuite("window.Barthann", reference.SourceMATLAB, "barthann", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, barthannCall)
	tukeySciPySuite     = newWindowSuite("window.Tukey", reference.SourceSciPy, "tukey", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, tukeyCall)
	tukeyMATLABSuite    = newWindowSuite("window.Tukey", reference.SourceMATLAB, "tukey", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, tukeyCall)
	hanningSciPySuite   = newWindowSuite("window.Hanning", reference.SourceSciPy, "hanning", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, hanningCall)
	hanningMATLABSuite  = newWindowSuite("window.Hanning", reference.SourceMATLAB, "hanning", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, hanningCall)

	// MATLAB has no cosine or lanczos
	cosineSciPySuite  = newWindowSuite("window.Cosine", reference.SourceSciPy, "cosine", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, cosineCall)
	lanczosSciPySuite = newWindowSuite("window.Lanczos", reference.SourceSciPy, "lanczos", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, lanczosCall)

	// exponential family
	// measured: within 1 ULP of SciPy and MATLAB (gausswin, alpha converted from std);
	// exp of a single well-conditioned argument leaves little to differ, set ~10x above
	gaussianSciPySuite  = newWindowSuite("window.Gaussian", reference.SourceSciPy, "gaussian", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, gaussianCall)
	gaussianMATLABSuite = newWindowSuite("window.Gaussian", reference.SourceMATLAB, "gaussian", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, gaussianCall)

	// MATLAB has no general_gaussian or exponential
	generalGaussianSciPySuite = newWindowSuite("window.GeneralGaussian", reference.SourceSciPy, "general_gaussian", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, generalGaussianCall)
	exponentialSciPySuite     = newWindowSuite("window.Exponential", reference.SourceSciPy, "exponential", reference.Tolerance{Rel: 1e-14, Scaled: 1e-14}, exponentialCall)
)

func hannCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Hann(p.N, opts...)
	return
}

func hammingCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Hamming(p.N, opts...)
	return
}

func blackmanCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Blackman(p.N, opts...)
	return
}

func blackmanHarrisCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.BlackmanHarris(p.N, opts...)
	return
}

func nuttallCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Nuttall(p.N, opts...)
	return
}

func flatTopCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.FlatTop(p.N, opts...)
	return
}

func generalCosineCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.GeneralCosine(p.N, p.A, opts...)
	return
}

func generalHammingCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.GeneralHamming(p.N, p.Alpha, opts...)
	return
}
func boxcarCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Boxcar(p.N, opts...)
	return
}

func triangCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Triang(p.N, opts...)
	return
}

func bartlettCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Bartlett(p.N, opts...)
	return
}

func parzenCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Parzen(p.N, opts...)
	return
}

func bohmanCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Bohman(p.N, opts...)
	return
}

func barthannCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Barthann(p.N, opts...)
	return
}

func tukeyCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Tukey(p.N, p.Alpha, opts...)
	return
}

func hanningCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Hanning(p.N, opts...)
	return
}

func cosineCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Cosine(p.N, opts...)
	return
}

func lanczosCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Lanczos(p.N, opts...)
	return
}
func gaussianCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.Gaussian(p.N, p.Std, opts...)
	return
}

func generalGaussianCall(p windowParams, opts []window.Option) (w []float64, err error) {
	w, err = window.GeneralGaussian(p.N, p.P, p.Sig, opts...)
	return
}

func exponentialCall(p windowParams, opts []window.Option) (w []float64, err error) {
	if p.Center != nil {
		opts = append(opts, window.Center(*p.Center))
	}
	w, err = window.Exponential(p.N, p.Tau, opts...)
	return
}
