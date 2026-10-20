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

// windowParams is every window's parameters; each window reads its own
type windowParams struct {
	N     int       `json:"n"`
	Sym   bool      `json:"sym"`
	Alpha float64   `json:"alpha"` // general_hamming
	A     []float64 `json:"a"`     // general_cosine
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
