package conformance

import (
	"strings"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/internal/reference"
)

// lowCutoffFilters are the (b, a) filters whose cutoffs are so low that their
// poles sit next to z = 1. They're the original MATLAB app's own designs.
// Rounding differences between implementations are magnified by roughly the
// square of 1/(1 - pole), so these filters agree with SciPy and MATLAB far less
// closely than the rest. Their vectors get suites of their own, so their looser
// tolerances don't loosen the checks on every other filter.
var lowCutoffFilters = []string{"butter2_lp0.001", "butter2_hp0.0008"}

// isLowCutoff reports whether a vector file is for one of lowCutoffFilters.
// Vector file names start with the filter's name, followed by "_" or ".json".
func isLowCutoff(file string) (ok bool) {
	for _, filter := range lowCutoffFilters {
		if strings.HasPrefix(file, filter+"_") || strings.HasPrefix(file, filter+".") {
			return true
		}
	}
	return false
}

// isWellConditioned reports whether a vector file is for any filter but
// lowCutoffFilters.
func isWellConditioned(file string) (ok bool) {
	ok = !isLowCutoff(file)
	return
}

// The low-cutoff suites check the same functions as their namesakes, on
// lowCutoffFilters only. Each is set ~10x above its worst on either architecture.
var (
	// Filter measured: worst 9.5e-13 against SciPy and 6.2e-13 against MATLAB,
	// both on the 0.0008 high-pass rejecting a constant input, so its output
	// cancels to ~0.
	filterLowCutoffSciPySuite  = newSuiteOf("Filter (low cutoff)", reference.SourceSciPy, "filter", tolerance(1e-11), filterOutputs, isLowCutoff)
	filterLowCutoffMATLABSuite = newSuiteOf("Filter (low cutoff)", reference.SourceMATLAB, "filter", tolerance(1e-11), filterOutputs, isLowCutoff)

	// FilterState measured: worst 9.9e-14 against SciPy and 1.1e-13 against MATLAB.
	filterStateLowCutoffSciPySuite  = newSuiteOf("FilterState (low cutoff)", reference.SourceSciPy, "lfilter_state", tolerance(1e-12), filterStateOutputs, isLowCutoff)
	filterStateLowCutoffMATLABSuite = newSuiteOf("FilterState (low cutoff)", reference.SourceMATLAB, "lfilter_state", tolerance(1e-11), filterStateOutputs, isLowCutoff)

	// FiltFilt measured: worst 4.4e-12 against SciPy and 1.5e-11 against MATLAB
	// (with MATLABPadLen), both on the 0.001 low-pass. Its starting state solves a
	// nearly singular system, and SciPy and MATLAB differ from each other by as much.
	filtfiltLowCutoffSciPySuite  = newSuiteOf("FiltFilt (low cutoff)", reference.SourceSciPy, "filtfilt", tolerance(1e-10), filtfiltOutputs(), isLowCutoff)
	filtfiltLowCutoffMATLABSuite = newSuiteOf("FiltFilt (low cutoff)", reference.SourceMATLAB, "filtfilt", tolerance(1e-9), filtfiltOutputs(dsp.MATLABPadLen()), isLowCutoff)

	// Freqz measured: worst 1.5e-12 against SciPy and 1.1e-11 against MATLAB.
	// At DC the 0.001 low-pass's gain is exactly 1, and tukey gives 1 + 1.6e-12
	// where MATLAB gives 1 + 1.3e-11, so much of the difference is MATLAB's.
	freqzLowCutoffSciPySuite  = newSuiteOf("Freqz (low cutoff)", reference.SourceSciPy, "freqz", tolerance(1e-10), freqzOutputs, isLowCutoff)
	freqzLowCutoffMATLABSuite = newSuiteOf("Freqz (low cutoff)", reference.SourceMATLAB, "freqz", tolerance(1e-9), freqzOutputs, isLowCutoff)
)
