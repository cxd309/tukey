package conformance

import (
	"slices"

	"github.com/cxd309/tukey/internal/reference"
)

// Suites is every conformance suite, run by this package's tests and by cmd/tolerances
// SciPy's first, then MATLAB's.
// The window suites come from windowChecks
var Suites = slices.Concat(
	[]Suite{
		butterSciPySuite,
		butterSOSSciPySuite,
		filterSciPySuite,
		filterStateSciPySuite,
		filterZiSciPySuite,
		filtfiltSciPySuite,
		freqzSciPySuite,
		sosFilterSciPySuite,
		sosFilterZiSciPySuite,
		sosFiltFiltSciPySuite,
		sosFreqzSciPySuite,
		zpk2sosSciPySuite,
		forwardSciPySuite,
		inverseSciPySuite,
		forwardRealSciPySuite,
		inverseRealSciPySuite,
		fftOptionsSciPySuite,
		fftFreqSciPySuite,
		fftShiftSciPySuite,
	},
	windowSuites(reference.SourceSciPy),

	// MATLAB, for the categories it can reproduce
	// see each SciPy suite for why any without a MATLAB counterpart are missing
	[]Suite{
		butterMATLABSuite,
		filterMATLABSuite,
		filterStateMATLABSuite,
		filtfiltMATLABSuite,
		freqzMATLABSuite,
		sosFilterMATLABSuite,
		sosFreqzMATLABSuite,
		forwardMATLABSuite,
		inverseMATLABSuite,
		forwardRealMATLABSuite,
		inverseRealMATLABSuite,
		fftOptionsMATLABSuite,
		fftShiftMATLABSuite,
	},
	windowSuites(reference.SourceMATLAB),
)
