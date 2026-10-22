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
		filterLowCutoffSciPySuite,
		filterStateSciPySuite,
		filterStateLowCutoffSciPySuite,
		filterZiSciPySuite,
		filtfiltSciPySuite,
		filtfiltLowCutoffSciPySuite,
		freqzSciPySuite,
		freqzLowCutoffSciPySuite,
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
		detrendSciPySuite,
		interp1SciPySuite,
		fillMissingSciPySuite,
		rmsSciPySuite,
		findPeaksSciPySuite,
		peakProminencesSciPySuite,
		peakWidthsSciPySuite,
		hilbertSciPySuite,
		envelopeSciPySuite,
		findSignalDirectSciPySuite,
		findSignalFastSciPySuite,
		findSignalHybridSciPySuite,
	},
	windowSuites(reference.SourceSciPy),

	// MATLAB, for the categories it can reproduce
	// see each SciPy suite for why any without a MATLAB counterpart are missing
	[]Suite{
		butterMATLABSuite,
		filterMATLABSuite,
		filterLowCutoffMATLABSuite,
		filterStateMATLABSuite,
		filterStateLowCutoffMATLABSuite,
		filtfiltMATLABSuite,
		filtfiltLowCutoffMATLABSuite,
		freqzMATLABSuite,
		freqzLowCutoffMATLABSuite,
		sosFilterMATLABSuite,
		sosFreqzMATLABSuite,
		forwardMATLABSuite,
		inverseMATLABSuite,
		forwardRealMATLABSuite,
		inverseRealMATLABSuite,
		fftOptionsMATLABSuite,
		fftShiftMATLABSuite,
		detrendMATLABSuite,
		interp1MATLABSuite,
		fillMissingMATLABSuite,
		rmsMATLABSuite,
		findPeaksMATLABSuite,
		peakProminencesMATLABSuite,
		peakWidthsMATLABSuite,
		hilbertMATLABSuite,
		envelopeMATLABSuite,
		findSignalDirectMATLABSuite,
		findSignalFastMATLABSuite,
		findSignalHybridMATLABSuite,
	},
	windowSuites(reference.SourceMATLAB),
)
