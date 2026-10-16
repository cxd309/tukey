package conformance

// Suites is every conformance suite, run by this package's tests and by cmd/tolerances
var Suites = []Suite{
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

	// MATLAB, for the categories it can reproduce; see each SciPy suite for why
	// any without a MATLAB counterpart are missing
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
}
