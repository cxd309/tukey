package conformance

// Suites is every conformance suite, run by this package's tests and by cmd/tolerances
var Suites = []Suite{
	butterSuite,
	butterSOSSuite,
	filterSuite,
	filterStateSuite,
	filterZiSuite,
	filtfiltSuite,
	freqzSuite,
	sosFilterSuite,
	sosFilterZiSuite,
	sosFiltFiltSuite,
	sosFreqzSuite,
	zpk2sosSuite,
}
