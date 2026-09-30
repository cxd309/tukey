package conformance

// Suites is every conformance suite, run by this package's tests and by cmd/tolerances
var Suites = []Suite{
	butterSuite,
	butterSOSSuite,
	filterSuite,
	filtfiltSuite,
	sosFilterSuite,
	sosFilterZiSuite,
	sosFiltFiltSuite,
	zpk2sosSuite,
	filterStateSuite,
	filterZiSuite,
}
