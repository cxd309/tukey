package dsp_test

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"slices"

	"github.com/cxd309/tukey/dsp"
	"github.com/cxd309/tukey/window"
)

func ExampleButter() {
	// a 2nd-order lowpass at 150Hz, for a signal sampled at 1000 Hz
	f, err := dsp.Butter(2, dsp.Lowpass(150.0/(1000.0/2)))
	if err != nil {
		log.Fatal(err)
	}
	b, a, err := f.BA()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("b = %.4f\na = %.4f\n", b, a)
	// Output:
	// b = [0.1311 0.2622 0.1311]
	// a = [1.0000 -0.7478 0.2722]
}

func ExampleFiltFilt() {
	// a narrow pulse at sample 50
	x := make([]float64, 100)
	x[50] = 1

	f, err := dsp.Butter(2, dsp.Lowpass(0.1))
	if err != nil {
		log.Fatal(err)
	}
	b, a, err := f.BA()
	if err != nil {
		log.Fatal(err)
	}
	causal, err := dsp.Filter(b, a, x)
	if err != nil {
		log.Fatal(err)
	}
	zeroPhase, err := dsp.FiltFilt(b, a, x)
	if err != nil {
		log.Fatal(err)
	}

	// Filter delays the pulse, and FiltFilt smooths it but leaves it where it was
	fmt.Println("Filter peak at sample", peak(causal))
	fmt.Println("FiltFilt peak at sample", peak(zeroPhase))
	// Output:
	// Filter peak at sample 54
	// FiltFilt peak at sample 50
}

func ExampleHzToNyquist() {
	// a 4th-order lowpass at 100 Hz, for a signal sampled at 1 kHz
	const fs = 1000.0
	wn, err := dsp.HzToNyquist(100.0, fs)
	if err != nil {
		log.Fatal(err)
	}
	f, err := dsp.Butter(4, dsp.Lowpass(wn))
	if err != nil {
		log.Fatal(err)
	}
	sos, err := f.SOS()
	if err != nil {
		log.Fatal(err)
	}

	// its response, labelled in Hz, is -3 dB at the cutoff
	w, h, err := dsp.SOSFreqz(sos, 500)
	if err != nil {
		log.Fatal(err)
	}
	hz, err := dsp.RadPerSampleToHz(w, fs)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("wn = %g\n", wn)
	fmt.Printf("%.0f Hz: %.2f dB\n", hz[100], 20*math.Log10(cmplx.Abs(h[100])))
	// Output:
	// wn = 0.2
	// 100 Hz: -3.01 dB
}

func ExampleWelch() {
	// two tones, at 62.5 Hz and 187.5 Hz, sampled at 1 kHz for 4 s
	const fs = 1000.0
	x := make([]float64, 4000)
	for i := range x {
		t := float64(i) / fs
		x[i] = math.Sin(2*math.Pi*62.5*t) + 0.5*math.Sin(2*math.Pi*187.5*t)
	}

	// ScalingSpectrum reads each tone's power, A²/2, at its peak
	w, err := window.Hann(256, window.Periodic())
	if err != nil {
		log.Fatal(err)
	}
	spectrum, err := dsp.Welch(x, fs, w, dsp.ScalingSpectrum())
	if err != nil {
		log.Fatal(err)
	}
	for k, f := range spectrum.Frequencies {
		if f == 62.5 || f == 187.5 {
			fmt.Printf("%.1f Hz: %.3f\n", f, spectrum.Power[k])
		}
	}
	// Output:
	// 62.5 Hz: 0.500
	// 187.5 Hz: 0.125
}

func ExampleFindPeaks() {
	x := []float64{0, 2, 1, 3, 2.5, 2.8, 0, 4, 1, 0}

	// every local maximum, then only those standing at least 1.5 above their surroundings
	all, err := dsp.FindPeaks(x)
	if err != nil {
		log.Fatal(err)
	}
	prominent, err := dsp.FindPeaks(x, dsp.MinPeakProminence(1.5))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("all: %+v\n", all)
	fmt.Printf("prominent: %+v\n", prominent)
	// Output:
	// all: [{Index:1 Height:2} {Index:3 Height:3} {Index:5 Height:2.8} {Index:7 Height:4}]
	// prominent: [{Index:3 Height:3} {Index:7 Height:4}]
}

func ExampleEnvelope() {
	// a 50 Hz tone whose amplitude swings between 0.5 and 1.5 twice a second, sampled at 1 kHz for 1 s
	const fs = 1000.0
	x := make([]float64, 1000)
	for i := range x {
		t := float64(i) / fs
		x[i] = (1 + 0.5*math.Cos(2*math.Pi*2*t)) * math.Cos(2*math.Pi*50*t)
	}

	upper, _, err := dsp.Envelope(x, dsp.EnvelopeAnalytic())
	if err != nil {
		log.Fatal(err)
	}
	for _, i := range []int{0, 125, 250} {
		fmt.Printf("t = %.3f s: %.2f\n", float64(i)/fs, upper[i])
	}
	// Output:
	// t = 0.000 s: 1.50
	// t = 0.125 s: 1.00
	// t = 0.250 s: 0.50
}

func ExamplePercentile() {
	x := []float64{12, 3, 7, 1, 9, 15, 4, 8, 2, 10}

	median, err := dsp.Percentile(x, 50)
	if err != nil {
		log.Fatal(err)
	}
	upper, err := dsp.Percentile(x, 75)
	if err != nil {
		log.Fatal(err)
	}
	// MATLAB's prctile places a percentile differently between the sorted samples
	upperHazen, err := dsp.Percentile(x, 75, dsp.PercentileHazen())
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("median:", median)
	fmt.Println("75th:", upper)
	fmt.Println("75th, same as MATLAB:", upperHazen)
	// Output:
	// median: 7.5
	// 75th: 9.75
	// 75th, same as MATLAB: 10
}

func ExampleConvolve() {
	// a 101-tap lowpass at a fifth of Nyquist, by the window method
	w, err := window.Hamming(101)
	if err != nil {
		log.Fatal(err)
	}
	h, err := dsp.FIRWin(w, dsp.Lowpass(0.2))
	if err != nil {
		log.Fatal(err)
	}

	// filter a slow tone with a fast one on top, keeping x's length
	x := make([]float64, 1000)
	slow := make([]float64, len(x))
	for i := range x {
		slow[i] = math.Sin(math.Pi * 0.02 * float64(i))
		x[i] = slow[i] + math.Sin(math.Pi*0.6*float64(i))
	}
	y, err := dsp.Convolve(x, h, dsp.ConvolveSame())
	if err != nil {
		log.Fatal(err)
	}

	// away from the ends, only the slow tone is left, less the passband's ripple
	worst := 0.0
	for i := 100; i < 900; i++ {
		worst = math.Max(worst, math.Abs(y[i]-slow[i]))
	}
	fmt.Printf("y matches the slow tone to within %.0e\n", worst)
	// Output:
	// y matches the slow tone to within 3e-03
}

func ExampleFindSignal() {
	// a short template hidden in a longer recording
	recording := []float64{0, 1, 0, -1, 0, 2, 5, 3, -2, 0, 1, 0}
	template := []float64{2, 5, 3}

	match, err := dsp.FindSignal(recording, template)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("samples %d to %d, distance %g\n", match.Start, match.Stop, match.Distance)
	// Output:
	// samples 5 to 7, distance 0
}

// peak is the index of v's largest value.
func peak(v []float64) (index int) {
	index = slices.Index(v, slices.Max(v))
	return
}
