package fft_test

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"

	"github.com/cxd309/tukey/fft"
)

func ExampleForwardReal() {
	// one second of a 50 Hz tone, sampled at 1 kHz
	const fs, n = 1000.0, 1000
	x := make([]float64, n)
	for i := range x {
		x[i] = math.Sin(2 * math.Pi * 50 * float64(i) / fs)
	}

	spectrum, err := fft.ForwardReal(x)
	if err != nil {
		log.Fatal(err)
	}
	freqs := fft.RealFreq(n, 1/fs) // d = 1/fs gives Hz

	peak := 0
	for k := range spectrum {
		if cmplx.Abs(spectrum[k]) > cmplx.Abs(spectrum[peak]) {
			peak = k
		}
	}
	fmt.Printf("peak at %.0f Hz\n", freqs[peak])
	// Output:
	// peak at 50 Hz
}

func ExampleShift() {
	// Forward's bins run from 0 up, then the negative frequencies;
	// Shift puts them in ascending order, with 0 in the centre
	freqs := fft.Freq(8, 1.0/8)
	fmt.Println(freqs)
	fmt.Println(fft.Shift(freqs))
	// Output:
	// [0 1 2 3 -4 -3 -2 -1]
	// [-4 -3 -2 -1 0 1 2 3]
}
