package conformance

import (
	"fmt"

	"github.com/cxd309/tukey/fft"
	"github.com/cxd309/tukey/internal/reference"
)

type fftOptionVector struct {
	reference.Meta
	Params struct {
		Transform string            `json:"transform"`
		X         reference.Complex `json:"x"`
		Norm      string            `json:"norm"`
		N         *int              `json:"n"` // N for fft/ifft/rfft (null: unset); the length for irfft
	} `json:"params"`
	Output struct {
		Y reference.Complex `json:"y"`
	} `json:"output"`
}

// fftOptionsSuite checks Norm and N on each transform against scipy.fft
// measured: worst 3 ULP (6.4e-16, ifft at n = 1009 with norm='forward'),
// no worse than the transforms themselves, set to match them
var fftOptionsSuite = newSuite("fft options", reference.SourceSciPy, "fft_options",
	reference.Tolerance{Rel: 1e-14, Scaled: 1e-14},
	func(v fftOptionVector) (outputs []Output, err error) {
		x, err := v.Params.X.Values()
		if err != nil {
			return nil, err
		}
		want, err := v.Output.Y.Values()
		if err != nil {
			return nil, err
		}
		norm, err := normalisation(v.Params.Norm)
		if err != nil {
			return nil, err
		}
		opts := []fft.Option{fft.Norm(norm)}
		if v.Params.N != nil && v.Params.Transform != "irfft" {
			opts = append(opts, fft.N(*v.Params.N))
		}

		var y []complex128
		switch v.Params.Transform {
		case "fft":
			y, err = fft.Forward(x, opts...)
		case "ifft":
			y, err = fft.Inverse(x, opts...)
		case "rfft":
			y, err = fft.ForwardReal(realParts(x), opts...)
		case "irfft":
			var r []float64
			r, err = fft.InverseReal(x, *v.Params.N, opts...)
			y = asComplex(r)
		default:
			err = fmt.Errorf("unknown transform %q", v.Params.Transform)
		}
		if err != nil {
			return nil, err
		}
		outputs = []Output{{Name: "y", GotComplex: y, WantComplex: want}}
		return
	})

// normalisation translates SciPy's norm names
func normalisation(name string) (norm fft.Normalisation, err error) {
	switch name {
	case "backward":
		norm = fft.NormBackward
	case "ortho":
		norm = fft.NormOrtho
	case "forward":
		norm = fft.NormForward
	default:
		err = fmt.Errorf("unknown norm %q", name)
	}
	return
}

// realParts and asComplex convert between real values and complex ones with zero
// imaginary parts, as fft_options stores real inputs and outputs
func realParts(x []complex128) (r []float64) {
	r = make([]float64, len(x))
	for i, v := range x {
		r[i] = real(v)
	}
	return
}

func asComplex(r []float64) (x []complex128) {
	x = make([]complex128, len(r))
	for i, v := range r {
		x[i] = complex(v, 0)
	}
	return
}

type fftFreqVector struct {
	reference.Meta
	Params struct {
		N int     `json:"n"`
		D float64 `json:"d"`
	} `json:"params"`
	Output struct {
		Freq  []float64 `json:"freq"`
		RFreq []float64 `json:"rfreq"`
	} `json:"output"`
}

// fftFreqSuite checks fft.Freq and fft.RealFreq against scipy.fft.fftfreq and rfftfreq
// measured: exact
// Each value is an integer times 1/(n·d), computed as NumPy does,
// and IEEE 754 rounds a single division and multiplication identically everywhere,
// so anything but an exact match is a bug
var fftFreqSuite = newSuite("fft.Freq", reference.SourceSciPy, "fftfreq",
	reference.Tolerance{},
	func(v fftFreqVector) (outputs []Output, err error) {
		outputs = []Output{
			{Name: "freq", Got: fft.Freq(v.Params.N, v.Params.D), Want: v.Output.Freq},
			{Name: "rfreq", Got: fft.RealFreq(v.Params.N, v.Params.D), Want: v.Output.RFreq},
		}
		return
	})

type fftShiftVector struct {
	reference.Meta
	Params struct {
		X []float64 `json:"x"`
	} `json:"params"`
	Output struct {
		Shift  []float64 `json:"shift"`
		IShift []float64 `json:"ishift"`
	} `json:"output"`
}

// fftShiftSuite checks fft.Shift and fft.InverseShift against scipy.fft.fftshift and ifftshift
// they only reorder values, so they must match exactly
var fftShiftSuite = newSuite("fft.Shift", reference.SourceSciPy, "fftshift",
	reference.Tolerance{},
	func(v fftShiftVector) (outputs []Output, err error) {
		outputs = []Output{
			{Name: "shift", Got: fft.Shift(v.Params.X), Want: v.Output.Shift},
			{Name: "ishift", Got: fft.InverseShift(v.Params.X), Want: v.Output.IShift},
		}
		return
	})
