package dsp

import (
	"encoding/json"
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

type butterVector struct {
	Description  string `json:"description"`
	SciPyVersion string `json:"scipy_version"`
	NumPyVersion string `json:"numpy_version"`
	Params       struct {
		Order int             `json:"order"`
		Wn    json.RawMessage `json:"wn"`
		BType string          `json:"btype"`
	} `json:"params"`
	Output struct {
		B []float64 `json:"b"`
		A []float64 `json:"a"`
	} `json:"output"`
}

// parseWn handles the golden vectors' wn field, which is a scalar for
// low/high-pass and a [low, high] pair for band-pass/stop.
func parseWn(raw json.RawMessage) (wn []float64, err error) {
	var scalar float64
	if err := json.Unmarshal(raw, &scalar); err == nil {
		return []float64{scalar}, nil
	}
	var pair []float64
	if err := json.Unmarshal(raw, &pair); err == nil {
		return pair, nil
	}
	return nil, fmt.Errorf("wn is neither a float nor a [low, high] pair: %s", raw)
}

// bandFromVector builds the Band a golden vector describes
func bandFromVector(btype string, wn []float64) (band Band, err error) {
	switch btype {
	case "lowpass":
		band = Lowpass(wn[0])
	case "highpass":
		band = Highpass(wn[0])
	case "bandpass":
		band = Bandpass(wn[0], wn[1])
	case "bandstop":
		band = Bandstop(wn[0], wn[1])
	default:
		err = fmt.Errorf("unknown btype %q", btype)
	}
	return
}

// butterTolerance was measured against SciPy (see testdata/generate.py):
// most cases agree within a few ULPs (~1e-16), worst is bandstop 0.01-0.99
// at ~7e-15, where prewarp near Nyquist makes the pole quadratic ill-conditioned.
// Set ~10x above the worst for headroom across platforms.
var butterTolerance = tolerance{rel: 1e-13, abs: 1e-13}

func TestButterGoldenVectors(t *testing.T) {
	runGolden(t, "butter", func(t *testing.T, v butterVector) {
		wn, err := parseWn(v.Params.Wn)
		require.NoError(t, err)
		band, err := bandFromVector(v.Params.BType, wn)
		require.NoError(t, err)

		gotB, gotA, err := Butter(v.Params.Order, band)
		require.NoError(t, err, v.Description)

		assertAllClose(t, "b", gotB, v.Output.B, butterTolerance)
		assertAllClose(t, "a", gotA, v.Output.A, butterTolerance)
	})
}

func TestButterRejectsInvalidInput(t *testing.T) {
	cases := []struct {
		name    string
		order   int
		band    Band
		wantErr string
	}{
		{"order 0", 0, Lowpass(0.3), "order must be >= 1"},
		{"cutoff 0", 2, Lowpass(0), "Lowpass(0): frequencies must be in (0, 1)"},
		{"cutoff 1", 2, Highpass(1), "Highpass(1): frequencies must be in (0, 1)"},
		{"cutoff NaN", 2, Lowpass(math.NaN()), "frequencies must be in (0, 1)"},
		{"high edge 1", 2, Bandstop(0.2, 1), "Bandstop(0.2, 1): frequencies must be in (0, 1)"},
		{"low >= high", 2, Bandpass(0.5, 0.2), "Bandpass(0.5, 0.2): low edge must be below high edge"},
		{"zero Band", 2, Band{}, "zero Band"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, _, err := Butter(c.order, c.band)
			assert.ErrorContains(t, err, c.wantErr)
		})
	}
}
