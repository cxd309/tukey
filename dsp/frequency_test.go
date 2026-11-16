package dsp

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// frequencyUnits are every unit, with Nyquist in it at fs = 1000
var frequencyUnits = []struct {
	name    string
	unit    FrequencyUnit
	nyquist float64
}{
	{"Nyquist", UnitNyquist, 1},
	{"rad/sample", UnitRadPerSample, math.Pi},
	{"cycles/sample", UnitCyclesPerSample, 0.5},
	{"Hz", UnitHz, 500},
}

// Nyquist in one unit is Nyquist in every other, and 0 is 0.
func TestConvertFrequencyNyquist(t *testing.T) {
	for _, from := range frequencyUnits {
		for _, to := range frequencyUnits {
			g, err := ConvertFrequency(from.nyquist, from.unit, to.unit, 1000)
			require.NoError(t, err)
			assert.InEpsilon(t, to.nyquist, g, 1e-15, "%s to %s", from.name, to.name)

			g, err = ConvertFrequency(0.0, from.unit, to.unit, 1000)
			require.NoError(t, err)
			assert.Equal(t, 0.0, g, "%s to %s", from.name, to.name)
		}
	}
}

// Without π, a conversion is one division by a power of 2 times fs,
// so it's exactly f/(fs/2) or f/fs.
func TestConvertFrequencyIsCorrectlyRounded(t *testing.T) {
	const fs = 44100.0
	for _, f := range []float64{50, 1000.0 / 3, 12345.678, 22049.9} {
		wn, err := HzToNyquist(f, fs)
		require.NoError(t, err)
		assert.Equal(t, f/(fs/2), wn, "f %g", f)

		cycles, err := ConvertFrequency(f, UnitHz, UnitCyclesPerSample, fs)
		require.NoError(t, err)
		assert.Equal(t, f/fs, cycles, "f %g", f)

		back, err := ConvertFrequency(cycles, UnitCyclesPerSample, UnitNyquist, fs)
		require.NoError(t, err)
		assert.Equal(t, 2*cycles, back, "f %g", f)
	}
}

func TestConvertFrequencySlice(t *testing.T) {
	w := []float64{0, math.Pi / 4, math.Pi / 2}
	f, err := RadPerSampleToHz(w, 1000)
	require.NoError(t, err)
	assert.InDeltaSlice(t, []float64{0, 125, 250}, f, 1e-12)
	assert.Equal(t, []float64{0, math.Pi / 4, math.Pi / 2}, w, "w is unchanged")

	for _, x := range [][]float64{nil, {}} {
		g, err := HzToNyquist(x, 1000)
		require.NoError(t, err)
		assert.NotNil(t, g)
		assert.Empty(t, g)
	}
}

// Each direct function is ConvertFrequency between its two units.
func TestFrequencyShortcuts(t *testing.T) {
	const fs = 48000.0
	x := []float64{0, 0.1, 0.5, 1, 3, 1000}
	convert := func(from, to FrequencyUnit) []float64 {
		g, err := ConvertFrequency(x, from, to, fs)
		require.NoError(t, err)
		return g
	}
	withFs := func(g []float64, err error) []float64 {
		require.NoError(t, err)
		return g
	}
	assert.Equal(t, convert(UnitHz, UnitNyquist), withFs(HzToNyquist(x, fs)))
	assert.Equal(t, convert(UnitNyquist, UnitHz), withFs(NyquistToHz(x, fs)))
	assert.Equal(t, convert(UnitHz, UnitRadPerSample), withFs(HzToRadPerSample(x, fs)))
	assert.Equal(t, convert(UnitRadPerSample, UnitHz), withFs(RadPerSampleToHz(x, fs)))
	assert.Equal(t, convert(UnitNyquist, UnitRadPerSample), NyquistToRadPerSample(x))
	assert.Equal(t, convert(UnitRadPerSample, UnitNyquist), RadPerSampleToNyquist(x))
}

// A single value gives a single value.
func TestFrequencyShortcutsTakeOneValue(t *testing.T) {
	wn, err := HzToNyquist(50.0, 1000)
	require.NoError(t, err)
	assert.Equal(t, 0.1, wn)
	assert.Equal(t, math.Pi/2, NyquistToRadPerSample(0.5))
}

func TestConvertFrequencyRejectsInvalidInput(t *testing.T) {
	for _, fs := range []float64{0, -1, math.NaN(), math.Inf(1)} {
		_, err := ConvertFrequency(1.0, UnitHz, UnitNyquist, fs)
		assert.ErrorIs(t, err, ErrInvalidParameter, "from Hz, fs %g", fs)
		_, err = ConvertFrequency(1.0, UnitNyquist, UnitHz, fs)
		assert.ErrorIs(t, err, ErrInvalidParameter, "to Hz, fs %g", fs)
	}
	// fs is only used for Hz
	_, err := ConvertFrequency(1.0, UnitNyquist, UnitRadPerSample, math.NaN())
	assert.NoError(t, err)

	for _, unit := range []FrequencyUnit{-1, UnitHz + 1} {
		_, err = ConvertFrequency(1.0, unit, UnitNyquist, 1000)
		assert.ErrorIs(t, err, ErrInvalidParameter, "from unit %d", unit)
		_, err = ConvertFrequency(1.0, UnitNyquist, unit, 1000)
		assert.ErrorIs(t, err, ErrInvalidParameter, "to unit %d", unit)
	}
}
