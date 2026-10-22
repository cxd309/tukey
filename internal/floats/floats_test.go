package floats

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestSum(t *testing.T) {
	cases := []struct {
		name   string
		values []float64
		want   float64
	}{
		{"empty", nil, 0},
		{"exact", []float64{1, 2, 3}, 6},
		// a plain sum gives 0, as each 1 is lost against 1e100
		{"small lost against large", []float64{1, 1e100, 1, -1e100}, 2},
		// 0.1 isn't exact in binary, so a plain sum of ten drifts to 0.9999999999999999
		{"many small", []float64{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, 1},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			var s Sum
			for _, v := range c.values {
				s.Add(v)
			}
			assert.Equal(t, c.want, s.Total())
		})
	}
}
