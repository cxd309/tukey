package reference

import (
	"encoding/json"
	"math"
)

// Floats is a vector of floats as stored in reference vector JSON,
// where NaN, which JSON can't represent, is null
type Floats []float64

// UnmarshalJSON decodes a JSON array of numbers and nulls, nulls as NaN
func (f *Floats) UnmarshalJSON(data []byte) (err error) {
	var values []*float64
	if err = json.Unmarshal(data, &values); err != nil {
		return
	}
	*f = make(Floats, len(values))
	for i, v := range values {
		if v == nil {
			(*f)[i] = math.NaN()
			continue
		}
		(*f)[i] = *v
	}
	return
}
