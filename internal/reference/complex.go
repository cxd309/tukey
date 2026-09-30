package reference

import "fmt"

// Complex is a vector of complex values as stored in reference vector JSON
// {"re":[...], "im":[...]}
type Complex struct {
	Re []float64 `json:"re"`
	Im []float64 `json:"im"`
}

// Values combines Re and Im into complex values
func (c Complex) Values() (values []complex128, err error) {
	if len(c.Re) != len(c.Im) {
		return nil, fmt.Errorf("re has %d values but im has %d", len(c.Re), len(c.Im))
	}
	values = make([]complex128, len(c.Re))
	for i := range c.Re {
		values[i] = complex(c.Re[i], c.Im[i])
	}
	return
}
