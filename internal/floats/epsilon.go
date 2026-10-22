package floats

// Epsilon is the gap between 1 and the next float64 (2^-52):
// the relative size of one rounding, and so the unit for expressing relative error in ULPs.
// It can undercount by up to 2x, since the relative gap between neighbouring float64s
// varies between 2^-53 and 2^-52 within each power of two.
const Epsilon = 0x1p-52
