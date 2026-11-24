"""Reference vectors for tukey/waveform: scipy.signal's waveforms"""

from waveform.chirp import generate_chirp
from waveform.square import generate_square
from waveform.sweep_poly import generate_sweep_poly

# every waveform generator, each writes 1 or more categories
GENERATORS = [generate_chirp, generate_sweep_poly, generate_square]
