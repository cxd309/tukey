"""Reference vectors for tukey/waveform: scipy.signal's waveforms"""

from waveform.chirp import generate_chirp

# every waveform generator, each writes 1 or more categories
GENERATORS = [generate_chirp]
