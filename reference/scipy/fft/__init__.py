"""
Reference vectors for tukey/fft:
    - transforms
    - options
    - frequency helpers
    - shift helpers
"""

from fft.options import generate_fft_options
from fft.transforms import generate_fft

# every fft generator, each write at least 1 category
GENERATORS = [generate_fft, generate_fft_options]
