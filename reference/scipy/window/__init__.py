"""Reference vectors for tukey/window: scipy.signal.windows"""

from window.windows import generate_windows

# every window generator, each writes 1 or more categories
GENERATORS = [generate_windows]
