"""Build tukey's logo, icons and social images as SVGs.

The shape is a Tukey window with alpha 0.5.
The word is set in Ubuntu Regular and converted to outlines,
so the SVGs don't require access to the font to display.

    uv run build.py
"""

import hashlib
import math
import urllib.request
from pathlib import Path

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

HERE = Path(__file__).parent

FONT_URL = "https://github.com/google/fonts/raw/main/ufl/ubuntu/Ubuntu-Regular.ttf"
FONT_SHA256 = "3128df86a31805618436d0ae5651ba4285d0c9de0a39057d025f64ee33bceb64"
FONT_PATH = HERE / "fonts" / "Ubuntu-Regular.ttf"

BLUE = "#00ADD8"  # Go blue
DARK = "#0D1117"  # GitHub dark background
TEXT = "#E6EDF3"  # GitHub dark text
MUTED = "#8B949E"  # GitHub dark secondary text

# The logo's window: 220 wide and 96 tall, with the word centred inside.
LOGO_X, LOGO_TOP, LOGO_WIDTH, LOGO_HEIGHT = 10.0, 14.0, 220.0, 96.0
LOGO_WORD_SIZE, LOGO_WORD_BASELINE = 55.0, 96.0

TAGLINE = "Signal processing for Go"
PACKAGES = "fft    window    dsp"


def main():
    font = Font(load_font())
    logo = logo_path(font)

    write("logo.svg", logo_svg(logo))
    write(
        "icon-word.svg",
        square(512, place(logo, 512, 420, (512 - 420 * LOGO_HEIGHT / LOGO_WIDTH) / 2)),
    )
    write(
        "icon-mark.svg",
        square(
            512,
            place(
                window_path(closed=True),
                512,
                464,
                (512 - 464 * LOGO_HEIGHT / LOGO_WIDTH) / 2,
            ),
            rounded=False,
        ),
    )
    write(
        "social.svg",
        canvas(
            1280,
            640,
            place(logo, 1280, 600, 110)
            + font.text(TAGLINE, 52, 640, 460, TEXT)
            + font.text(PACKAGES, 38, 640, 545, MUTED),
        ),
    )
    write(
        "social-square.svg",
        canvas(
            1024,
            1024,
            place(logo, 1024, 760, 230)
            + font.text(TAGLINE, 66, 512, 690, TEXT)
            + font.text(PACKAGES, 48, 512, 790, MUTED),
        ),
    )


def load_font():
    """
    Font is gitignored, so on first run download a local copy.
    This will check if the font file exists at the path and if not then downloads it.
    Compared with sha256 to check if it has changed since last download so logo doesn't change unexpectedly.
    """
    if not FONT_PATH.exists():
        print(f"downloading {FONT_URL}")
        FONT_PATH.parent.mkdir(exist_ok=True)
        with urllib.request.urlopen(FONT_URL, timeout=60) as response:
            FONT_PATH.write_bytes(response.read())
    digest = hashlib.sha256(FONT_PATH.read_bytes()).hexdigest()
    if digest != FONT_SHA256:
        raise SystemExit(
            f"{FONT_PATH} has SHA-256 {digest}, want {FONT_SHA256}. Delete it to download it again."
        )
    return FONT_PATH


class Font:
    """Ubuntu Regular, for turning text into outline paths."""

    def __init__(self, path):
        self.font = TTFont(path)
        self.glyphs = self.font.getGlyphSet()
        self.cmap = self.font.getBestCmap()
        self.upem = self.font["head"].unitsPerEm
        self.kern = {}
        for table in self.font["kern"].kernTables:
            self.kern.update(table.kernTable)

    def outline(self, text, size, centre_x, baseline):
        """Returns text's outline as path data, centred on centre_x."""
        scale = size / self.upem
        names = [self.cmap[ord(c)] for c in text]
        steps = [
            self.font["hmtx"][name][0] + self.kern.get((name, after), 0)
            for name, after in zip(names, [*names[1:], None], strict=True)
        ]
        x = centre_x - sum(steps) * scale / 2
        pen = SVGPathPen(self.glyphs)
        for name, step in zip(names, steps, strict=True):
            self.glyphs[name].draw(
                TransformPen(pen, (scale, 0, 0, -scale, x, baseline))
            )
            x += step * scale
        return pen.getCommands()

    def text(self, text, size, centre_x, baseline, fill):
        return (
            f'<path fill="{fill}" d="{self.outline(text, size, centre_x, baseline)}"/>'
        )


def tukey(t, alpha=0.5):
    """A Tukey window at t from 0 to 1: cosine tapers over alpha/2 at each end, flat between."""
    edge = min(t, 1 - t)
    if edge < alpha / 2:
        return 0.5 * (1 - math.cos(2 * math.pi * edge / alpha))
    return 1.0


def window_path(closed, points=240):
    """Returns the window's outline as path data, in logo units."""
    base = LOGO_TOP + LOGO_HEIGHT
    commands = []
    for i in range(points + 1):
        t = i / points
        command = "L" if i else "M"
        commands.append(
            f"{command}{LOGO_X + t * LOGO_WIDTH:.2f} {base - tukey(t) * LOGO_HEIGHT:.2f}"
        )
    return " ".join(commands) + (" Z" if closed else "")


def logo_path(font):
    """Returns the window with the word as holes, as path data for an even-odd fill."""
    word = font.outline(
        "tukey", LOGO_WORD_SIZE, LOGO_X + LOGO_WIDTH / 2, LOGO_WORD_BASELINE
    )
    return f"{window_path(closed=True)} {word}"


def logo_svg(path):
    """The logo alone, transparent, so the page shows through the word in light and dark themes."""
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{LOGO_X:g} {LOGO_TOP:g} {LOGO_WIDTH:g} {LOGO_HEIGHT:g}" role="img" aria-label="tukey">\n'
        f'  <path fill="{BLUE}" fill-rule="evenodd" d="{path}"/>\n'
        "</svg>\n"
    )


def place(path, canvas_width, width, top):
    """Returns path, in logo units, scaled to width and centred across the canvas with its top at top."""
    scale = width / LOGO_WIDTH
    x = (canvas_width - width) / 2 - LOGO_X * scale
    y = top - LOGO_TOP * scale
    return f'<path fill="{BLUE}" fill-rule="evenodd" transform="translate({x:.2f} {y:.2f}) scale({scale:.4f})" d="{path}"/>'


def canvas(width, height, body, rounded=True):
    """An SVG of the given size on GitHub's dark background, with squircle corners if rounded."""
    if rounded:
        background = f'<path fill="{DARK}" d="{squircle(width, height)}"/>'
    else:
        background = f'<rect width="{width}" height="{height}" fill="{DARK}"/>'
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n'
        f"  {background}\n"
        f"  {body}\n"
        "</svg>\n"
    )


def square(size, body, rounded=True):
    return canvas(size, size, body, rounded)


def squircle(width, height, radius=0.2237, exponent=5, points=24):
    """
    Returns a rectangle with superellipse corners as path data.
    Each corner is a quarter of |x|^n + |y|^n = 1,
    so its curvature eases in rather than starting abruptly.
    radius is a fraction of the shorter side.
    """
    r = radius * min(width, height)
    centres = [
        (width - r, r, -90),
        (width - r, height - r, 0),
        (r, height - r, 90),
        (r, r, 180),
    ]
    commands = []
    for cx, cy, start in centres:
        for i in range(points + 1):
            theta = math.radians(start + 90 * i / points)
            c, s = math.cos(theta), math.sin(theta)
            x = cx + r * math.copysign(abs(c) ** (2 / exponent), c)
            y = cy + r * math.copysign(abs(s) ** (2 / exponent), s)
            commands.append(f"{'L' if commands else 'M'}{x:.2f} {y:.2f}")
    return " ".join(commands) + " Z"


def write(name, svg):
    (HERE / name).write_text(svg, encoding="ascii")
    print(f"wrote images/{name}")


if __name__ == "__main__":
    main()
