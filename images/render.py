"""Render the SVGs in images/ to PNGs in images/png/, which git ignores.

    uv run render.py                      # every image at its default width
    uv run render.py social               # one image
    uv run render.py icon-mark --width 64 # one image at another width

The height follows from each SVG's proportions.
"""

import argparse
from pathlib import Path

import resvg_py

HERE = Path(__file__).parent
OUT = HERE / "png"

# Default widths in pixels: the size each image is uploaded at.
WIDTHS = {
    "logo": 1100,
    "icon-mark": 512,
    "icon-word": 512,
    "social": 1280,
    "social-square": 1024,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "names",
        nargs="*",
        choices=[[], *WIDTHS],
        help="images to render, all by default",
    )
    parser.add_argument(
        "--width", type=int, help="width in pixels, each image's default if not given"
    )
    args = parser.parse_args()

    OUT.mkdir(exist_ok=True)
    for name in args.names or WIDTHS:
        width = args.width or WIDTHS[name]
        png = bytes(
            resvg_py.svg_to_bytes(svg_path=str(HERE / f"{name}.svg"), width=width)
        )
        path = OUT / f"{name}-{width}.png"
        path.write_bytes(png)
        print(f"wrote images/png/{path.name}, {len(png) // 1024} KB")


if __name__ == "__main__":
    main()
