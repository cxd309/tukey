<a href="../README.md"><img src="logo.svg" alt="tukey" width="120" align="right"></a>

# Images

tukey's logo, icons and social images, built by `build.py`.

| File                | Size     | For                                             |
| ------------------- | -------- | ----------------------------------------------- |
| `logo.svg`          | 220x96   | README, transparent so works with light or dark |
| `icon-mark.svg`     | 512      | small icons: under ~128 px                      |
| `icon-word.svg`     | 512      | large icons: above ~128 px                      |
| `social.svg`        | 1280x640 | GitHub social preview and link cards            |
| `social-square.svg` | 1024     | square social and thumbnails                    |

## The design

The shape is a Tukey window (named after John Tukey) with an alpha of 0.5.
The word is cut out, so it will match the background theme.

The word uses [Ubuntu Regular](https://design.ubuntu.com/font) converted to outlines,
so a copy of the font is not required to render the SVG.

The highlight colour is the Go blue: `#00ADD8`.
Otherwise use GitHub dark theme colours.

`icon-word`, `social` and `social-square` all use 'squircle' rounded corners.

## Rebuilding

Requires uv.

```sh
just images
```

On first run, downloads Ubuntu Regular from [Google Fonts repository](https://github.com/google/fonts/tree/main/ufl/ubuntu) into `fonts/`.
This is gitignored but is checked against a SHA-256 so if the font is updated then the logo will not change unexpectedly.

## Rendering PNGs

```sh
just images-png                         # every image at its default width
just images-png social                  # one image
just images-png icon-mark --width 64    # one image at another width
```

PNGs go in `png/`, which is gitignored, named with their width, e.g. `social-1280.png`.
The height follows from each image's proportions.

## Licences

Images are covered by the same licence as the rest of the repo: [LICENSE](../LICENSE).

Ubuntu font is not redistributed here, usage follows the [Ubuntu Font Licence](https://ubuntu.com/legal/font-licence).
