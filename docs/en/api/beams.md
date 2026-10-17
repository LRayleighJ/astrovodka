# astrovodka.beams — Telescope Beam Models

Models the telescope beam. Built-in: `GaussianBeam` (elliptical 2D Gaussian).

The implementation is MARTINI's: it is imported from `martini.beams` and re-exported here.

## GaussianBeam

### `__init__(self, bmaj=15 arcsec, bmin=None, bpa=0 deg, truncate=4.0)`

| Name | Type | Description |
|------|------|-------------|
| `bmaj` | `Quantity[arcsec]` | Major axis FWHM (default: 15") |
| `bmin` | `Quantity[arcsec] \| None` | Minor axis FWHM. `None` = circular |
| `bpa` | `Quantity[deg]` | Position angle, E of N (default: 0°) |
| `truncate` | `float` | Truncation radius in FWHM (default: 4.0) |

> Area = π × bmaj × bmin / (4 × ln2). The beam is integrated over each pixel using RectBivariateSpline.

### `needs_pad(self)` → `tuple[int, int]`

Returns (pad_x, pad_y) = half kernel size. Used before convolution.

---

> **Custom beams:** Inherit from `_BaseBeam`, implement `f_kernel()` (returns callable), `kernel_size_px()` (returns half-size), and optionally `init_beam_header()`.
