# astrovodka.noise — Noise Models

Generate noise cubes. Built-in: `GaussianNoise`.

## GaussianNoise

### `__init__(self, rms=1 Jy/beam, seed=0)`

| Name | Type | Description |
|------|------|-------------|
| `rms` | `Quantity[Jy/beam]` | Post-convolution RMS (default: 1) |
| `seed` | `int` | Random seed (default: 0) |

### `generate(self, datacube, beam)` → `Quantity[Jy/arcsec²]`

| Name | Type | Description |
|------|------|-------------|
| `datacube` | `DataCube` | Array shape + pixel scale |
| `beam` | `_BaseBeam` | FWHM for RMS scaling |

> Noise is added **before** convolution. RMS scaled: `rms_pre = rms_post × 2.19568 × √(π × σ_maj_px × σ_min_px)`.

### `reset_rng(self)`

Re-seeds with original seed for reproducible realizations.

---

> **Custom noise:** Inherit from `_BaseNoise`, implement `generate(datacube, beam)`.
