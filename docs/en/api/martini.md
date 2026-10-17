# astrovodka.martini — Core Pipeline

The `Martini` class orchestrates the mock-observing pipeline by connecting source, DataCube, beam, noise, spectral model, and SPH kernel.

These classes subclass `martini.martini.Martini` and `martini.martini.GlobalProfile`. The pipeline itself is inherited unchanged.

## Martini

### `__init__(self, *, source, datacube, beam=None, noise=None, sph_kernel, spectral_model, quiet=False)`

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `source` | `SPHSource` | HI-emitting source object (required) |
| `datacube` | `DataCube` | Data cube grid definition (required) |
| `beam` | `_BaseBeam \| None` | Telescope beam. `None` = skip convolution (default: None) |
| `noise` | `_BaseNoise \| None` | Noise model. `None` = skip noise (default: None) |
| `sph_kernel` | `_BaseSPHKernel` | SPH smoothing kernel (required) |
| `spectral_model` | `_BaseSpectrum` | 21-cm spectral line model (required) |
| `quiet` | `bool` | Suppress stdout (default: False) |

**Returns:** None

> **Notes:**
> - All parameters are keyword-only.
> - During init, particles are pruned (spatially, spectrally, by mass) to skip non-contributing ones.
> - `DataCube._array` evolves: `Jy/pix²` after insertion → `Jy/beam` after convolution.

---

### `insert_source_in_cube(self, skip_validation=False, progressbar=None, ncpu=1)`

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `skip_validation` | `bool` | Skip SPH kernel accuracy checks (default: False) |
| `progressbar` | `bool \| None` | Show tqdm bar. `None` = auto |
| `ncpu` | `int` | Parallel processes. Requires `multiprocess` (default: 1) |

**Returns:** None

> **Notes:**
> - Most expensive step — each pixel loops over nearby particles.
> - Parallelism uses `multiprocess` (dill), NOT `multiprocessing` (pickle).
> - Per-pixel parallelism: each worker processes one pixel independently.
> - Mass conservation reported: flux-summed mass should match input HI mass within ~1%.

---

### `convolve_beam(self)`

**Returns:** None

> **Notes:**
> - Requires `beam` to be set. Warns and skips if None.
> - Uses `scipy.signal.fftconvolve` for FFT-based convolution.
> - After convolution, padding removed; array units change to `Jy/beam`.

---

### `add_noise(self)`

**Returns:** None

> **Notes:**
> - Requires both `beam` and `noise` set. Warns and skips otherwise.
> - Noise added **before** convolution. RMS auto-scaled: `rms_pre = rms_post × 2.19568 × √(π × σ_maj × σ_min)`.

---

### `write_fits(self, filename, overwrite=True, obj_name='MOCK', channels=None, dtype=None)`

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `filename` | `str` | Output path. `.fits` appended if not present |
| `overwrite` | `bool` | Allow overwrite (default: True) |
| `obj_name` | `str` | OBJECT header (max 16 chars) |
| `dtype` | `str \| np.dtype \| None` | Output type. `None` = keep native |

**Returns:** None

> **Notes:** Instrument header = `ASTROVODKA`. Axes transposed: `(RA,Dec,Chan)` → FITS `(Chan,Dec,RA)`.

---

### `write_hdf5(self, filename, overwrite=True, memmap=False, compact=False, channels=None)`

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `memmap` | `bool` | Create in-memory file, return handle |
| `compact` | `bool` | Omit coordinate arrays to save space |

**Returns:** `h5py.File | None`

---

### `write_beam_fits(self, filename, overwrite=True, channels=None)`

**Returns:** None. Raises `ValueError` if beam is None.

---

### `preview(self, max_points=5000, fig=1, lim=None, vlim=None, point_scaling='auto', title='', save=None)`

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `max_points` | `int` | Max particles per panel (random subsample) |
| `lim` | `str \| Quantity[deg] \| None` | Spatial limit. `'datacube'` = cube bounds |
| `vlim` | `str \| Quantity[km/s] \| None` | Velocity limit |
| `point_scaling` | `str` | `'auto'` (mass+smoothing) or `'fixed'` |
| `save` | `str \| None` | Save figure to file |

**Returns:** `matplotlib.figure.Figure`

> 3-panel figure: sky projection + 2 PV diagrams. Red boxes = DataCube extent.

---

### `reset(self)`

**Returns:** None. Re-initializes DataCube with zeros, keeping source + sub-module configs.

---

## GlobalProfile

Simplified class for spatially-integrated spectra only. No beam, noise, or SPH kernel needed.

### `__init__(self, *, source, spectral_model, n_channels, channel_width, spectral_centre=0 km/s, quiet=False)`

| Name | Type | Description |
|------|------|-------------|
| `source` | `SPHSource` | HI source |
| `spectral_model` | `_BaseSpectrum` | Spectral line model |
| `n_channels` | `int` | Number of channels |
| `channel_width` | `Quantity[km/s \| Hz]` | Channel width |
| `spectral_centre` | `Quantity` | Spectral centre (default: 0 km/s) |
| `quiet` | `bool` | Suppress output |

> Uses 1×1 pixel DataCube + DiracDeltaKernel. All particles contribute.

### `spectrum` (property)

**Returns:** `Quantity[Jy]` — 1D flux density array. Lazy-evaluated.

### `plot_spectrum(self, fig=1, title='', channels='velocity', show_vsys=True, save=None)`

| Name | Type | Description |
|------|------|-------------|
| `channels` | `str` | `'velocity'` or `'frequency'` |
| `show_vsys` | `bool` | Draw line at systemic velocity |
| `save` | `str \| None` | Save figure to file |

**Returns:** `Figure`
