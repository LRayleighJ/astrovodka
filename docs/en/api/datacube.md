# astrovodka.datacube — DataCube

Defines the output data cube grid: pixels, channels, and WCS.

This class subclasses `martini.datacube.DataCube`. Everything except `spectral_unit` comes from MARTINI.

## Module Constants

- **`HIfreq`** = `1.420405751e9 * U.Hz` — HI 21-cm rest frequency.

## DataCube

### `__init__(self, *, n_px_x, n_px_y, n_channels, px_size, channel_width, spectral_centre=0 km/s, ra=0 deg, dec=0 deg, stokes_axis=False, coordinate_frame=ICRS(), specsys='icrs', spectral_unit='MHz')`

| Name | Type | Description |
|------|------|-------------|
| `n_px_x` | `int` | Pixel count along RA. Even integers preferred |
| `n_px_y` | `int` | Pixel count along Dec |
| `n_channels` | `int` | Number of spectral channels |
| `px_size` | `Quantity[arcsec]` | Angular pixel size |
| `channel_width` | `Quantity[km/s \| Hz]` | Channel width |
| `spectral_centre` | `Quantity` | Centre of spectral axis (default: 0 km/s) |
| `ra` | `Quantity[deg]` | RA of cube centroid (default: 0) |
| `dec` | `Quantity[deg]` | Dec of cube centroid (default: 0) |
| `stokes_axis` | `bool` | Include Stokes axis |
| `coordinate_frame` | `BaseRADecFrame` | Frame (default: ICRS) |
| `specsys` | `str` | Spectral ref: `'icrs'`, `'gcrs'`, `'lsrk'` |
| `spectral_unit` | `str \| None` | `None` (default) adopts the units of `channel_width`; or force `'MHz'` / `'km/s'` |

> **Notes:**
> - All parameters keyword-only.
> - `spectral_unit=None` (default) reproduces MARTINI: a velocity `channel_width` gives a
>   velocity axis (WCS CTYPE3=`VRAD`), a frequency one gives `FREQ`.
> - `spectral_unit='MHz'` converts a velocity `channel_width` to frequency. The conversion
>   is *differential* (`Δf = f_rest·Δv/c`, i.e. 20 km/s → 94.76 kHz near the HI line), not
>   the absolute frequency at that velocity.
> - `channel_width` must have frequency or velocity units, else `ValueError`.
> - The resolved unit is stored as `_spectral_unit` and is reused by `copy()`,
>   `save_state()`/`load_state()`, `from_wcs()` and `Martini.reset()`, so a cube's spectral
>   axis never silently changes type.
> - Internal `_array`: shape `(n_px_x, n_px_y, n_channels)`, unit `Jy/pix²`.

---

### `from_wcs(cls, input_wcs, specsys=None)`

| Name | Type | Description |
|------|------|-------------|
| `input_wcs` | `astropy.wcs.WCS` | Existing WCS (e.g., from FITS header) |
| `specsys` | `str \| None` | Spectral ref override |

**Returns:** `DataCube`

> Extracts all dimensions from WCS. Requires square pixels. Warns if RESTFRQ missing.

---

## Properties

| Property | Returns | Description |
|----------|---------|-------------|
| `wcs` | `WCS` | Lazy-built WCS |
| `channel_mids` | `SpectralCoord` | Channel centres (native unit) |
| `channel_edges` | `SpectralCoord` | Channel edges (n_channels+1 values) |
| `velocity_channel_mids` | `Quantity[km/s]` | Centres in km/s |
| `velocity_channel_edges` | `Quantity[km/s]` | Edges in km/s |
| `frequency_channel_mids` | `Quantity[Hz]` | Centres in Hz |
| `frequency_channel_edges` | `Quantity[Hz]` | Edges in Hz |
| `spatial_slices` | `Iterator[Quantity]` | Iterate (ny,nx) slices |
| `spectra` | `Iterator[Quantity]` | Iterate (nchan,) per pixel |

## Methods

### `save_state(self, filename, overwrite=False)`

Save array + WCS to HDF5. Requires `h5py`.

### `load_state(cls, filename)`

Restore from `save_state()`. Returns `DataCube`.

### `copy(self)` → `DataCube`

Deep copy.

### `add_pad(self, pad: tuple[int,int])` / `drop_pad(self)`

Padding management for beam convolution. Called internally.
