# astrovodka.datacube — DataCube

Defines the output data cube grid: pixels, channels, and WCS.

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
| `spectral_unit` | `str` | **`'MHz'`** (default), `'km/s'`, or `None` (auto) |

> **Notes:**
> - All parameters keyword-only.
> - `spectral_unit='MHz'` converts velocity to frequency. WCS CTYPE3=`FREQ`.
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
