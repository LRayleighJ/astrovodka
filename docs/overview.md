# astrovodka — Synthetic HI Observations with Cosmological Distances & K13 HI Modeling

**astrovodka** creates mock resolved HI (21 cm) spectral line observations (data cubes)
from IllustrisTNG SPH simulations. It improves upon [MARTINI](https://github.com/kyleaoman/martini)
in two key areas:

1. **Proper cosmological distances** via `astropy.cosmology.Planck15` (replaces the simple linear Hubble law)
2. **K13 volumetric HI/H2 model** (Krumholz 2013, Diemer+2018 implementation) with adaptive per-galaxy FFT UV field computation

---

## Installation

```bash
pip install -e /path/to/astrovodka.0.1/

# Optional: parallel processing, HDF5 output, plotting
pip install -e "/path/to/astrovodka.0.1/[parallel,hdf5_output,plotting]"
```

**Requirements:** Python ≥ 3.10, numpy, scipy, astropy, tqdm

---

## Quick Start

```python
import numpy as np
import astropy.units as U
from astrovodka import Martini, DataCube
from astrovodka.sources import SPHSource
from astrovodka.beams import GaussianBeam
from astrovodka.spectral_models import GaussianSpectrum
from astrovodka.sph_kernels import CubicSplineKernel

# 1. Create a simple rotating disk source at z=0.01
N = 1000
r = np.random.rand(N) * 10 * U.kpc
phi = np.random.rand(N) * 2 * np.pi
xyz_g = U.Quantity([np.zeros(N)*U.kpc, r*np.cos(phi), r*np.sin(phi)])
vxyz_g = U.Quantity([np.zeros(N)*U.km/U.s, -100*np.sin(phi)*U.km/U.s, 100*np.cos(phi)*U.km/U.s])

source = SPHSource(
    redshift=0.01,               # ★ NEW: cosmological redshift mode
    ra=30*U.deg, dec=45*U.deg,
    T_g=np.ones(N)*1e4*U.K,
    mHI_g=np.ones(N)*1e8*U.Msun,
    xyz_g=xyz_g, vxyz_g=vxyz_g,
    hsm_g=np.ones(N)*U.kpc,
)

# 2. Set up the data cube (spectral axis defaults to MHz)
datacube = DataCube(
    n_px_x=128, n_px_y=128, n_channels=64,
    px_size=2*U.arcsec,
    channel_width=20*U.km/U.s,         # can be velocity...
    spectral_centre=source.vsys,        # systemic velocity
    spectral_unit="MHz",                # ...but output axis is MHz (default)
)

# 3. Configure sub-modules
beam = GaussianBeam(bmaj=8*U.arcsec, bmin=6*U.arcsec, bpa=30*U.deg)
spectral_model = GaussianSpectrum(sigma="thermal")  # temperature-dependent width
sph_kernel = CubicSplineKernel()                     # adaptive fallback included

# 4. Run the pipeline
M = Martini(source=source, datacube=datacube, beam=beam,
            spectral_model=spectral_model, sph_kernel=sph_kernel)
M.insert_source_in_cube()   # particles → pixels (can use ncpu=N for parallel)
M.convolve_beam()           # convolve with telescope beam
M.write_fits("output.fits") # save to FITS

# 5. Access channel information
print(datacube.frequency_channel_mids)  # in Hz
print(datacube.velocity_channel_edges)  # in km/s
```

### TNG Quick Start

```python
from astrovodka.sources import TNGSource

source = TNGSource(
    basePath="/data/.../TNG/TNG-100/",  # path to TNG data
    snapNum=99,                          # snapshot number
    subID=400547,                        # subhalo ID
    observed_redshift=0.05,             # place at custom cosmological distance
    load_subhalo_only=True,             # fast: only subhalo particles
    uv_mode="global",                   # "global" or "adaptive_fft"
)

# source.D_L, source.D_A, source.m_H2_g are available
# source.mHI_g is the K13-computed HI mass per particle
print(f"D_L = {source.D_L:.1f}, D_A = {source.D_A:.1f}")
print(f"M_HI = {source.input_mass:.2e}, M_H2 = {source.m_H2_g.sum():.2e}")
```

---

## Module Reference

### 1. Core Pipeline — `Martini`

The `Martini` class orchestrates the mock-observing pipeline. It connects a source,
data cube, beam, noise model, spectral model, and SPH kernel, then executes the
pipeline steps in order.

| Method | Description |
|--------|-------------|
| `Martini(source, datacube, beam, noise, sph_kernel, spectral_model, quiet=False)` | Constructor. All parameters are keyword-only. `beam` and `noise` may be `None`. |
| `insert_source_in_cube(skip_validation=False, progressbar=None, ncpu=1)` | Populate the data cube with flux from source particles. Use `ncpu > 1` for parallel execution (requires `multiprocess`). Set `skip_validation=True` to bypass SPH kernel accuracy checks. |
| `convolve_beam()` | Convolve the data cube with the telescope beam. Requires `beam` to be set. |
| `add_noise()` | Add noise to the data cube. Requires both `beam` and `noise` to be set. |
| `write_fits(filename, overwrite=True, obj_name="MOCK")` | Write the data cube to a FITS file. Auto-appends `.fits` suffix. |
| `write_beam_fits(filename, overwrite=True)` | Write the beam image to a FITS file. |
| `write_hdf5(filename, overwrite=True, memmap=False, compact=False)` | Write cube + beam to HDF5 (requires `h5py`). `compact=True` omits coordinate arrays. |
| `preview(max_points=5000, fig=1, lim=None, vlim=None, point_scaling="auto", title="", save=None)` | Preview source particles overlaid with cube extent. Returns matplotlib Figure. |
| `reset()` | Re-initialize the data cube with zeros (keep source + sub-module config). |

#### `GlobalProfile` — Spectrum-Only Mode

A simplified version of `Martini` for producing only a spatially-integrated global spectrum.
No beam, noise, or SPH kernel required.

```python
from astrovodka import GlobalProfile

gp = GlobalProfile(
    source=source,
    spectral_model=GaussianSpectrum(sigma=7*U.km/U.s),
    n_channels=64,
    channel_width=10*U.km/U.s,
    spectral_centre=source.vsys,
)

# Spectrum is auto-evaluated on first access:
gp.spectrum        # U.Quantity[U.Jy] array
gp.plot_spectrum()  # quick matplotlib plot
```

| Method | Description |
|--------|-------------|
| `GlobalProfile(source, spectral_model, n_channels, channel_width, spectral_centre=0, quiet=False)` | Constructor. No beam, noise, or SPH kernel needed. |
| `insert_source_in_spectrum()` | Explicitly populate the spectrum. |
| `plot_spectrum(fig=1, title="", channels="velocity", show_vsys=True, save=None)` | Plot the spectrum as flux vs velocity or frequency. |
| `reset()` | Re-initialize the spectrum with zeros. |

---

### 2. DataCube — `DataCube`

Defines the output data cube grid: pixels, channels, and WCS.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_px_x`, `n_px_y` | int | — | Pixel count along RA and Dec axes |
| `n_channels` | int | — | Number of spectral channels |
| `px_size` | Quantity[arcsec] | — | Angular pixel size |
| `channel_width` | Quantity | — | Channel width (velocity or frequency) |
| `spectral_centre` | Quantity | 0 km/s | Centre of spectral axis |
| `ra`, `dec` | Quantity[deg] | 0 deg | Cube centroid sky coordinates |
| `stokes_axis` | bool | False | Include Stokes axis |
| `coordinate_frame` | BaseRADecFrame | ICRS() | Astropy coordinate frame |
| `specsys` | str | "icrs" | Spectral reference system |
| `spectral_unit` | str | **"MHz"** | Spectral axis unit: `"MHz"`, `"km/s"`, or `None` (auto) |

**Properties:**

| Property | Returns | Description |
|----------|---------|-------------|
| `wcs` | `astropy.wcs.WCS` | World Coordinate System (lazy-built) |
| `_array` | Quantity | Raw data cube array (shape: `(nx, ny, nchan)`) |
| `channel_mids` | SpectralCoord | Channel centres (native unit) |
| `channel_edges` | SpectralCoord | Channel edges (native unit) |
| `velocity_channel_mids` | Quantity[km/s] | Channel centres in km/s |
| `velocity_channel_edges` | Quantity[km/s] | Channel edges in km/s |
| `frequency_channel_mids` | Quantity[Hz] | Channel centres in Hz |
| `frequency_channel_edges` | Quantity[Hz] | Channel edges in Hz |
| `spatial_slices` | iterator | Iterate over (ny, nx) spatial slices |
| `spectra` | iterator | Iterate over (nchan,) spectra per pixel |

**Class Methods:**

| Method | Description |
|--------|-------------|
| `DataCube.from_wcs(input_wcs, specsys=None)` | Create DataCube from an existing astropy WCS (e.g., from a FITS header) |
| `DataCube.load_state(filename)` | Restore DataCube state from HDF5 |
| `save_state(filename, overwrite=False)` | Save DataCube state to HDF5 |

---

### 3. Sources

#### 3.1 `SPHSource` — Generic SPH Source

Abstracts an HI-emitting object from SPH particle arrays. **New:** supports cosmological
redshift mode.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `redshift` | float or None | None | ★ Cosmological redshift (use with `cosmology`) |
| `distance` | Quantity[Mpc] or None | None | Physical distance (legacy mode, use with `h`) |
| `cosmology` | Cosmology or None | None | Astropy cosmology (default: Planck15) |
| `vpeculiar` | Quantity[km/s] | 0 | Peculiar velocity |
| `rotation` | Rotation or None | None | scipy Rotation for source orientation |
| `L_coords` | L_coords or None | None | Angular momentum alignment tuple |
| `ra`, `dec` | Quantity[deg] | 0 | Source centroid |
| `h` | float | 0.7 | Hubble parameter (legacy mode only) |
| `T_g` | Quantity[K] | — | Gas temperature per particle |
| `mHI_g` | Quantity[Msun] | — | HI mass per particle |
| `xyz_g` | Quantity[kpc] | — | Particle positions; x = line-of-sight |
| `vxyz_g` | Quantity[km/s] | — | Particle velocities |
| `hsm_g` | Quantity[kpc] | — | SPH smoothing lengths (FWHM) |
| `coordinate_frame` | BaseRADecFrame | ICRS() | Coordinate frame |

**Key attributes:** `vsys`, `vhubble`, `distance`, `angular_diameter_distance`, `skycoords`, `pixcoords`

**Methods:** `rotate()`, `translate()`, `boost()`, `preview()`, `save_current_rotation()`, `save/load_affine_transformations()`

> **Redshift vs Distance mode:** Specify EITHER `redshift` (new) OR `distance` + `h` (legacy), not both.
> Redshift mode uses `Planck15` by default to compute `D_L` (flux) and `D_A` (angular scale).

#### 3.2 `TNGSource` — IllustrisTNG Source

Loads TNG simulation data from local disk and computes K13 HI/H2 masses.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `basePath` | str | — | Path to TNG output directory |
| `snapNum` | int | — | Snapshot number (99=z~0, 67=z~0.5, 50=z~1.0) |
| `subID` | int | — | Subhalo ID |
| `cosmology` | Cosmology or None | None | Cosmology (default: Planck15) |
| `uv_mode` | str | "global" | UV field: `"global"` (fast) or `"adaptive_fft"` (accurate) |
| `observed_redshift` | float or None | None | ★ Override observed redshift for distance |
| `load_subhalo_only` | bool | **True** | Load subhalo-only (fast) or full FOF group |
| `vpeculiar` | Quantity[km/s] | 0 | Peculiar velocity |
| `rotation` | Rotation or None | None | Source rotation |
| `L_coords` | L_coords or None | None | Angular momentum alignment |
| `ra`, `dec` | Quantity[deg] | 0 | Source centroid |
| `coordinate_frame` | BaseRADecFrame | ICRS() | Coordinate frame |

**Key attributes:** `source_redshift`, `snapshot_hubble`, `D_L`, `D_A`, `m_H2_g`

> **`observed_redshift`:** Decouples K13 physics (uses snapshot z for unit conversions) from
> cosmological distance (uses observed z for `D_L`, `D_A`). This lets you place a z~0 galaxy
> at any cosmological distance.

---

### 4. Sub-Modules

#### 4.1 `GaussianBeam` — Telescope Beam

```python
GaussianBeam(bmaj=15*U.arcsec, bmin=None, bpa=0*U.deg, truncate=4.0)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `bmaj` | Quantity[arcsec] | 15" | Beam major axis FWHM |
| `bmin` | Quantity[arcsec] | =bmaj | Beam minor axis FWHM |
| `bpa` | Quantity[deg] | 0° | Position angle (E of N) |
| `truncate` | float | 4.0 | Truncation radius in FWHM units |

#### 4.2 `GaussianNoise` — Noise Model

```python
GaussianNoise(rms=1*U.Jy/U.beam, seed=0)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `rms` | Quantity[Jy/beam] | 1 Jy/beam | Post-convolution RMS noise |
| `seed` | int | 0 | Random seed (reproducible) |

#### 4.3 Spectral Models

```python
# Gaussian line (fixed width or thermal)
GaussianSpectrum(sigma=7*U.km/U.s, ncpu=None, spec_dtype=np.float64)
# sigma="thermal" → uses sqrt(kB*T/m_p) per particle

# Dirac delta line
DiracDeltaSpectrum(ncpu=None, spec_dtype=np.float64)
```

#### 4.4 SPH Kernels

All kernel classes include **adaptive fallback**: when the pixel-to-smoothing-length ratio
makes the primary approximation inaccurate, they automatically fall back to
`DiracDeltaKernel` then `GaussianKernel(truncate=6)`.

| Kernel | min_valid_size | Typical simulation |
|--------|---------------|-------------------|
| `CubicSplineKernel()` | 1.16 px | Gadget/Gadget2 (classic) |
| `WendlandC2Kernel()` | 1.51 px | EAGLE |
| `WendlandC6Kernel()` | 1.29 px | Magneticum |
| `QuarticSplineKernel()` | 1.24 px | Colibre (SPHENIX) |
| `GaussianKernel(truncate=3.0)` | 2.34 px | General purpose |
| `DiracDeltaKernel()` | ≤0.5 px | Point-like particles |

> Use the **underscore-prefixed** versions (e.g., `_CubicSplineKernel`) to use the raw kernel
> without adaptive fallback.

---

### 5. Cosmology Utilities — `cosmology_utils`

```python
from astrovodka.cosmology_utils import *
```

| Function | Signature | Description |
|----------|-----------|-------------|
| `get_cosmology` | `(cosmo=None)` | Return cosmology instance (default Planck15) |
| `luminosity_distance` | `(redshift, cosmo=None)` | D_L [Mpc] for flux conversion |
| `angular_diameter_distance` | `(redshift, cosmo=None)` | D_A [Mpc] for angular scale |
| `comoving_distance` | `(redshift, cosmo=None)` | Comoving distance [Mpc] |
| `hubble_parameter` | `(redshift, cosmo=None)` | H(z) [km/s/Mpc] |
| `systemic_velocity_from_redshift` | `(redshift)` | Relativistic Doppler velocity [km/s] |
| `observed_frequency` | `(redshift, rest_freq=None)` | f_obs = f_rest/(1+z) [Hz] |
| `rho_code_to_physical` | `(rho_code, redshift, hubble_h)` | TNG density → g/cm³ |
| `rho_code_to_nH` | `(rho_code, redshift, hubble_h, X_H=0.76)` | TNG density → n_H [cm⁻³] |
| `calc_Jeans_length_physical` | `(u_code, rho_code, redshift, hubble_h)` | Jeans length [kpc] |
| `calc_Sigma_neutral_physical` | `(u_code, rho_code, f_neutral_XH, redshift, hubble_h)` | Surface density [M_sun/kpc²] |

---

### 6. K13 Model Functions — `_k13_model`

```python
from astrovodka.sources._k13_model import *
```

| Function | Signature | Description |
|----------|-----------|-------------|
| `calc_UMW_global_Faucher11` | `(redshift)` | Global UV background [Draine units] |
| `get_f_neutral_H` | `(nH_code, sfr_code)` | Estimate f_neutral from TNG fields |
| `f_mol_K13_vol` | `(rho_code, u_code, Z, nH_code, sfr_code, X_H, U_MW, redshift=0, hubble_h=0.6774, mgas=None)` | K13 H2 fraction → `(f_mol, f_HI, f_H2)` |
| `compute_HI_mass_particle` | `(gas_data, U_MW, redshift=0, hubble_h=0.6774)` | Per-particle HI/H2 masses [code units] |
| `compute_galaxy_fft_uv` | `(gas_data, redshift=0, hubble_h=0.6774, ngrid=128, f_esc=0.1, box_expand=1.1, return_grid=False)` | Adaptive FFT UV field. If `return_grid=True`, returns `(per_particle, grid, edges)` |

---

### 7. Data Products — `SubhaloDataProducts`

Container for 3D gridded data products of a TNG subhalo: UV field, DM/HI/H2 mass density
in comoving real-space and redshift-space. Grid resolution is configurable (default 128³).

```python
from astrovodka import SubhaloDataProducts
```

**Class Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `from_subhalo` | `(basePath, snapNum, subID, ngrid=128, uv_mode="global", load_subhalo_only=True, cosmology=None)` | Load TNG data and compute all 3D grids |
| `load_hdf5` | `(filename)` | Load products from HDF5 |

**Instance Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `save_hdf5` | `(filename)` | Save all grids to compressed HDF5 |

**Attributes (all `(ngrid, ngrid, ngrid)` float32 arrays):**

| Attribute | Unit | Description |
|-----------|------|-------------|
| `uv_field` | Draine units | UV radiation field (from global UVB or adaptive FFT) |
| `dm_real` | M_sun/(ckpc/h)³ | Comoving DM density (real space) |
| `dm_rsd` | M_sun/(ckpc/h)³ | Comoving DM density (redshift space) |
| `hi_real` | M_sun/(ckpc/h)³ | Comoving HI density (real space) |
| `hi_rsd` | M_sun/(ckpc/h)³ | Comoving HI density (redshift space) |
| `h2_real` | M_sun/(ckpc/h)³ | Comoving H2 density (real space) |
| `h2_rsd` | M_sun/(ckpc/h)³ | Comoving H2 density (redshift space) |
| `grid_edges` | tuple of 3 arrays | Bin edges [ckpc/h] for each axis |
| `metadata` | dict | redshift, hubble_h, scale_factor, etc. |

**Redshift-space formula:** `s_x = x + v_los / (a * H(a))` where `v_los` is the
line-of-sight peculiar velocity, `a = 1/(1+z)`, and `H(a)` from Planck15.

**Usage:**

```python
# Compute products
dp = SubhaloDataProducts.from_subhalo(
    basePath="/data/.../TNG/TNG-100/",
    snapNum=99, subID=8, ngrid=128,
    uv_mode="adaptive_fft",
)

# Save to HDF5
dp.save_hdf5("subhalo_8_products.hdf5")

# Load later
dp2 = SubhaloDataProducts.load_hdf5("subhalo_8_products.hdf5")

# Access
print(dp.dm_real.shape)  # (128, 128, 128)
print(dp.hi_rsd.sum())   # total HI mass density in redshift space
```

---

### 8. Coordinate Helpers

#### `L_coords` — Angular Momentum Alignment

```python
from astrovodka import L_coords

L_coords(incl=0*U.deg, pa=270*U.deg, az_rot=0*U.deg)
```

A `NamedTuple` specifying source orientation based on angular momentum.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `incl` | Quantity[deg] | 0° | Inclination angle (0° = face-on) |
| `pa` | Quantity[deg] | 270° | Position angle on sky (E of N, receding side) |
| `az_rot` | Quantity[deg] | 0° | Azimuthal rotation around the angular momentum axis |

**Algorithm:** First, the angular momentum vector of the central 1/3 of particles (HI-mass-weighted)
is rotated to the x-axis (face-on). Then `az_rot` rotates around L, `incl` tilts the disc, and
`pa` sets the sky position angle. All rotations are right-handed and extrinsic.

---

## Complete Example: TNG Subhalo with Data Products

```python
import astropy.units as U
from astrovodka import Martini, DataCube, SubhaloDataProducts
from astrovodka.sources import TNGSource
from astrovodka.beams import GaussianBeam
from astrovodka.spectral_models import GaussianSpectrum
from astrovodka.sph_kernels import CubicSplineKernel

basePath = "/data/.../TNG/TNG-100/"
snapNum = 99
subID = 8

# --- Compute 3D data products ---
dp = SubhaloDataProducts.from_subhalo(
    basePath=basePath, snapNum=snapNum, subID=subID,
    ngrid=64, uv_mode="global",
)
dp.save_hdf5(f"sub{subID}_products.hdf5")
print(f"DM real: {dp.dm_real.sum():.2e}")

# --- Create mock observation ---
source = TNGSource(
    basePath=basePath, snapNum=snapNum, subID=subID,
    observed_redshift=0.05, load_subhalo_only=True,
)

datacube = DataCube(
    n_px_x=128, n_px_y=128, n_channels=64,
    px_size=5*U.arcsec,
    channel_width=20*U.km/U.s,
    spectral_centre=source.vsys,  # automatically uses MHz by default
)

M = Martini(source=source, datacube=datacube,
            beam=GaussianBeam(bmaj=15*U.arcsec),
            spectral_model=GaussianSpectrum(sigma="thermal"),
            sph_kernel=CubicSplineKernel())

M.insert_source_in_cube(progressbar=False)
M.convolve_beam()
M.write_fits(f"sub{subID}_cube.fits")

# Check spectral axis is in MHz
print(f"Channel 0: {datacube.frequency_channel_edges[0].to(U.MHz):.4f}")
```

---

## Differences from MARTINI (astromartini)

| Feature | MARTINI | astrovodka |
|---------|---------|------------|
| Distance | `v = H0 * D` (linear) | `D_L`, `D_A` from `Planck15` |
| HI model | Pressure-based (Leroy+2008) | K13 volumetric (Diemer+2018) |
| TNG loading | Web API + JupyterLab | Local disk |
| Subhalo loading | FOF group only | FOF group or subhalo-only (default) |
| H2 masses | Not exposed | `source.m_H2_g` |
| Spectral axis | Velocity or frequency | **MHz by default** (`spectral_unit`) |
| Redshift override | — | `observed_redshift` parameter |
| 3D data products | — | `SubhaloDataProducts` (UV, DM, HI, H2) |
| Redshift space | — | RSD grid output |
