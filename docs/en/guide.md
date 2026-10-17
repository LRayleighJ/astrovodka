# User guide

How astrovodka works, and how to drive it. For installation see
[installation.md](installation.md); for per-class signatures see the
[API reference](api/index.md); for the HI/H₂ models in depth see
[hi-models.md](hi-models.md).

- [The pipeline](#the-pipeline)
- [Sources](#sources)
- [Two distances](#two-distances)
- [The data cube](#the-data-cube)
- [Sub-modules](#sub-modules)
- [A complete worked example](#a-complete-worked-example)
- [Common tasks](#common-tasks)
- [Pitfalls](#pitfalls)
- [Differences from MARTINI](#differences-from-martini)

---

## The pipeline

astrovodka is object-oriented: `Martini` holds a **source**, a **data cube**,
and three sub-modules, and runs them in order.

```
    source ────────┐
    datacube ──────┤
    beam ──────────┤     ┌──────────────────────┐
    noise ─────────┼────▶│       Martini        │
    spectral_model ┤     └──────────────────────┘
    sph_kernel ────┘        │      │       │
                            ▼      ▼       ▼
                   insert_source_in_cube()  →  particles onto pixels
                            │
                      add_noise()           →  (optional) noise per pixel
                            │
                    convolve_beam()         →  (optional) telescope beam
                            │
                      write_fits()          →  a cube you can open in DS9
```

Each step is a separate method, so you can stop after any of them. Stopping
after `insert_source_in_cube()` gives a noiseless, unconvolved cube — which is
what you want when checking the physics, because there is then no noise to
compete with the signal and no beam to redistribute it.

## Sources

A source abstracts "the thing that emits 21 cm line radiation". Two are
implemented.

### `SPHSource` — any SPH data

```python
SPHSource(
    redshift=0.01, ra=30*U.deg, dec=45*U.deg,
    T_g=..., mHI_g=..., xyz_g=..., vxyz_g=..., hsm_g=...,
)
```

The particle arrays are ordinary `astropy` quantities:

| Array | Unit | Meaning |
|---|---|---|
| `xyz_g` | kpc | Positions. **`x` is the line of sight**; `y`-`z` is the sky plane |
| `vxyz_g` | km/s | Velocities, same axis convention |
| `mHI_g` | M☉ | HI mass per particle |
| `T_g` | K | Temperature, used by thermal line widths |
| `hsm_g` | kpc | SPH smoothing length **as an FWHM**, not a kernel radius |

That last point deserves emphasis. Simulations disagree about what a "smoothing
length" is — some store the radius at which the kernel reaches zero, some store
the FWHM, some store something else again. MARTINI chose FWHM as its
convention. If you feed it a kernel radius instead, every particle will be
smoothed over a different-sized region than intended and nothing will raise an
error.

### `TNGSource` — IllustrisTNG, from local disk

```python
TNGSource(
    basePath="/path/to/TNG/TNG-100/",
    snapNum=99, subID=113,
    hi_model="K13_vol",
    observed_redshift=0.02,
)
```

`TNGSource` reads the snapshot directly off the filesystem, runs the chosen
HI/H₂ model, and passes the result to `SPHSource`. It adds:

- **`hi_model`** — one of the four models (see [hi-models.md](hi-models.md)).
  Validated *before* any disk I/O, so a typo fails immediately rather than after
  loading several GB.
- **`observed_redshift`** — see [Two distances](#two-distances).
- **`load_subhalo_only`** (default `True`) — load just the subhalo's particles
  rather than its entire FOF group. Much faster, and usually what you want.
- **DM centre-of-mass centring**, rather than the catalog's `SubhaloPos`.

After construction, these are the attributes worth knowing:

| Attribute | Meaning |
|---|---|
| `m_HI_total` | Total HI mass [M☉] |
| `m_H2_g` | Per-particle H₂ mass |
| `f_mol_g`, `f_neutral_g` | Per-particle fractions |
| `u_mw_g` | Per-particle UV field, in Draine units |
| `hi_diagnostics` | Model diagnostics, including `f_gas_in_map` |
| `D_L`, `D_A`, `vsys` | Cosmological distances and systemic velocity |
| `source_redshift` | The snapshot's own redshift |

## Two distances

This is the part that most distinguishes astrovodka from MARTINI, and getting it
wrong is expensive, so it is worth being explicit.

**`D_L` and `D_A` are not interchangeable, and they differ by `(1+z)²`.**

| Distance | Used for | Getting it wrong costs |
|---|---|---|
| `D_L`, luminosity distance | Converting an HI mass to a flux; normalising the cube | `(1+z)⁴` |
| `D_A`, angular diameter distance | Every angular size: pixel scale, source extent, beam | `(1+z)²` |

astrovodka computes both from `astropy.cosmology.Planck15` (or a cosmology you
pass in) and applies each where it belongs.

### The two redshifts

A TNG snapshot has one redshift, and it is not always the redshift you want to
observe at. A `z ≈ 0` snapshot has a vanishing angular diameter distance, so a
galaxy in it subtends essentially no angle and cannot be resolved at all.

`observed_redshift` separates the two roles:

| Quantity | Comes from | Because |
|---|---|---|
| Unit conversions, UV background, the HI/H₂ model | the **snapshot's** redshift | that is the physical state of the gas |
| `D_L`, `D_A`, `vsys`, angular sizes, flux | `observed_redshift` | that is the geometry of the observation |

```python
source = TNGSource(basePath=..., snapNum=99, subID=113, observed_redshift=0.05)
print(source.source_redshift)   # ~0, the snapshot
print(source.D_L, source.D_A)   # computed at z = 0.05
```

### The systemic velocity convention

`vsys` uses the **radio convention**, `c·z/(1+z)`, which is what the WCS spectral
axis expects. This matters more than it sounds: at z = 0.5 the alternative
(`H(z)·D_C`) is 74% too large, exceeds `c` beyond z ≈ 0.79, and produces negative
frequencies at high redshift.

## The data cube

```python
DataCube(
    n_px_x=256, n_px_y=256, n_channels=64,
    px_size=5*U.arcsec,
    channel_width=20*U.km/U.s,
    spectral_centre=source.vsys,
)
```

Two things about the spectral axis are worth knowing.

**The unit follows `channel_width` by default.** Give a velocity width and you
get a velocity axis; the WCS records the conversion. Pass
`spectral_unit="MHz"` to force a frequency axis instead.

**Widths are converted differentially.** A velocity width crossing to a
frequency axis is converted as `Δf = f_rest·Δv/c`, not by transforming the
endpoints and subtracting — the two differ, and using the wrong one silently
discards almost all the signal (see [Pitfalls](#pitfalls)).

The spectral axis always spans the line, not just the gas. `Martini` prunes
particles whose line profile does not fit inside the channel range, and the
profile is `±4σ` wide, so the axis must be roughly `8σ` wider than the particle
velocity range. For a star-forming subhalo the hottest cells reach σ ≈ 100 km/s,
which means about 800 km/s of slack. Get this wrong and you lose mass silently —
the recovery line drops below 100% and nothing else complains.

The axis can be read back in either unit regardless of which one it was built in:

```python
datacube.velocity_channel_edges   # km/s
datacube.frequency_channel_edges  # Hz
```

## Sub-modules

| Component | Options | Notes |
|---|---|---|
| **Beam** | `GaussianBeam(bmaj, bmin, bpa, truncate)` | `bmin` defaults to `bmaj` |
| **Noise** | `GaussianNoise(rms, seed)` | Added *before* convolution, so the beam correlates it |
| **Spectral model** | `GaussianSpectrum(sigma=...)`, `DiracDeltaSpectrum()` | `sigma="thermal"` uses `sqrt(k_B T / m_p)` per particle |
| **SPH kernel** | `CubicSplineKernel`, `WendlandC2Kernel`, `WendlandC6Kernel`, `QuarticSplineKernel`, `GaussianKernel`, `DiracDeltaKernel` | Pick the one matching your simulation's kernel |

The SPH kernels are **adaptive**: each falls back to `DiracDeltaKernel` and then
`GaussianKernel(truncate=6)` when the pixel-to-smoothing-length ratio makes the
primary approximation inaccurate, and says so. Validation warns when the ratio
is marginal; `skip_validation=True` on `insert_source_in_cube()` silences it.

## A complete worked example

```python
import astropy.units as U
from astrovodka import Martini, DataCube, SubhaloDataProducts
from astrovodka.sources import TNGSource
from astrovodka.beams import GaussianBeam
from astrovodka.spectral_models import GaussianSpectrum
from astrovodka.sph_kernels import CubicSplineKernel

basePath, snapNum, subID = "/path/to/TNG/TNG-100/", 99, 113

# --- 3D data products: UV field, and DM/HI/H2 densities in real and redshift space
dp = SubhaloDataProducts.from_subhalo(
    basePath=basePath, snapNum=snapNum, subID=subID,
    ngrid=64, hi_model="K13_vol",
)
dp.save_hdf5(f"sub{subID}_products.hdf5")

# --- The source
source = TNGSource(
    basePath=basePath, snapNum=snapNum, subID=subID,
    hi_model="K13_vol", observed_redshift=0.05, load_subhalo_only=True,
)

# --- The cube. The spectral axis adopts channel_width's units (km/s here).
datacube = DataCube(
    n_px_x=128, n_px_y=128, n_channels=64,
    px_size=5 * U.arcsec,
    channel_width=20 * U.km / U.s,
    spectral_centre=source.vsys,
)

# --- Run it
M = Martini(source=source, datacube=datacube,
            beam=GaussianBeam(bmaj=15 * U.arcsec),
            spectral_model=GaussianSpectrum(sigma="thermal"),
            sph_kernel=CubicSplineKernel())
M.insert_source_in_cube(progressbar=False)
M.convolve_beam()
M.write_fits(f"sub{subID}_cube.fits")

# --- Read the axis back in either unit
print(f"channel 0: {datacube.frequency_channel_edges[0].to(U.MHz):.4f} MHz")
print(f"channel 0: {datacube.velocity_channel_edges[0].to(U.km/U.s):.1f} km/s")
```

`insert_source_in_cube()` prints the mass it put in the cube as a fraction of
the source total. **Read that line.** It should say 100%.

## Common tasks

### Change the angular resolution

Rebuild the cube with a different beam and pixel size. Interferometers resolve
structure but miss extended emission; single dishes do the opposite.

```python
datacube = DataCube(..., px_size=beam_fwhm / 3, ...)   # ~3 pixels per beam
```

### Put a galaxy at several redshifts

```python
for z in (0.0, 0.02, 0.05, 0.1):
    source = TNGSource(basePath=..., snapNum=99, subID=113, observed_redshift=z)
    ...
```

Build a fresh `TNGSource` each time. Building a cube **prunes the source in
place**, so one source object is good for one observation — see
[Pitfalls](#pitfalls).

### Process many subhalos

See `examples/en/batch_subhalos.py`. The short version: read the group catalog,
filter on gas-cell count, then loop.

```python
import h5py
with h5py.File(basePath + "groups_099/fof_subhalo_tab_099.0.hdf5", "r") as f:
    n_gas = f["SubhaloLenType"][:, 0]
ids = [i for i, n in enumerate(n_gas) if 64 <= n <= 5000]
```

**Do not** take the first N of the catalog — it is sorted by descending mass, so
the first entries are cluster cores. Filter first.

### Match an existing observation

```python
datacube = DataCube.from_wcs(fits.getheader("observation.fits")[0])
```

The cube then inherits the real pixel grid, channel width and spectral
reference system.

### A spectrum instead of a cube

`GlobalProfile` skips the cube entirely — no beam, no noise, no SPH kernel:

```python
from astrovodka import GlobalProfile
gp = GlobalProfile(source=source, spectral_model=GaussianSpectrum(sigma=7*U.km/U.s),
                   n_channels=64, channel_width=10*U.km/U.s,
                   spectral_centre=source.vsys)
gp.plot_spectrum()
```

## Pitfalls

These four cost real debugging time. All are inherited from MARTINI.

### 1. `Martini` prunes the source in place

`Martini.__init__` calls `_prune_particles`, which removes out-of-cube particles
from the **source object itself** via `source.apply_mask` — not just from that
one observation. Build a small throwaway cube to check something, then a big
science cube from the same source, and the second one silently gets a truncated
source. The only symptom is that the recovery line drops below 100%.

**One source, one observation.** Rebuild the source for each cube. (The pruning
does not mask `m_H2_g` along with `mHI_g`, so a pruned source is also left
internally inconsistent.)

### 2. `write_fits` and `write_hdf5` drop the padding

Both call `datacube.drop_pad()`, because the padding exists for convolution.
Write an intermediate product mid-pipeline and the next `convolve_beam()` raises
`ValueError: datacube padding insufficient for beam convolution`. Restore it:

```python
M.datacube.add_pad(M.beam.needs_pad())
```

### 3. A spectral model cannot be shared between cubes of different shape

`GaussianSpectrum` caches its `spectra` array, and `Martini` only re-initialises
when that array is `None`. Reusing one instance for a 32-channel cube and then a
42-channel one raises `cannot reshape array of size ...`.

**Build one spectral model per `Martini`.** (`M.init_spectra()` is the escape
hatch if you have already hit it.)

### 4. The spectral axis must span the line, not just the gas

Covered under [The data cube](#the-data-cube): allow about `8σ` beyond the
particle velocity range, or about 800 km/s for a star-forming subhalo.

## Differences from MARTINI

The observation machinery is MARTINI's: astrovodka imports the `astromartini`
package, subclasses it where the behaviour below differs, and re-exports the
rest unchanged. What differs:

| | MARTINI | astrovodka |
|---|---|---|
| **Distances** | Linear Hubble law `v = H₀D` | `Planck15`: `D_L` for flux, `D_A` for angular scales |
| **Systemic velocity** | `h·100·D` | Radio convention `c·z/(1+z)` |
| **HI mass** | Pressure-based (Leroy+2008) | Four selectable Diemer+2018 models |
| **H₂ mass** | Not exposed | `source.m_H2_g`, per particle |
| **UV field** | — | Adaptive per-galaxy FFT field, floored at the UV background |
| **Redshift override** | — | `observed_redshift`, decoupling physics from geometry |
| **TNG access** | Web API + JupyterLab | Local disk |
| **Subhalo loading** | Full FOF group | Subhalo-only (default) or the full group |
| **Centring** | Catalog `SubhaloPos` | DM centre of mass |
| **3D products** | — | `SubhaloDataProducts`: UV, DM, HI, H₂; real and redshift space |
| **Sources** | EAGLE, Simba, FIRE, Magneticum, Colibre, SWIFT, … | TNG only |

Both `distance` + `h` (MARTINI's interface) and `redshift` + `cosmology` are
accepted by `SPHSource`. In distance mode `D_L == D_A == distance`, so that path
behaves exactly as upstream does.

### What is not implemented

- **Only TNG sources are re-exported.** Upstream's EAGLE, Simba, FIRE,
  Magneticum, Colibre, SWIFT and `CombinedSource` classes remain available from
  `martini.sources`, but `astrovodka.sources` deliberately does not wrap sources
  the project cannot test.
- **Multi-galaxy observations in one cube**, which is what `CombinedSource` is
  for.
- **Memory chunking** for large subhalos — the 3D grids are built whole.
