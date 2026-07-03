# astrovodka

Synthetic HI (21 cm) data cubes from IllustrisTNG cosmological SPH simulations,
with proper cosmological distance handling and K13 volumetric HI/H2 modeling.

**astrovodka** is built upon [MARTINI](https://github.com/kyleaoman/martini)
(`astromartini`), improving its distance calculation and gas mass modeling for
high-redshift applications.

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPLv3-green.svg)](LICENSE)

---

## Installation

### Requirements

- Python ≥ 3.10
- numpy, scipy, astropy, tqdm

### Install from source

```bash
git clone <repo-url>
cd astrovodka.0.1
pip install -e .
```

### Dependencies

All dependency libraries:

| Package | Required | Used for |
|---------|----------|----------|
| `numpy` | ✓ | Array operations throughout |
| `scipy` | ✓ | FFT convolution (`fftconvolve`), interpolation (`RegularGridInterpolator`), optimization (`fsolve`), special functions (`erf`) |
| `astropy` | ✓ | Cosmology (`Planck15`), units (`Quantity`), coordinates (`SkyCoord`, `SpectralCoord`), WCS (`astropy.wcs`), FITS I/O (`astropy.io.fits`), constants |
| `tqdm` | ✓ | Progress bars |
| `multiprocess` | — | Multi-core parallel source insertion (`ncpu > 1`) |
| `h5py` | — | HDF5 output (`write_hdf5`, `save_state`, `SubhaloDataProducts`) |
| `matplotlib` | — | Preview figures (`preview`, `plot_spectrum`) |

```bash
# Minimal install (required only)
pip install -e .

# With all optional dependencies
pip install -e ".[dev]"
```

---

## Key Features

- **Mock HI data cubes** from SPH particle data (positions, velocities, HI masses, temperatures)
- **IllustrisTNG support** — load subhalos directly from local TNG snapshot files
- **Telescope beam convolution** (Gaussian, extensible)
- **Gaussian noise** addition with configurable RMS
- **Configurable spectral line models** (Gaussian with thermal broadening, Dirac delta)
- **Multiple SPH smoothing kernels** (cubic spline, Wendland C2/C6, quartic spline, Gaussian) with adaptive accuracy fallback
- **FITS and HDF5 output** with full WCS headers
- **Global spectrum mode** (spatially integrated, no beam/noise)
- **3D gridded data products** (UV field, DM/HI/H2 density in real and redshift space)

### Quick Example

```python
import astropy.units as U
from astrovodka import Martini, DataCube
from astrovodka.sources import TNGSource
from astrovodka.beams import GaussianBeam
from astrovodka.spectral_models import GaussianSpectrum
from astrovodka.sph_kernels import CubicSplineKernel

# Load a TNG subhalo with K13 HI modeling
source = TNGSource(
    basePath="/data/.../TNG/TNG-100/",
    snapNum=99, subID=400547,
    observed_redshift=0.05,    # place at custom distance
    load_subhalo_only=True,     # fast: subhalo particles only
)

# Build data cube (spectral axis defaults to MHz)
datacube = DataCube(
    n_px_x=256, n_px_y=256, n_channels=64,
    px_size=5 * U.arcsec,
    channel_width=20 * U.km / U.s,
    spectral_centre=source.vsys,
)

# Run pipeline
M = Martini(source=source, datacube=datacube,
            beam=GaussianBeam(bmaj=15 * U.arcsec),
            spectral_model=GaussianSpectrum(sigma="thermal"),
            sph_kernel=CubicSplineKernel())
M.insert_source_in_cube()
M.convolve_beam()
M.write_fits("output.fits")
```

---

## What's New Compared to MARTINI (`astromartini`)

| Feature | MARTINI | astrovodka |
|---------|---------|------------|
| **Distance calculation** | Linear Hubble law `v = H₀ × D` | `astropy.cosmology.Planck15`: proper D<sub>L</sub>, D<sub>A</sub>, H(z) |
| **HI mass model** | Pressure-based (Leroy+2008) | **K13 volumetric model** (Krumholz 2013, Diemer+2018) with adaptive per-galaxy FFT UV field |
| **High-redshift support** | z ≈ 0 only | Works at any redshift |
| **Spectral axis** | Velocity or frequency | **MHz by default** (`spectral_unit='MHz'`) |
| **Redshift override** | — | `observed_redshift` parameter: place a z≈0 galaxy at any cosmological distance |
| **H₂ masses** | Not exposed | `source.m_H2_g` per-particle H₂ mass from K13 |
| **TNG data access** | Web API + JupyterLab | **Local disk** (offline, fast) |
| **Subhalo loading** | Entire FOF group | **Subhalo-only** (fast, default) or full group |
| **3D data products** | — | `SubhaloDataProducts`: UV field, DM/HI/H₂ density grids in real + redshift space, with HDF5 I/O |

---

## Credits & Citation

### MARTINI (underlying framework)

astrovodka is built upon and includes modified code from **MARTINI**
([github.com/kyleaoman/martini](https://github.com/kyleaoman/martini)).
If your use of astrovodka leads to a publication, please cite the original MARTINI papers:

- Oman, K. A. 2024, *J. Open Source Softw.*, 9, 6860. [doi:10.21105/joss.06860](https://doi.org/10.21105/joss.06860)
- Oman, K. A., et al. 2019, *MNRAS*, 482, 821. [doi:10.1093/mnras/sty2687](https://doi.org/10.1093/mnras/sty2687)
- Oman, K. A. 2019, *ASCL*, ascl:1911.005. [ascl.net/1911.005](https://ascl.net/1911.005)

### MARTINI License

MARTINI is distributed under the **GNU General Public License v3 (GPLv3)**.
See the [MARTINI LICENSE](https://github.com/kyleaoman/martini/blob/main/LICENSE.md) for details.
astrovodka inherits this license for all derivative code.

### astrovodka

For questions regarding astrovodka-specific features (cosmological distances, K13 HI model,
observed_redshift, data products), please contact:

📧 **liuzr@bao.ac.cn**

---

## Documentation

Full API documentation: see [`docs/html/index.html`](docs/html/index.html) or [`docs/md/index.md`](docs/md/index.md).

Example notebook: [`example_vodka.ipynb`](example_vodka.ipynb) — processes 50 TNG-100 subhalos
at three redshifts with K13 mass comparison and data cube visualization.
