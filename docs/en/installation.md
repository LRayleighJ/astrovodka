# Installation and environment setup

This page covers everything needed to get astrovodka running: supported Python
versions, the three supported installation routes, the optional dependencies and
what each one unlocks, how to work offline on a cluster, how to lay out your
IllustrisTNG data, and how to check that the install is sound.

If you just want the shortest path, it is:

```bash
git clone https://github.com/LRayleighJ/astrovodka.git
cd astrovodka
conda env create -f environment.yml
conda activate astrovodka
```

and then jump to [Verifying the installation](#verifying-the-installation).

- [Requirements](#requirements)
- [Installation](#installation)
  - [Route 1: conda (recommended)](#route-1-conda-recommended)
  - [Route 2: virtualenv and pip](#route-2-virtualenv-and-pip)
  - [Route 3: core only, no optional dependencies](#route-3-core-only-no-optional-dependencies)
- [Optional dependencies](#optional-dependencies)
- [Offline and air-gapped installs](#offline-and-air-gapped-installs)
- [IllustrisTNG data](#illustristng-data)
- [Verifying the installation](#verifying-the-installation)
- [Troubleshooting](#troubleshooting)
- [Python version support](#python-version-support)

---

## Requirements

### Software

| Requirement | Minimum | Notes |
|---|---|---|
| Python | 3.9 | See [Python version support](#python-version-support) |
| pip | 21.3 | Needed for `pyproject.toml`-only projects |
| numpy | 1.20 | Core dependency |
| scipy | 1.6 | FFT convolution, interpolation, special functions |
| astropy | 5.0 | Cosmology, units, coordinates, WCS, FITS I/O |
| tqdm | 4.40 | Progress bars |
| h5py | 3.0 | *Optional but required for TNG data* — see below |
| matplotlib | 3.4 | Optional: figures and tutorials |
| multiprocess | 0.70 | Optional: `ncpu > 1` source insertion |

`numpy`, `scipy`, `astropy` and `tqdm` are installed automatically with the
package. The rest are opt-in extras, except that **you need `h5py` for any TNG
workflow** — it is what reads the snapshot files, not just what writes HDF5
output.

Note that astrovodka uses **`multiprocess`**, not the standard library
`multiprocessing`. MARTINI parallelises with `dill`-based serialisation, which
the stdlib module does not provide. Installing `multiprocess` is only necessary
if you pass `ncpu > 1`.

### Hardware

The package itself is light, but two features are memory-hungry:

- **Loading a TNG subhalo.** Cost scales with the subhalo's particle count.
  The tutorial's target, TNG-100 snapshot 99 subhalo 113, has 8,169 gas cells
  and loads in seconds. Cluster cores are a different matter: subhalo 0 of the
  same snapshot has roughly 38 million gas cells and needs well over 10 GB.
  Avoid low subhalo IDs, which are the most massive objects (the catalog is
  sorted by descending mass).
- **The 3D data products and the adaptive FFT UV field** build a full
  `ngrid³` grid (128³ by default), which is on the order of 10 MB per grid per
  array. Subhalos with more than ~10⁷ gas cells would benefit from chunking,
  which is not implemented.

Disk space for the software itself is negligible (a few MB). Disk space for
TNG data is not — see [IllustrisTNG data](#illustristng-data).

---

## Installation

astrovodka is not on PyPI. Install it from a checkout. All three routes below
start from:

```bash
git clone https://github.com/LRayleighJ/astrovodka.git
cd astrovodka
```

### Route 1: conda (recommended)

`environment.yml` in the repository root pins the stack that the project is
developed and verified against, and installs astrovodka itself in editable mode:

```bash
conda env create -f environment.yml
conda activate astrovodka
```

To use a different environment name, edit the `name:` line first, or override
it on the command line:

```bash
conda env create -f environment.yml -n my-env
```

To update an existing environment after pulling new commits:

```bash
conda activate astrovodka
conda env update -f environment.yml --prune
```

### Route 2: virtualenv and pip

Any Python 3.9+ environment works:

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

`[dev]` pulls in h5py, matplotlib, multiprocess, pytest and pytest-cov. For a
leaner install, name only the extras you need:

```bash
python -m pip install -e ".[hdf5]"                 # TNG work, nothing else
python -m pip install -e ".[hdf5,parallel]"        # + multi-core insertion
python -m pip install -e ".[hdf5,parallel,plotting]"   # + figures
```

If you prefer a plain requirements file over extras:

```bash
python -m pip install -r requirements-dev.txt
python -m pip install -e . --no-deps
```

Editable (`-e`) is recommended while the code is moving. For a fixed
installation use `pip install .` instead, but note that the examples and
tutorials assume the repository is present.

### Route 3: core only, no optional dependencies

The physics layer — `astrovodka.himodel`, which computes per-cell HI and H2
masses from a TNG gas cell array — depends only on numpy, scipy and astropy. It
can be used on its own, with no h5py, no matplotlib, and no data cube anywhere
in sight:

```bash
python -m pip install -e .
```

```python
from astrovodka.himodel import HIModelConfig, compute_hi_mass

cfg = HIModelConfig(h=0.6774, z=0.0, boxsize_ckpc_per_h=75000.0, ngrid=64)
result = compute_hi_mass(gas, model="K13_vol", cfg=cfg)
print(result.m_hi_msun, result.m_h2_msun)
```

The package is laid out so that this works: `astrovodka.himodel` imports
nothing from the observation layer, and `import astrovodka` does not import
h5py or matplotlib. A test enforces both (`tests/test_selfcontained.py`).

---

## Optional dependencies

| Extra | Package | What it unlocks | Needed for |
|---|---|---|---|
| `hdf5` | h5py | Reading TNG snapshot and group catalogs; `write_hdf5`, `save_state`, `SubhaloDataProducts.save_hdf5` | **Any TNG workflow** |
| `parallel` | multiprocess | `insert_source_in_cube(ncpu=N)`, multi-core spectrum initialisation | Large cubes |
| `plotting` | matplotlib | `Martini.preview()`, `GlobalProfile.plot_spectrum()`, the tutorial figures | Figures |
| `testing` | pytest, pytest-cov | The test suite | Development |
| `all` | all three runtime extras | — | Convenience |
| `dev` | all of the above + testing | — | Development |

`h5py` is listed as optional because the physics layer genuinely does not need
it. In practice, though, every route that touches actual TNG data requires it,
so treat it as mandatory unless you are only using `astrovodka.himodel`.

---

## Offline and air-gapped installs

Cluster login nodes often have no outbound network. Two things need attention.

**1. Build isolation.** By default pip creates a temporary virtual environment
to build the package and downloads `setuptools` into it. If you already have a
recent setuptools, skip that:

```bash
python -m pip install -e . --no-build-isolation
```

`pyproject.toml` requires `setuptools>=64`; check yours with
`python -c "import setuptools; print(setuptools.__version__)"`.

**2. Dependencies.** Either pre-download wheels on a networked machine,

```bash
# on a networked machine, matching the target Python
pip download astrovodka[hdf5,parallel,plotting] -d wheels/ --no-deps
pip download numpy scipy astropy tqdm h5py matplotlib multiprocess pytest -d wheels/
```

then copy `wheels/` across and install from it:

```bash
python -m pip install --no-index --find-links wheels/ -e . --no-build-isolation
```

or build a conda environment file on a networked machine and transfer it
(`conda env export > environment.lock.yml`, then `conda env create -f
environment.lock.yml`).

**3. Conda offline.** `conda env create --offline` works if the packages are
already in the local package cache. If your cluster provides TNG data but not a
scientific Python stack, `conda create -n astrovodka --clone base` followed by
`pip install -e . --no-deps` is often the fastest way in, provided numpy, scipy
and astropy are already present.

---

## IllustrisTNG data

astrovodka reads TNG data **from local disk**. It does not use the TNG web API
and does not need network access at run time, but it does need the files.

### Where to get it

- **TNG-100 / TNG-300 public data release:** <https://www.tng-project.org/data/>
  (registration required). Download the `groups_<snap>` and `snapdir_<snap>`
  directories for the snapshots you want.
- **TNG-50:** same site, smaller volumes, good for testing.
- If you only want the group catalogs (subhalo positions and masses) and not
  the particle data, `groups_<snap>` alone is enough to *choose* subhalos, but
  not to build cubes.

### Expected layout

Point `basePath` at the simulation directory, i.e. the one that *contains*
`groups_<snap>` and `snapdir_<snap>`:

```
/path/to/TNG/TNG-100/
├── groups_099/
│   ├── fof_subhalo_tab_099.0.hdf5      # group/subhalo catalog
│   └── ...
└── snapdir_099/
    ├── snap_099.0.hdf5                 # particle data, chunked
    └── ...
```

The snapshot number is zero-padded to three digits. Directory names must match
this pattern; `TNGSource` builds the paths from `basePath` and `snapNum`.

### Size

For TNG-100, one snapshot's particle data is on the order of **1–2 TB** and the
group catalog around **13 GB**. You do not need every snapshot — for most work
one or two is enough. Consider a symlink farm if the data lives on a scratch
filesystem with a purge policy.

### Telling astrovodka where the data is

There is no global configuration file. Pass the path explicitly:

```python
source = TNGSource(basePath="/path/to/TNG/TNG-100/", snapNum=99, subID=113)
```

The example scripts and both tutorial notebooks read it from an environment
variable instead, so that no path is baked into the repository:

```bash
export TNG_BASE_PATH=/path/to/TNG/TNG-100/
python examples/en/astrovodka_demo.py --snapNum 99 --subID 113
```

---

## Verifying the installation

### 1. It imports

```bash
python -c "import astrovodka; print(astrovodka.__version__)"
python -c "from astrovodka.himodel import METHODS; print(METHODS)"
```

The second should print `('GD14_vol', 'GD14_map', 'K13_vol', 'K13_map')`.

### 2. The optional pieces are there

```bash
python -c "import h5py, matplotlib, multiprocess; print('optional deps OK')"
```

A missing `h5py` will not stop `import astrovodka` — by design — but it will
fail the moment you load a TNG subhalo.

### 3. The test suite passes

```bash
python -m pytest tests/ -q
```

The suite is synthetic: **it needs no TNG data, no network, and no GPU**, and
takes about 15 seconds. It covers the four HI models against analytic
expectations, the cosmology helpers, the data-cube spectral-axis conversions,
the 3D data products, and the self-containment rules described above.

### 4. The examples run

Without any TNG data:

```bash
python examples/en/astrovodka_demo.py
```

This builds an analytic rotating disc at z = 0.01 and writes a FITS cube to
`examples/outputs/`. It prints a mass-conservation check; the line

```
Mass in cube (assuming distance ... ): 1.00e+11 solMass
    [100% of initial source mass]
```

should read 100%.

With TNG data:

```bash
export TNG_BASE_PATH=/path/to/TNG/TNG-100/
python examples/en/astrovodka_demo.py --snapNum 99 --subID 113 --hiModel K13_vol
```

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: No module named 'astrovodka'` | Not installed, or installed into a different interpreter | `python -m pip install -e .` using the *same* `python` you run the script with |
| `ModuleNotFoundError: No module named 'h5py'` when constructing `TNGSource` | h5py is an optional extra | `pip install -e ".[hdf5]"` |
| `ModuleNotFoundError: No module named 'six'` | Older revisions of the vendored TNG I/O used `six` | Updated in 1.0.0; update your checkout |
| `ModuleNotFoundError: No module named 'typing_extensions'` on Python 3.9/3.10 | `typing.Self` is only in the 3.11 stdlib | `pip install typing_extensions`, or reinstall the package so the conditional dependency is resolved |
| `ModuleNotFoundError: No module named 'multiprocess'` when using `ncpu > 1` | Stdlib `multiprocessing` is not a substitute | `pip install -e ".[parallel]"` |
| `ImportError: numpy.core.multiarray failed to import` | numpy/scipy built against mismatched ABIs | `pip install --force-reinstall numpy scipy` in a fresh environment |
| `ValueError: datacube padding insufficient for beam convolution` | `write_fits` calls `drop_pad()`, so writing mid-pipeline removes the padding | `M.datacube.add_pad(M.beam.needs_pad())` before the next `convolve_beam()` |
| `FileNotFoundError: .../snapdir_099/...` | `basePath` points at `snapdir_099/` itself, or at the parent of `TNG-100/` | `basePath` must be the directory *containing* `groups_<snap>` and `snapdir_<snap>` |
| Loading a subhalo takes minutes and hundreds of GB | You are loading a cluster core | Use a subhalo ID in the 100s or higher; see [Hardware](#hardware) |
| `RuntimeError: unknown hi_model` | A typo, or a name from an older revision | One of `K13_vol`, `K13_map`, `GD14_vol`, `GD14_map` |

---

## Python version support

**Python 3.9 is the floor and the version this release is tested on**, together
with 3.10 through 3.13 in CI. `pyproject.toml` declares `requires-python =
">=3.9"`.

Two details are worth knowing:

- `astrovodka/datacube.py` uses `typing.Self`, which only reaches the standard
  library in Python 3.11. On 3.9 and 3.10 that falls back to `typing_extensions`,
  which is declared as a conditional dependency, so it is installed
  automatically on those versions.
- `types.EllipsisType` (3.10+) is imported defensively, with a
  `type(Ellipsis)` fallback, so the annotation works on 3.9.

Because astropy 7 requires Python 3.11 and later, on Python 3.9 pip will
resolve an older astropy (6.x), which is supported. The declared floor for
astropy is 5.0.

The vendored MARTINI observation machinery is not restyled here, but it is
3.9-clean and CI runs the whole suite on every supported version.
