# Changelog

All notable changes to astrovodka are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] — 2026-10-18

### Changed

- **astrovodka now depends on and calls MARTINI instead of containing a copy of
  it.** Following the MARTINI authors' request, the vendored observation
  machinery was removed and `astromartini>=2.1.9` became a required dependency.
  `astrovodka.datacube.DataCube`, `astrovodka.martini.Martini`/`GlobalProfile`,
  `astrovodka.sources.SPHSource` and the two spectral models are now thin
  subclasses of their MARTINI counterparts; `beams`, `noise`, `sph_kernels`,
  `L_coords` and `sources._illustris_tools` re-export MARTINI's classes.
  Deleting `sources/_L_align.py` and `sources/_extra_cartesian_transforms.py`
  removes roughly 2,500 lines of upstream code from the repository.

  astrovodka's own behaviour is unchanged: the redshift distance mode (`D_L` for
  flux, `D_A` for angular scales, radio-convention `vsys`), the `spectral_unit`
  cube option and the luminosity-distance flux normalisation all remain, and are
  covered by the existing cosmology and mass-recovery tests. A rotation adapter
  detects MARTINI's 2.1.16 API change at import time, so every released version
  from 2.1.9 onward works; Python 3.9 remains the floor.

- `tests/test_selfcontained.py` → `tests/test_upstream_dependency.py`: the
  invariant it enforces is inverted. It now checks that `astromartini` is
  declared and imported, that every re-exported class is an upstream object or
  a subclass of one rather than a reimplementation, and that `himodel` imports
  neither MARTINI nor the observation layer.

- `typing_extensions` is now declared for Python 3.9/3.10 because MARTINI's
  `datacube` imports `typing.Self` with a `typing_extensions` fallback without
  declaring it; `requests` is a core dependency because MARTINI 2.1.16+ imports
  it when `martini.sources` loads.

- `six` is listed alongside `h5py` in the `hdf5` extra, because MARTINI's
  IllustrisTNG I/O helpers import it.

### Fixed

- **`SubhaloDataProducts`' redshift-space grids retained 0.5% of their mass.**
  Two mistakes compounded, and both were silent:

  1. the line-of-sight shift was applied to *absolute* velocities rather than to
     velocities relative to the object, and
  2. the shifted particles were then histogrammed into the *real-space* box.

  At z = 0 one km/s is 10 ckpc/h, so a galaxy's bulk motion alone displaces every
  particle by several Mpc — far outside a box ~40 ckpc/h across. Measured on
  subhalo 113 (TNG-100 snapshot 99, ngrid = 64): `dm_rsd` summed to 3.4e8
  against `dm_real`'s 5.9e10, i.e. **0.58% of the mass in 82 non-zero cells out
  of 9384**, and `hi_rsd`/`hi_real` was 0.5%. The tutorial's
  `sub113_..._UV_and_RSD.png` panel, which plotted `dm_rsd - dm_real`, was
  therefore showing a real map minus almost nothing.

  The shift is now taken in the object's own frame (the DM centre of mass, the
  same reference the box is centred on), and redshift space gets its own box,
  `grid_edges_rsd`, whose line-of-sight axis is sized to hold the result —
  2722 ckpc/h against 41 ckpc/h for subhalo 113, which is the fingers-of-God
  effect at its real size, not an artefact. `metadata['f_mass_in_rsd_box_dm']`
  (and `_hi`, `_h2`) now report what the box keeps and are 1.000.

  **The real-space grids are untouched.** They are built before this code runs,
  from the same expressions as before; verified bit-identical against the
  previous commit (`dm_real`, `hi_real`, `h2_real` sums and all six grid bounds
  unchanged), because a redshift-space choice must not reshape real-space data.
  The consequence is that `*_rsd` and `*_real` no longer share a box and so can
  no longer be subtracted cell by cell — the tutorial panel now compares their
  line-of-sight profiles instead, which shows the same physics with the
  bookkeeping visible.

  The box sizing is shared with `GroupDataProducts` through
  `data_products._rsd_los_halfwidth`, so there is one implementation of it
  rather than two.

### Added

- **`astrovodka.group_products` — FOF group / cluster-scale data products.**
  `SubhaloDataProducts` grids one subhalo; `GroupDataProducts.from_group()`
  grids a whole FOF group — every subhalo in it plus the group's unbound fuzz —
  and splits the gas into **ISM** (bound to a subhalo and within
  `min(k_gas R_half,gas, k_star R_half,star)` of its centre), **CGM** (bound but
  outside that radius) and **IGM** (the fuzz). Grids cover dark matter, gas and
  HI/H2 by class, with per-cell arrays and the member subhalo table alongside
  them so that a re-classification or a zoomed re-grid needs no second pass over
  the snapshot. Shipped with `examples/{en,zh}/group_cluster_viz.py` and
  `tests/test_group_products.py`.

  Membership comes from the group catalog offsets and is exact, which took three
  details to get right, all of them documented where they bite:

  1. **Lengths must come from `SubhaloLenType`, never from differencing
     `Subhalo/SnapByType`.** `SnapByType` is a cumulative start-offset table, so
     for the *last* subhalo of a group its next row belongs to the next group
     and the difference silently absorbs that group's fuzz. On TNG-100 snapshot
     99 group 59 the last subhalo has `SubhaloLenType == 0` yet a 319,249-cell
     step. With the catalog lengths, every block's summed mass reproduces
     `SubhaloMassType` at a ratio of 1.0000 (68 of 68 gas-bearing subhalos,
     worst deviation 2e-4); `from_group` now checks this and refuses to continue
     if it fails.
  2. **A gas-only ISM radius does not survive group scale.** A group's central
     subhalo contains the whole hot halo, so its `R_half,gas` is 263 ckpc/h
     against 19 ckpc/h for `R_half,star` (group 59): `k x R_half,gas` at k=0.5
     already calls 14.7% of the group's gas mass ISM. `ism_radius()` therefore
     takes the minimum of the two.
  3. **Redshift space needs the group's bulk velocity removed.** At z=0 one km/s
     is 10 ckpc/h, so a group velocity of a few hundred km/s displaces every
     particle by Mpc and the map comes out empty; subtracting it first, group 59
     spans 17 cMpc along the line of sight against 1.9 cMpc on the sky. That
     geometry is why `compute_rsd` defaults to False.

  **Real space is classified first, and only then transformed.** The classes are
  assigned on the real-space coordinates and carried into redshift space as a
  mask, so nothing is ever classified in velocity space. Redshift space also
  gets its own grid box (`grid_edges_rsd`) rather than sharing the real-space
  one, which has two consequences, both verified on group 59: the per-cell
  arrays, the class masses and every real-space grid are **bit-identical with
  `compute_rsd` on or off**, and each class keeps its mass exactly through the
  transform (real-space and redshift-space grid masses agree at 1.000000 for
  ISM, CGM and IGM). Sharing one box would instead have stretched the
  real-space line-of-sight axis by 9x, degrading the real-space products to
  describe a transform that is not about them.

  HI is computed **per subhalo**, each in its own adaptive FFT box — the
  prescription the volumetric models were calibrated under — and only for
  gas-bearing subhalos (68 of group 59's 1710 members), with the fuzz modelled
  separately against the group-wide UV field.

  The bilingual tutorials gained a matching section (**6.3**), which uses the
  subhalo they follow as the way in: subhalo 113 turns out to live in group 0,
  the most massive cluster in the box (17,185 subhalos, 4.5e7 gas cells), so the
  section demonstrates the cheap way to check a group's size before loading it
  and then works through a group that is affordable in a notebook. The section
  asserts its own results (the classes partition the group exactly, and
  `f_HI(ISM) > f_HI(CGM) > f_HI(IGM)`), so a change in the physics shows up as a
  failure there rather than as a quietly different figure.

- `himodel.physics.calc_temperature(u_code, X_H, xe)` — the TNG internal-energy
  to temperature conversion, extracted from `TNGSource` so the group products
  and the source share one expression. `TNGSource` now calls it; the two
  formulations agree to 6.7e-16 relative.

### Fixed

- `ruff check` passed locally on the narrow rule set but failed in CI with 15
  findings, all of them genuine: unused imports (`astropy.constants` in
  `data_products.py` and `sph_source.py`, `os` in `martini.py` and
  `tng_source.py`, `numpy` in `test_cosmology.py`, a `TYPE_CHECKING` block in
  `sph_source.py` whose names appeared only in docstrings), a dead `rho_g`
  computation in `tng_source.py`, and three unused local variables.

  The `sph_source.py` case was worth checking rather than deleting: `L_coords`
  is imported from `..L_coords` and also used as a parameter name in
  `__init__` and `rotate`, so ruff reported it as both unused and redefined.
  Every use is the parameter — the rotation itself is done by `L_align` from
  `._L_align` — so the import was a leftover, and removing it is safe. Verified
  by exercising `rotate(L_coords=...)` directly.

  In `tests/test_cosmology.py`, the unused `expected` was the symptom of a weak
  test: it computed `cz/H0` and then asserted a 20%-wide range instead of
  comparing against it. The assertion now checks agreement with the computed
  value to within 0.5 Mpc (~1%).

## [1.0.0] — 2026-10-15

The first release packaged as a complete, self-contained project. The physics
is that of 0.3.0; what changes here is everything around it — packaging,
documentation, examples and packaging correctness.

### Added

- **Bilingual documentation.** `docs/en/` and `docs/zh/` are parallel trees
  covering installation, a user guide, the HI/H₂ models, the API reference and
  the verification record. `README.md` (English) and `README.zh-CN.md` (中文)
  are cross-linked.
- **Bilingual examples.** `examples/en/` and `examples/zh/` each hold a
  tutorial notebook, a demo script and a batch script. The two tutorials share
  byte-identical code, so they produce identical numbers.
- `docs/build_html.py` — builds a single-file-per-page HTML rendering of the
  documentation from the Markdown sources, using only the standard library.
  A test asserts the committed HTML matches the sources, which removes the
  drift the two hand-maintained doc trees had accumulated.
- `docs/en/verification.md` — the validation record (bit-for-bit agreement with
  the reference implementation, end-to-end checks on real TNG-100 subhalos),
  previously kept outside the repository.
- `examples/*/batch_subhalos.py` — loop over a subhalo catalog, replacing the
  older root-level `example_vodka.ipynb`.
- `LICENSE` (GPLv3) and `NOTICE`. The licence file was previously referenced by
  the README but absent, which for a GPLv3-derived fork is a compliance defect.
  `NOTICE` records the provenance of each part of the package.
- `CITATION.cff`, `CONTRIBUTING.md`, `environment.yml`, `requirements.txt`,
  `requirements-dev.txt`, `MANIFEST.in`, and a GitHub Actions CI workflow
  testing Python 3.9–3.13 plus a core-dependencies-only import check.
- `[tool.pytest.ini_options]` with `pythonpath`, so `pytest` works from a fresh
  clone without an editable install.
- `[tool.setuptools.packages.find]`, so a new top-level directory can never be
  mistaken for a package.
- `examples/outputs/` as the single, git-ignored destination for example
  products, replacing the mix of repo-root and `examples/examplefiles/` paths.

### Changed

- **Installation works.** `pyproject.toml` declared `requires-python = ">=3.10"`
  but omitted `typing_extensions`, which `astrovodka/datacube.py` needs on
  3.9 and 3.10 (`typing.Self` only enters the stdlib in 3.11). On a clean
  Python 3.10 environment `import astrovodka` therefore failed outright, and
  `pip install -e .` was blocked on the Python 3.9 environment the project is
  actually developed on. The floor is now `>=3.9`, matching what is tested, and
  `typing_extensions` is declared as a conditional dependency.
- The `hdf5_output` extra is renamed `hdf5`, and its documentation corrected:
  h5py is required to *read* TNG snapshots, not only to write HDF5 output.
  Both the old and new names are equally valid before PEP 685 normalisation;
  the new name avoids the underscore.
- Dependencies now carry lower bounds (`numpy>=1.20`, `scipy>=1.6`,
  `astropy>=5.0`, `tqdm>=4.40`) instead of being unpinned.
- `pyproject.toml` gained `license`, `authors`, `maintainers`, `keywords`,
  `[project.urls]`, a fuller classifier set, and ruff/mypy/coverage config.
- Example and tutorial paths are no longer hard-coded. TNG data comes from
  `$TNG_BASE_PATH`, outputs go to `examples/outputs/` resolved relative to the
  checkout, and the notebooks locate the package via `$ASTROVODKA_ROOT` or by
  walking up from the working directory, instead of a `sys.path.insert` naming
  one machine.
- `astrovodka.noise` imports `numpy.random.Generator` (public API since
  numpy 1.17) rather than the private `numpy.random._generator`.
- The `pyproject.toml` description no longer says "K13 HI modeling", which
  stopped being accurate when the model set grew to four.

### Removed

- `astrovodka/sources/_k13_model.py` and `tests/test_k13_model.py` (deleted in
  0.3.0; their bytecode and build metadata lingered, and are now gone too).
- `example_vodka.ipynb` — superseded by `examples/*/batch_subhalos.py`.
- Build and run artifacts from the repository: `astrovodka.egg-info/`
  (which still recorded version 0.1.0), `.pytest_cache/` (which still listed
  `test_k13_model` node ids), `__pycache__/` trees containing orphaned
  bytecode for `_k13_model`, an 8 MB `astrovodka_demo.fits` at the repository
  root, the 63 MB `examples/examplefiles/` output directory, the `examplefigure/`
  directory, and the empty tracked file `examples/examplesrequirement.md`.
- Dependency on `six`, which was imported unguarded on the TNG read path but
  never declared. Its single use, `six.string_types`, is just `str` on Python 3.

## [0.3.0] — 2026-10-15

### Changed

- **The HI model was replaced.** `astrovodka/sources/_k13_model.py` is deleted,
  not deprecated, and `astrovodka/himodel/` takes its place with four
  keyword-selectable Diemer+2018 models: `GD14_vol`, `GD14_map`, `K13_vol`
  (default) and `K13_map`. Selection is by keyword on `compute_hi_mass`,
  `TNGSource(hi_model=...)` and `SubhaloDataProducts.from_subhalo(hi_model=...)`.

### Fixed

- **The adaptive FFT UV field was silently the uniform background.** The 3D
  interpolation built its `RegularGridInterpolator` from one axis' edge array
  reused for all three coordinates. Each axis is centred on its own box centre,
  so for almost every particle two of the three fell outside the interpolator's
  bounds and `fill_value=0.0` returned zero. Measured on TNG-100 snapshot 99
  subhalo 122: raw FFT grid maximum 7.05, interpolated maximum 0.0000. The
  consequence was that `uv_mode="adaptive_fft"` and `uv_mode="global"` gave
  numerically identical masses, and molecular masses were wrong by up to 42%.

Because the old `"adaptive_fft"` was numerically identical to `"global"`, the
change of default costs nothing that previously worked.

## [0.2.0] — 2026-10-15

### Fixed

Four defects, each with a regression test:

- **`spectral_unit="MHz"` discarded the signal.** The velocity-to-frequency
  conversion for `channel_width` returned `f_rest` instead of the width, so a
  frequency-axis cube retained 0.007% of the flux. Widths are now converted
  differentially, and MHz and km/s cubes agree to a flux ratio of 1.0000.
- **`D_A` was computed but never used**, so every angular size was too small by
  `(1+z)²`. Sky placement now uses `D_A`; the flux normalisation explicitly
  uses `D_L`.
- **`vsys` used `H(z)·D_C`**, which is not the radio convention of the WCS
  spectral axis, is 74% too large at z = 0.5, exceeds `c` beyond z ≈ 0.79, and
  produced negative frequencies at high redshift. Now `c·z/(1+z)`.
- **`Martini.reset()` dropped `spectral_unit`**, silently flipping a velocity
  cube to a frequency axis. `copy()`, `save_state`/`load_state()` and
  `from_wcs()` had related inconsistencies; all reconstruction paths now
  round-trip the spectral unit.

## [0.1.0] — 2026-07-03

Initial fork of MARTINI, adding cosmological distances via
`astropy.cosmology.Planck15`, the radio-convention systemic velocity, local-disk
TNG access, `observed_redshift`, 3D data products and DM centre-of-mass
centring.

[Unreleased]: https://github.com/LRayleighJ/astrovodka/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/LRayleighJ/astrovodka/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/LRayleighJ/astrovodka/compare/v0.3.0...v1.0.0
[0.3.0]: https://github.com/LRayleighJ/astrovodka/releases/tag/v0.3.0
[0.2.0]: https://github.com/LRayleighJ/astrovodka/releases/tag/v0.2.0
[0.1.0]: https://github.com/LRayleighJ/astrovodka/releases/tag/v0.1.0
