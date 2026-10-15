# Changelog

All notable changes to astrovodka are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Nothing yet.

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

[Unreleased]: https://github.com/LRayleighJ/astrovodka/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/LRayleighJ/astrovodka/compare/v0.3.0...v1.0.0
[0.3.0]: https://github.com/LRayleighJ/astrovodka/releases/tag/v0.3.0
[0.2.0]: https://github.com/LRayleighJ/astrovodka/releases/tag/v0.2.0
[0.1.0]: https://github.com/LRayleighJ/astrovodka/releases/tag/v0.1.0
