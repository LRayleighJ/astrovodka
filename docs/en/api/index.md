# astrovodka API reference

Version 1.0.0 — mock resolved HI observations of IllustrisTNG galaxies, with
cosmological distances and four selectable Diemer+2018 HI/H2 models.

This is the per-module reference. For the concepts behind it, start with the
[user guide](../guide.md); for the HI/H2 physics, see
[the HI/H2 models](../hi-models.md); for installation, see
[installation and environment setup](../installation.md).

## Modules

| Module | Description |
|--------|-------------|
| [martini](martini.md) | Core pipeline — Martini and GlobalProfile classes |
| [datacube](datacube.md) | DataCube — WCS grid definition |
| [sph_source](sph_source.md) | SPHSource — Generic SPH source (redshift + distance modes) |
| [tng_source](tng_source.md) | TNGSource — IllustrisTNG + selectable HI model |
| [beams](beams.md) | GaussianBeam — Telescope beam models |
| [noise](noise.md) | GaussianNoise — Noise models |
| [spectral_models](spectral_models.md) | GaussianSpectrum, DiracDeltaSpectrum — 21-cm line models |
| [sph_kernels](sph_kernels.md) | SPH smoothing kernels (6 types with adaptive fallback) |
| [cosmology_utils](cosmology_utils.md) | Cosmological distance and unit conversion helpers |
| [himodel](himodel.md) | Four selectable HI/H2 models (GD14/K13, volume/projected) |
| [data_products](data_products.md) | SubhaloDataProducts — 3D gridded UV/DM/HI/H2 |
| [group_products](group_products.md) | GroupDataProducts — FOF-group gas/DM with ISM/CGM/IGM |
| [L_coords](L_coords.md) | Angular momentum alignment coordinates |
| [_illustris_tools](_illustris_tools.md) | TNG snapshot data I/O utilities |
