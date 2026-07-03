# astrovodka API Documentation

Version 0.2.0 — Synthetic HI observations with cosmological distances and K13 HI modeling.

## Modules

| Module | Description |
|--------|-------------|
| [martini](martini.md) | Core pipeline — Martini and GlobalProfile classes |
| [datacube](datacube.md) | DataCube — WCS grid definition |
| [sph_source](sph_source.md) | SPHSource — Generic SPH source (redshift + distance modes) |
| [tng_source](tng_source.md) | TNGSource — IllustrisTNG + K13 HI model |
| [beams](beams.md) | GaussianBeam — Telescope beam models |
| [noise](noise.md) | GaussianNoise — Noise models |
| [spectral_models](spectral_models.md) | GaussianSpectrum, DiracDeltaSpectrum — 21-cm line models |
| [sph_kernels](sph_kernels.md) | SPH smoothing kernels (6 types with adaptive fallback) |
| [cosmology_utils](cosmology_utils.md) | Cosmological distance and unit conversion helpers |
| [_k13_model](_k13_model.md) | K13 volumetric HI/H2 model functions |
| [data_products](data_products.md) | SubhaloDataProducts — 3D gridded UV/DM/HI/H2 |
| [L_coords](L_coords.md) | Angular momentum alignment coordinates |
| [_illustris_tools](_illustris_tools.md) | TNG snapshot data I/O utilities |
