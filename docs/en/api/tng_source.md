# astrovodka.sources.tng_source — TNGSource

Loads IllustrisTNG data from local disk and computes per-particle HI/H2 masses using
one of four [Diemer+2018 models](himodel.md). Inherits from `SPHSource`.

## TNGSource(SPHSource)

### `__init__(self, basePath, snapNum, subID, *, cosmology=None, hi_model='K13_vol', hi_config=None, uv_mode='adaptive_fft', observed_redshift=None, load_subhalo_only=True, vpeculiar=0 km/s, rotation=None, L_coords=None, ra=0 deg, dec=0 deg, coordinate_frame=ICRS())`

| Name | Type | Description |
|------|------|-------------|
| `basePath` | `str` | TNG data root. E.g. `/data/.../TNG/TNG-100/` |
| `snapNum` | `int` | Snapshot: 99 (z~0), 67 (z~0.5), 50 (z~1.0) |
| `subID` | `int` | Subhalo ID |
| `cosmology` | `Cosmology \| None` | Default: Planck15 |
| `hi_model` | `str` | `'K13_vol'` (default), `'K13_map'`, `'GD14_vol'`, `'GD14_map'`; case-insensitive, `'K13'`/`'GD14'` mean the volume variant |
| `hi_config` | `HIModelConfig \| dict \| None` | HI model settings; `h`/`z` always come from the snapshot |
| `uv_mode` | `str` | Volume models only: `'adaptive_fft'` (default) or `'global'` |
| `observed_redshift` | `float \| None` | Override observed redshift for distance. K13 uses snapshot z |
| `load_subhalo_only` | `bool` | **True** (fast, subhalo-only) or False (FOF group) |
| `vpeculiar` | `Quantity[km/s]` | Peculiar velocity (default: 0) |
| `rotation` | `Rotation \| None` | Source orientation |
| `L_coords` | `L_coords \| None` | Angular momentum alignment |
| `ra` / `dec` | `Quantity[deg]` | Source centroid (default: 0) |
| `coordinate_frame` | `BaseRADecFrame` | Frame (default: ICRS) |

> **Notes:**
> - First 3 params positional; rest keyword-only.
> - **`observed_redshift`**: decouples the HI-model physics (snapshot z) from distance (observed z).
> - **`load_subhalo_only=True`** is the default and much faster than FOF group loading.
> - The volume models use 8 TNG gas fields: Masses, Density, InternalEnergy, SFR, NeutralHydrogenAbundance, GFM_Metals, GFM_Metallicity, Coordinates. The projected (`*_map`) models do not need InternalEnergy but do need the stellar particles and the subhalo record.
> - `hi_model` is validated before anything is read from disk, so a typo costs nothing.

## Key Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `source_redshift` | `float` | Snapshot redshift |
| `D_L` | `Quantity[Mpc]` | Luminosity distance |
| `D_A` | `Quantity[Mpc]` | Angular diameter distance |
| `m_H2_g` | `Quantity[M_sun]` | H2 mass per particle |
| `hi_model` | `str` | The resolved model name |
| `hi_config` | `HIModelConfig` | The configuration actually used |
| `f_mol_g`, `f_neutral_g` | `ndarray` | Per-cell molecular / neutral fractions |
| `u_mw_g` | `ndarray` | Per-cell UV field (volume models) |
| `hi_diagnostics` | `dict` | Projection geometry, alignment, `f_gas_in_map` |
