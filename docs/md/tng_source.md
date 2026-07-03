# astrovodka.sources.tng_source — TNGSource

Loads IllustrisTNG data from local disk and computes per-particle HI/H2 masses using the K13 (Krumholz 2013) volumetric model. Inherits from `SPHSource`.

## TNGSource(SPHSource)

### `__init__(self, basePath, snapNum, subID, *, cosmology=None, uv_mode='global', observed_redshift=None, load_subhalo_only=True, vpeculiar=0 km/s, rotation=None, L_coords=None, ra=0 deg, dec=0 deg, coordinate_frame=ICRS())`

| Name | Type | Description |
|------|------|-------------|
| `basePath` | `str` | TNG data root. E.g. `/data/.../TNG/TNG-100/` |
| `snapNum` | `int` | Snapshot: 99 (z~0), 67 (z~0.5), 50 (z~1.0) |
| `subID` | `int` | Subhalo ID |
| `cosmology` | `Cosmology \| None` | Default: Planck15 |
| `uv_mode` | `str` | `'global'` (UVB, fast) or `'adaptive_fft'` (FFT, accurate) |
| `observed_redshift` | `float \| None` | Override observed redshift for distance. K13 uses snapshot z |
| `load_subhalo_only` | `bool` | **True** (fast, subhalo-only) or False (FOF group) |
| `vpeculiar` | `Quantity[km/s]` | Peculiar velocity (default: 0) |
| `rotation` | `Rotation \| None` | Source orientation |
| `L_coords` | `L_coords \| None` | Angular momentum alignment |
| `ra` / `dec` | `Quantity[deg]` | Source centroid (default: 0) |
| `coordinate_frame` | `BaseRADecFrame` | Frame (default: ICRS) |

> **Notes:**
> - First 3 params positional; rest keyword-only.
> - **`observed_redshift`**: decouples K13 physics (snapshot z) from distance (observed z).
> - **`load_subhalo_only=True`** is the default and much faster than FOF group loading.
> - K13 uses 8 TNG fields: Masses, Density, InternalEnergy, SFR, NeutralHydrogenAbundance, GFM_Metals, GFM_Metallicity, Coordinates.

## Key Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `source_redshift` | `float` | Snapshot redshift |
| `D_L` | `Quantity[Mpc]` | Luminosity distance |
| `D_A` | `Quantity[Mpc]` | Angular diameter distance |
| `m_H2_g` | `Quantity[M_sun]` | H2 mass per particle (K13) |
