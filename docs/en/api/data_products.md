# astrovodka.data_products — 3D Gridded Data Products

Container for 3D gridded data products: UV field, DM/HI/H2 mass density in comoving real-space and redshift-space, with HDF5 I/O.

## SubhaloDataProducts

### `__init__(self, ngrid=128)`

| Name | Type | Description |
|------|------|-------------|
| `ngrid` | `int` | Grid resolution per dimension (default: 128) |

---

### `from_subhalo(cls, basePath, snapNum, subID, ngrid=128, hi_model='K13_vol', hi_config=None, uv_mode='adaptive_fft', load_subhalo_only=True, cosmology=None)` → `SubhaloDataProducts`

| Name | Type | Description |
|------|------|-------------|
| `basePath` | `str` | TNG data root |
| `snapNum` | `int` | Snapshot number |
| `subID` | `int` | Subhalo ID |
| `ngrid` | `int` | Grid resolution (default: 128) |
| `hi_model` | `str` | `'K13_vol'` (default), `'K13_map'`, `'GD14_vol'`, `'GD14_map'` |
| `hi_config` | `HIModelConfig \| dict \| None` | HI model settings |
| `uv_mode` | `str` | Volume models only: `'adaptive_fft'` (default) or `'global'` |
| `load_subhalo_only` | `bool` | Subhalo-only (True) or FOF group |
| `cosmology` | `Cosmology \| None` | For RSD H(z) |

> Loads gas (PartType0) + DM (PartType1), plus stars for the projected models. Computes HI/H2 masses with the selected model. The `ngrid` argument also sets the UV grid resolution. DM mass from `Header.MassTable[1]`.
> Grid edges in **comoving ckpc/h**.

---

### `save_hdf5(self, filename)`

Saves all 7 grids + edges + metadata. LZF compression, chunked storage.

### `load_hdf5(cls, filename)` → `SubhaloDataProducts`

Restores all grids + metadata.

---

## Data Grid Attributes

All grids are `(ngrid, ngrid, ngrid)` float32 arrays:

| Attribute | Unit | Description |
|-----------|------|-------------|
| `uv_field` | Draine units | UV radiation field |
| `dm_real` | M<sub>sun</sub>/(ckpc/h)³ | Comoving DM density (real space) |
| `dm_rsd` | M<sub>sun</sub>/(ckpc/h)³ | Comoving DM density (redshift space) |
| `hi_real` | M<sub>sun</sub>/(ckpc/h)³ | Comoving HI density (real space) |
| `hi_rsd` | M<sub>sun</sub>/(ckpc/h)³ | Comoving HI density (redshift space) |
| `h2_real` | M<sub>sun</sub>/(ckpc/h)³ | Comoving H2 density (real space) |
| `h2_rsd` | M<sub>sun</sub>/(ckpc/h)³ | Comoving H2 density (redshift space) |
| `grid_edges` | ckpc/h | Tuple of 3 edge arrays (length ngrid+1), **real space** |
| `grid_edges_rsd` | ckpc/h | The same for the **redshift-space** box, which is longer along the line of sight |
| `metadata` | — | redshift, hubble_h, grid bounds, `f_mass_in_rsd_box_*`, etc. |

> **Redshift-space formula:** s<sub>x</sub> = x + v<sub>los</sub> / (a × H(a))
> where Δx = v<sub>los</sub> / (a×H(a)) × 1000 × h [ckpc/h], with v<sub>los</sub> measured
> **relative to the DM centre of mass** — the same reference the box is centred on
> (`metadata['rsd_velocity_frame']`).

> **The two boxes are different, so `*_rsd` and `*_real` cannot be subtracted.**
> At z = 0 one km/s is 10 ckpc/h, so the redshift-space extent along the line of
> sight is set by the velocity range, not by the geometry: for subhalo 113 it is
> 2722 ckpc/h against 41 ckpc/h in real space. Δx is a *displacement*, so a
> galaxy-sized box cannot hold it without stretching the real-space grids too —
> which would be a redshift-space choice corrupting real-space data.
> Compare the two grids through their line-of-sight profiles, or use
> `f_mass_in_rsd_box_dm` / `_hi` / `_h2` (1.000 when the box is sized to fit) to
> check nothing was lost.
>
> Before 1.0.0 the shift used *absolute* velocities and kept the real-space box,
> so the `*_rsd` grids retained **0.5%** of their mass. Anything that assumed
> `dm_rsd − dm_real` was therefore comparing a real map against almost nothing.
