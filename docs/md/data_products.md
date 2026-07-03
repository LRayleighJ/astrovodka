# astrovodka.data_products — 3D Gridded Data Products

Container for 3D gridded data products: UV field, DM/HI/H2 mass density in comoving real-space and redshift-space, with HDF5 I/O.

## SubhaloDataProducts

### `__init__(self, ngrid=128)`

| Name | Type | Description |
|------|------|-------------|
| `ngrid` | `int` | Grid resolution per dimension (default: 128) |

---

### `from_subhalo(cls, basePath, snapNum, subID, ngrid=128, uv_mode='global', load_subhalo_only=True, cosmology=None)` → `SubhaloDataProducts`

| Name | Type | Description |
|------|------|-------------|
| `basePath` | `str` | TNG data root |
| `snapNum` | `int` | Snapshot number |
| `subID` | `int` | Subhalo ID |
| `ngrid` | `int` | Grid resolution (default: 128) |
| `uv_mode` | `str` | `'global'` or `'adaptive_fft'` |
| `load_subhalo_only` | `bool` | Subhalo-only (True) or FOF group |
| `cosmology` | `Cosmology \| None` | For RSD H(z) |

> Loads gas (PartType0) + DM (PartType1). Computes K13 HI/H2 masses. DM mass from `Header.MassTable[1]`.
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
| `grid_edges` | ckpc/h | Tuple of 3 edge arrays (length ngrid+1) |
| `metadata` | — | redshift, hubble_h, grid bounds, etc. |

> **Redshift-space formula:** s<sub>x</sub> = x + v<sub>los</sub> / (a × H(a))
> where Δx = v<sub>los</sub> / (a×H(a)) × 1000 × h [ckpc/h].
