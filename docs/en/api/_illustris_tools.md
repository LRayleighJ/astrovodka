# astrovodka.sources._illustris_tools — TNG Data I/O

Internal utilities for reading IllustrisTNG simulation data, re-exported from `martini.sources._illustris_tools` (itself derived from the `illustris_python` toolkit).

### `loadHeader(basePath, snapNum)` → `dict`

Loads group catalog header. Returns dict with: `HubbleParam`, `Time`, `Redshift`, `BoxSize`, etc.

### `loadSingle(basePath, snapNum, haloID=-1, subhaloID=-1)` → `dict`

| Name | Type | Description |
|------|------|-------------|
| `haloID` / `subhaloID` | `int` | Mutually exclusive |

Loads all group catalog fields for one halo/subhalo. Keys: `SubhaloPos`, `SubhaloVel`, `SubhaloGrNr`, `SubhaloLenType`.

### `loadSubset(basePath, snapNum, partType, *, fields, subset=None, mdi=None, sq=True, float32=False)` → `dict | ndarray`

| Name | Type | Description |
|------|------|-------------|
| `partType` | `str \| int` | `'gas'`(0), `'dm'`(1), `'star'`(4), `'bh'`(5) |
| `fields` | `list[str]` | Field names to load |
| `subset` | `dict \| None` | From `getSnapOffsets` |
| `mdi` | `list[int\|None] \| None` | Multi-dim index per field |
| `sq` | `bool` | Single field → return array instead of dict |
| `float32` | `bool` | Downcast float64 to float32 |

Handles multi-file snapshot chunks automatically.

### `getSnapOffsets(basePath, snapNum, id, type)` → `dict`

| Name | Type | Description |
|------|------|-------------|
| `id` | `int` | Halo or subhalo ID |
| `type` | `str` | `'Group'` or `'Subhalo'` |

Returns: `lenType`, `offsetType`, `snapOffsets`.

### `partTypeNum(partType: int | str)` → `int`

Maps: gas/cells→0, dm→1, star→4, bh→5, etc.

### `snapPath(basePath, snapNum, chunkNum=0)` → `str`

Path to snapshot HDF5 file.

### `gcPath(basePath, snapNum, chunkNum=0)` → `str`

Path to group catalog HDF5 file.

### `offsetPath(basePath, snapNum)` → `str`

Path to offsets HDF5 file.

### `getNumPart(header: dict)` → `ndarray`

Number of particles per type (6 types).
