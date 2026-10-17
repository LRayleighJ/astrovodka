# astrovodka.sources._illustris_tools — TNG 数据 I/O

读取 IllustrisTNG 模拟数据的内部工具，转发自 `martini.sources._illustris_tools`（其本身源自 `illustris_python` 工具包）。

### `loadHeader(basePath, snapNum)` → `dict`

读取群表（group catalog）头。返回的字典包含：`HubbleParam`、`Time`、`Redshift`、`BoxSize` 等。

### `loadSingle(basePath, snapNum, haloID=-1, subhaloID=-1)` → `dict`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `haloID` / `subhaloID` | `int` | 二者互斥 |

读取单个晕/子晕的全部群表字段。键包括：`SubhaloPos`、`SubhaloVel`、`SubhaloGrNr`、`SubhaloLenType`。

### `loadSubset(basePath, snapNum, partType, *, fields, subset=None, mdi=None, sq=True, float32=False)` → `dict | ndarray`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `partType` | `str \| int` | `'gas'`(0)、`'dm'`(1)、`'star'`(4)、`'bh'`(5) |
| `fields` | `list[str]` | 要读取的字段名 |
| `subset` | `dict \| None` | 来自 `getSnapOffsets` |
| `mdi` | `list[int\|None] \| None` | 每个字段的多维索引 |
| `sq` | `bool` | 单字段时 → 返回数组而不是字典 |
| `float32` | `bool` | 把 float64 降为 float32 |

自动处理由多个文件分块的快照。

### `getSnapOffsets(basePath, snapNum, id, type)` → `dict`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `id` | `int` | 晕或子晕 ID |
| `type` | `str` | `'Group'` 或 `'Subhalo'` |

返回：`lenType`、`offsetType`、`snapOffsets`。

### `partTypeNum(partType: int | str)` → `int`

映射关系：gas/cells→0，dm→1，star→4，bh→5，等等。

### `snapPath(basePath, snapNum, chunkNum=0)` → `str`

快照 HDF5 文件的路径。

### `gcPath(basePath, snapNum, chunkNum=0)` → `str`

群表 HDF5 文件的路径。

### `offsetPath(basePath, snapNum)` → `str`

偏移量 HDF5 文件的路径。

### `getNumPart(header: dict)` → `ndarray`

每种粒子类型的粒子数（共 6 种类型）。
