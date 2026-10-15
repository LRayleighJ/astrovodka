# astrovodka.data_products — 三维网格化数据产品

三维网格化数据产品的容器：UV 场，以及共动实空间和红移空间中的 DM/HI/H2 质量密度，带 HDF5 读写。

## SubhaloDataProducts

### `__init__(self, ngrid=128)`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `ngrid` | `int` | 每个维度的网格分辨率（默认：128） |

---

### `from_subhalo(cls, basePath, snapNum, subID, ngrid=128, hi_model='K13_vol', hi_config=None, uv_mode='adaptive_fft', load_subhalo_only=True, cosmology=None)` → `SubhaloDataProducts`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `basePath` | `str` | TNG 数据根目录 |
| `snapNum` | `int` | 快照号 |
| `subID` | `int` | 子晕 ID |
| `ngrid` | `int` | 网格分辨率（默认：128） |
| `hi_model` | `str` | 四种可选模型之一：`'K13_vol'`（默认）、`'K13_map'`、`'GD14_vol'`、`'GD14_map'` |
| `hi_config` | `HIModelConfig \| dict \| None` | HI 模型设置 |
| `uv_mode` | `str` | 仅体密度模型：`'adaptive_fft'`（默认）或 `'global'` |
| `load_subhalo_only` | `bool` | 仅子晕（True）或 FOF 群 |
| `cosmology` | `Cosmology \| None` | 用于红移空间畸变的 H(z) |

> 读取气体（PartType0）+ 暗物质（PartType1），投影模型还需读入恒星粒子。用四种可选模型中选定的那一个计算 HI/H2 质量。`ngrid` 参数同时决定 UV 网格的分辨率。暗物质质量取自 `Header.MassTable[1]`。
> 网格边界以**共动 ckpc/h** 为单位。

---

### `save_hdf5(self, filename)`

保存全部 7 个网格 + 边界 + 元数据。LZF 压缩、分块存储。

### `load_hdf5(cls, filename)` → `SubhaloDataProducts`

恢复全部网格与元数据。

---

## 数据网格属性

所有网格都是 `(ngrid, ngrid, ngrid)` 的 float32 数组：

| 属性 | 单位 | 说明 |
|-----------|------|-------------|
| `uv_field` | Draine 单位 | UV 辐射场 |
| `dm_real` | M<sub>sun</sub>/(ckpc/h)³ | 共动 DM 密度（实空间） |
| `dm_rsd` | M<sub>sun</sub>/(ckpc/h)³ | 共动 DM 密度（红移空间） |
| `hi_real` | M<sub>sun</sub>/(ckpc/h)³ | 共动 HI 密度（实空间） |
| `hi_rsd` | M<sub>sun</sub>/(ckpc/h)³ | 共动 HI 密度（红移空间） |
| `h2_real` | M<sub>sun</sub>/(ckpc/h)³ | 共动 H2 密度（实空间） |
| `h2_rsd` | M<sub>sun</sub>/(ckpc/h)³ | 共动 H2 密度（红移空间） |
| `grid_edges` | ckpc/h | 3 个边界数组组成的元组（长度为 ngrid+1） |
| `metadata` | — | redshift、hubble_h、网格范围等 |

> **红移空间公式：** s<sub>x</sub> = x + v<sub>los</sub> / (a × H(a))
> 其中 Δx = v<sub>los</sub> / (a×H(a)) × 1000 × h [ckpc/h]。
