# astrovodka.group_products — 星系群/团尺度的三维产品

把整个 **FOF group**（它包含的每个 subhalo，加上 group 中不归属任何 subhalo 的
fuzz）网格化，并将气体分为 ISM、CGM、IGM 三类。[data_products](data_products.md)
画的是单个星系，这里画的是它所属的整个 group。

| 名称 | 值 | 含义 |
|------|-----|------|
| `ISM` | `0` | 归属于某个 subhalo，且距其中心 < `k x R_ISM` |
| `CGM` | `1` | 归属于某个 subhalo，但在该半径之外（不做 R<sub>200c</sub> 截断） |
| `IGM` | `2` | 在 group 内但不属于任何 subhalo（即 fuzz） |

`CLASS_LABELS` 把上述整数映射回 `"ISM"`、`"CGM"`、`"IGM"`。

---

## 归属关系

胞元归属来自 group catalog 的 offsets，而不是空间上的猜测：

```
subhalo s 的气体胞元  =  [Subhalo/SnapByType[s, 0], + SubhaloLenType[s, 0])
group 的胞元          =  [Group/SnapByType[g, 0],  + GroupLenType[g, 0])
fuzz                  =  group 的胞元 减去 所有 subhalo 的胞元
```

> **长度必须取自 `SubhaloLenType`，绝不能用 offsets 数组做相邻差分。**
> `SnapByType` 是累积起始偏移表，而一个 group 的*最后一个* subhalo 的下一行已经属于
> 下一个 group，差分会把该 group 的 fuzz 悄悄吞掉——在 TNG-100 快照 99 的 group 59
> 上，最后一个 subhalo 的 `SubhaloLenType == 0`，跨度却是 319,249。
> `from_group` 会把每个块的质量求和与 `SubhaloMassType` 对照，不一致就直接报错；
> 正是这个检查让归属关系可信，而不是想当然。

fuzz 是实打实的成分：在 TNG-100 最大的 group 里占气体质量 1.5%，在 group 59 这个
质量上约占 10%，更小的 group 里可达约 24%。

---

## `ism_radius(rad_gas, rad_star, k_gas=1.0, k_star=2.0)` → `ndarray`

`min(k_gas x R_half,gas, k_star x R_half,star)`，两个半径取自
`SubhaloHalfmassRadType` 的第 0、4 列，单位与输入一致（ckpc/h）。

| 名称 | 类型 | 说明 |
|------|------|------|
| `rad_gas`, `rad_star` | `array-like` | 气体与恒星半质量半径 [ckpc/h] |
| `k_gas`, `k_star` | `float` | 系数（默认 1 与 2） |

> 取 `min` 是让这个定义在 group 尺度上仍然可用的关键。group 的中心 subhalo 因为包含
> 整个热晕，`R_half,gas` 有几百 ckpc/h——TNG-100 group 59 是 263 ckpc/h，而
> `R_half,star` 只有 19 ckpc/h——只用气体半径会把大半个 group 都划进 ISM。
> 非正或非有限的半径视为「不构成限制」；若两者都不可用，结果为 0，该 subhalo 不贡献
> 任何 ISM。

---

## `classify_ism_cgm_igm(owner, r_gal, r_ism, nH=None, T=None, nH_max=None, T_max=None)` → `ndarray[int8]`

| 名称 | 类型 | 说明 |
|------|------|------|
| `owner` | `ndarray of int` | 每个胞元所属的 subhalo 序号，fuzz 为 `-1` |
| `r_gal` | `ndarray` | 到所属 subhalo 中心的距离 [ckpc/h]；`owner < 0` 处忽略 |
| `r_ism` | `ndarray or float` | 每个胞元对应的 ISM 半径（已按 `owner` 索引），或标量 |
| `nH` | `ndarray` | 氢数密度 [cm⁻³]；给出 `nH_max` 时必需 |
| `T` | `ndarray` | 温度 [K]；给出 `T_max` 时必需 |
| `nH_max`, `T_max` | `float` | 可选的相位上限：不满足的胞元被降级出 ISM |

> 这两个上限**默认关闭**，实测就是原因：它们针对的「又密又热」的激波气体
> （`n_H > 0.1 cm^-3` 且 `T > 1e6 K`）在 group 59 中只占 ISM 质量的 0.000，
> 而一个普通的 `n_H < 100` 上限也才剔除 0.6%。

---

## `load_group_members(basePath, snapNum, groupID)` → `dict`

group 内每个 subhalo 的 catalog 记录，取自 subhalo 表中连续的
`[GroupFirstSub, + GroupNsubs)` 行。

返回 `group`（FOF group 自身的记录）、`first`、`pos` (n, 3) [ckpc/h]、
`len_type` (n, 6)、`mass_type` (n, 6) [1e10 M<sub>sun</sub>/h]、`mass`、
`rhalf_type` (n, 6) [ckpc/h]、`sfr`、`flag`，以及 `start` (n, 6)——各 subhalo
各粒子类型在快照数组中的起始偏移。

---

## `GroupDataProducts`

### `__init__(self, ngrid=256)`

### `from_group(cls, basePath, snapNum, groupID, *, ngrid=256, hi_model='K13_vol', hi_config=None, uv_mode='adaptive_fft', compute_hi=True, compute_rsd=False, k_gas=1.0, k_star=2.0, nH_max=None, T_max=None, igm_uv='group', ncpu=1, cosmology=None, box_halfsize_ckpc=None)` → `GroupDataProducts`

| 名称 | 类型 | 说明 |
|------|------|------|
| `basePath` | `str` | TNG 输出目录（含 `groups_<snapNum>/` 与 `postprocessing/offsets`） |
| `snapNum` | `int` | 快照号 |
| `groupID` | `int` | FOF group ID（即 `haloID`） |
| `ngrid` | `int` | 网格边长（默认 256；最大的 group 用 128） |
| `hi_model` | `str` | `'K13_vol'`（默认）或 `'GD14_vol'`；投影模型会被拒绝 |
| `hi_config` | `HIModelConfig \| dict \| None` | 模型设置；`h`、`z` 与盒长始终取自快照 |
| `uv_mode` | `str` | `'adaptive_fft'`（逐 subhalo）或 `'global'` |
| `compute_hi` | `bool` | 为 `False` 时完全跳过 HI/H2 模型 |
| `compute_rsd` | `bool` | 是否同时生成红移空间网格（默认关闭，原因见下） |
| `k_gas`, `k_star` | `float` | ISM 半径系数 |
| `nH_max`, `T_max` | `float \| None` | 可选的 ISM 相位上限 |
| `igm_uv` | `str` | fuzz 的 UV 场：`'group'`（默认）或 `'uvb'` |
| `ncpu` | `int` | 逐 subhalo HI 循环的进程数（需要 `fork`） |
| `cosmology` | `Cosmology \| None` | 用于 RSD 位移；默认 Planck15 |
| `box_halfsize_ckpc` | `float \| None` | 指定盒子半宽；否则取粒子范围的 1.05 倍 |

> **HI 是逐 subhalo 计算的**，每个都在自己的自适应 FFT 盒里——这正是体密度模型被
> 标定时所用的方案。对整个 group 只做一次 FFT 会把 UV 场抹到星系团尺度的约 30 kpc
> 网格上。只有含气体的 subhalo 需要计算（group 59 的 1710 个成员里只有 68 个），
> 所以循环远比 subhalo 总数便宜。fuzz 单独建模，用整组的 UV 场。
>
> **红移空间默认不生成。** 在 z=0 时 1 km/s 的视向速度对应 10 ckpc/h，因此成员速度
> 跨越几百 km/s 的 group，在红移空间中沿视线方向能延伸好几个 Mpc：group 59 需要
> 17 cMpc 的轴长，而天球面只有 1.9 cMpc。位移是在减去 group 整体运动之后施加的，
> 这也是生成模拟观测所用的参考系；若不做这个扣除，图会是空的
> （`SubhaloDataProducts` 的红移空间网格正是因此丢掉了约 99.5% 的质量）。

### 先在实空间分类

**分类在实空间中完成，之后才带入红移空间，绝不在红移空间中判定。**
`classify_ism_cgm_igm` 作用在实空间坐标上，红移空间网格则是把同一批胞元沿视线方向
平移后用它们的实空间 `class_id` 作掩膜得到。由此得到的三点结论都已在 TNG-100 group 59
上验证：

- `cells`、各类质量以及**每一个实空间网格，无论 `compute_rsd` 开或关都逐位相同**——
  实空间产品是 canonical 的。
- 每一类在变换中严格守恒：实空间与红移空间网格的质量之比对 ISM、CGM、IGM 均为
  1.000000，且都与 metadata 中的总量一致。
- 红移空间有**自己的网格盒子**，视线方向更长（`grid_edges_rsd`，metadata 中另有
  `box_halfsize_ckpc_rsd` 与 `box_bounds_ckpc_h_rsd`），因此开启 RSD 不会降低实空间
  的分辨率。两套盒子的格子体积不同，`f_mass_in_rsd_box_gas` 已计入这一点。

---

### `save_hdf5(self, filename)` / `load_hdf5(cls, filename)`

保存网格（LZF 压缩、分块）、逐胞元数组、subhalo 表与 metadata。
值为 `None` 的 metadata 会写成字符串 `"none"`；标量的 subhalo 条目读回时是长度为 1
的数组——HDF5 的标量数据集无法使用文件其余部分所用的压缩。

---

## 网格属性

所有网格都是 `(ngrid, ngrid, ngrid)` 的 float32，共用同一个盒子，`x` 沿视线方向，
单位 M<sub>sun</sub>/(ckpc/h)³：

| 属性 | 说明 |
|-----------|-------------|
| `dm_real`, `dm_rsd` | 暗物质总量 |
| `dm_igm_real` | 不属于任何 subhalo 的暗物质 |
| `gas_real`, `gas_rsd` | 气体总量 |
| `gas_ism_real`, `gas_cgm_real`, `gas_igm_real` | 按类别分的气体 |
| `hi_real`, `hi_rsd` | 中性氢 |
| `hi_ism_real`, `hi_ism_rsd`, `hi_cgm_real`, `hi_cgm_rsd`, `hi_igm_real`, `hi_igm_rsd` | 按类别分的 HI |
| `h2_real` | 分子氢 |
| `t_mw`, `nh_mw` | 气体质量加权的温度 [K] 与 n<sub>H</sub> [cm⁻³] |
| `grid_edges` | 实空间盒子的 3 个边界数组 [ckpc/h] |
| `grid_edges_rsd` | 红移空间盒子的边界数组（视线方向更长）；`compute_rsd` 关闭时为 `None` |

## `cells` — 逐胞元数组

`coords_ckpc_h`、`m_gas_msun`、`m_hi_msun`、`m_h2_msun`、`nH_cm3`、`T_K`、
`r_gal_ckpc_h`、`r_grp_over_r200c`、`owner`、`class_id`。

> 保留它们是为了在重新分类、做径向剖面或局部放大重网格时不必再读快照——group 尺度的
> 网格分辨不出星系，所以那些缩放图是直接对小体积内的胞元重新分箱得到的。

## `metadata`

group 身份（`groupID`、`snapNum`、`redshift`、`hubble_h`、`M200c_msun`、
`R200c_ckpc_h`、`R200c_kpc`、`group_pos_ckpc_h`、`Nsubs`、`N_sub_with_gas`）、
计数（`n_gas`、`n_dm`、`n_cells_ism/cgm/igm`）、质量收支（`m_gas_*_msun`、
`m_hi_*_msun`、`m_h2_total_msun`、`m_dm_total_msun`、`m_dm_igm_total_msun`）、
分类设置（`ism_k_gas`、`ism_k_star`、`ism_nH_max`、`ism_T_max`）、盒子
（`box_halfsize_ckpc`、`box_bounds_ckpc_h`）、模型（`hi_model`、`uv_mode`、
`igm_uv`、`vol_fine_grid`），以及开启 RSD 时的 `f_mass_in_rsd_box_gas`。

## `subhalos`

即 `load_group_members` 返回的表，便于把逐 subhalo 的量与网格对应起来，而不必再读
catalog。
