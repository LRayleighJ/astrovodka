# astrovodka.sources.tng_source — TNGSource

从本地磁盘读取 IllustrisTNG 数据，并用四种 [Diemer+2018 模型](himodel.md)之一
逐粒子计算 HI/H2 质量。继承自 `SPHSource`。

## TNGSource(SPHSource)

### `__init__(self, basePath, snapNum, subID, *, cosmology=None, hi_model='K13_vol', hi_config=None, uv_mode='adaptive_fft', observed_redshift=None, load_subhalo_only=True, vpeculiar=0 km/s, rotation=None, L_coords=None, ra=0 deg, dec=0 deg, coordinate_frame=ICRS())`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `basePath` | `str` | TNG 数据根目录。例如 `/data/.../TNG/TNG-100/` |
| `snapNum` | `int` | 快照号：99 (z~0)、67 (z~0.5)、50 (z~1.0) |
| `subID` | `int` | 子晕 ID |
| `cosmology` | `Cosmology \| None` | 默认：Planck15 |
| `hi_model` | `str` | 四种可选模型之一：`'K13_vol'`（默认）、`'K13_map'`、`'GD14_vol'`、`'GD14_map'`；不区分大小写，`'K13'`/`'GD14'` 表示体密度变体 |
| `hi_config` | `HIModelConfig \| dict \| None` | HI 模型设置；`h`/`z` 始终取自快照 |
| `uv_mode` | `str` | 仅体密度模型：`'adaptive_fft'`（默认）或 `'global'` |
| `observed_redshift` | `float \| None` | 用观测红移覆盖距离计算。K13 使用快照红移 z |
| `load_subhalo_only` | `bool` | **True**（快，仅子晕）或 False（FOF 群） |
| `vpeculiar` | `Quantity[km/s]` | 本动速度（默认：0） |
| `rotation` | `Rotation \| None` | 源的取向 |
| `L_coords` | `L_coords \| None` | 角动量对齐 |
| `ra` / `dec` | `Quantity[deg]` | 源的中心（默认：0） |
| `coordinate_frame` | `BaseRADecFrame` | 坐标系（默认：ICRS） |

> **说明：**
> - 前 3 个参数为位置参数，其余为仅关键字参数。
> - **`observed_redshift`**：把 HI 模型物理（快照 z）与距离（观测 z）解耦。
> - **`load_subhalo_only=True`** 是默认值，比加载 FOF 群快得多。
> - 体密度模型使用 8 个 TNG 气体场：Masses、Density、InternalEnergy、SFR、NeutralHydrogenAbundance、GFM_Metals、GFM_Metallicity、Coordinates。投影（`*_map`）模型不需要 InternalEnergy，但需要恒星粒子和子晕记录。
> - `hi_model` 会在从磁盘读取任何数据之前先校验，因此拼错也不会付出代价。

## 关键属性

| 属性 | 类型 | 说明 |
|-----------|------|-------------|
| `source_redshift` | `float` | 快照红移 |
| `D_L` | `Quantity[Mpc]` | 光度距离 |
| `D_A` | `Quantity[Mpc]` | 角直径距离 |
| `m_H2_g` | `Quantity[M_sun]` | 每个粒子的 H2 质量 |
| `hi_model` | `str` | 解析后的模型名 |
| `hi_config` | `HIModelConfig` | 实际使用的配置 |
| `f_mol_g`, `f_neutral_g` | `ndarray` | 逐单元的分子/中性比例 |
| `u_mw_g` | `ndarray` | 逐单元的 UV 场（体密度模型） |
| `hi_diagnostics` | `dict` | 投影几何、对齐情况、`f_gas_in_map` |
