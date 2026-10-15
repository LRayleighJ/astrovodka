# `astrovodka.himodel` — HI/H2 质量模型

用于给出 IllustrisTNG 气体单元中性氢与分子氢含量的四种可选模型，
全部按照
[Diemer+2018 (arXiv:1806.02341)](https://arxiv.org/abs/1806.02341) 的实现方式：

| 方法 | 面密度来源 | 参考文献 |
|--------|----------------------|-----------|
| `GD14_vol` | 逐单元 Jeans 近似，`Sigma = lambda_J * rho` | Gnedin & Draine (2014)，附录 C.3 |
| `GD14_map` | 面向观测者的 128² 投影 | 同上 + §2.3 |
| `K13_vol` | 逐单元 Jeans 近似 | Krumholz, McKee & Tumlinson (2013)，附录 C.4 |
| `K13_map` | 面向观测者的 128² 投影 | 同上 + §2.3 |

这四个模型都取代了 MARTINI 基于压强的划分方案，后者是在
近邻星系上标定的，且不随红移变化。

## 快速开始

```python
from astrovodka.himodel import HIModelConfig, compute_hi_mass

cfg = HIModelConfig(h=0.6774, z=0.0)        # h、z 取自快照头

res = compute_hi_mass(gas, model='K13_vol', cfg=cfg)
res.m_hi          # 每个气体单元，TNG 代码单位 (1e10 Msun/h)
res.m_h2
res.m_hi_msun     # 总量，Msun

# 投影模型还需要恒星粒子和子晕记录
res = compute_hi_mass(gas, model='K13_map', cfg=cfg, stars=stars, sub=sub)
```

## 用关键字选择模型

`model` 不区分大小写。只写模型名 `'K13'` 和 `'GD14'` 时表示
体密度变体，也接受双参数写法 `('GD14', 'map')`。

```python
from astrovodka.himodel import METHODS, resolve, is_known

METHODS           # ('GD14_vol', 'GD14_map', 'K13_vol', 'K13_map')
is_known('k13_map')          # True
resolve('K13')['label']      # 'K13_vol'
```

源类上也可以使用同一个关键字：

```python
from astrovodka.sources import TNGSource

source = TNGSource(basePath=..., snapNum=99, subID=113, hi_model='GD14_map')
source.m_HI_total     # 单位为 Msun
source.m_H2_g         # 逐粒子 H2 质量
source.f_mol_g        # 逐单元分子比例
```

## 函数

| 函数 | 说明 |
|----------|-------------|
| `compute_hi_mass(gas, model, cfg, *, stars, sub, uv_mode, uv)` | 在四种模型中任选一种，计算逐单元的 HI 与 H2 质量。主入口。 |
| `compute_galaxy_hi_h2_map(gas, stars, sub, cfg, model, ...)` | 单独使用投影路径：构建面向观测者的投影图，返回星系总量以及对齐诊断信息。 |
| `calc_UMW_global_Faucher11(redshift)` | 以 Draine 单位表示的 Faucher-Giguère+2011 UV 背景。 |
| `resolve(method)` / `is_known(method)` | 方法名查找。 |

如果需要直接使用更底层的组件：

| 模块 | 内容 |
|--------|----------|
| `himodel.physics` | 单位换算、`f_mol_GD14_vol`、`f_mol_K13_vol`、`get_f_neutral_H` |
| `himodel.models` | `f_mol_GD14_map`、`f_mol_K13_map`、`calc_rho_star_map`、`calc_rho_dm_nfw_map`、`MODEL_SPEC` |
| `himodel.uvfield` | `compute_galaxy_fft_uv`、`compute_galaxy_fft_uv_map`、`compute_align_vector`、`deposit_gaussian_maps` |
| `himodel.constants` | 物理常数与模型默认值（不含 `h`） |
| `himodel.config` | `HIModelConfig` |

## 单位

| 量 | 单位 |
|----------|------|
| `HiMassResult.m_hi`、`.m_h2` | TNG 代码单位 `1e10 Msun / h` — 与 `Masses` 字段相同的单位 |
| `HiMassResult.m_hi_msun`、`.m_h2_msun`、`.m_neutral_msun` | Msun |
| `HiMassResult.u_mw` | Draine 场单位（`U_MW = 1` 即银河系的值） |

模型内部的密度用 `HIModelConfig.h` 和 `.z` 换算；
`h` 是逐快照的量，因此刻意在任何地方都
不作为模块常量出现。

## 结果

`compute_hi_mass` 返回一个 `HiMassResult`：

| 属性 | 说明 |
|-----------|-------------|
| `m_hi`, `m_h2` | 逐单元质量，**长度完整** — 被模型排除的单元取零，而不会返回更短的数组 |
| `m_neutral` | `m_hi + m_h2`，恰好就是中性氢质量 |
| `f_mol`, `f_neutral` | 逐单元的比例 |
| `u_mw` | 逐单元 UV 场（体密度模型） |
| `diagnostics` | UV 模式、投影几何、对齐情况、`f_gas_in_map` |
| `cells` | 底层的 `Cells` 表 |

## 投影模型与落在框外的气体

投影模型会构建一张面向观测者的投影图，边长为
`L = 2 * max(2 R_half,gas, 6 R_half,*, 30 kpc)`，并从中读出每个单元的分子
比例。落在该框外的单元得到**零** HI 和 H2：模型的
“星系”*就是*投影框本身，因此该方法在星表中的总量
就是框内总量。

落在框内的气体质量比例会报告为 `f_gas_in_map`；
当它低于 0.9 时 `TNGSource` 会打印警告。若要改为保留框外的气体，
可以把它夹到边缘像素上：

```python
cfg = HIModelConfig(h=0.6774, z=0.0, map_clip_outside=True)
```

这样能守恒质量，代价是把气体堆在边界上。

## 配置

`HIModelConfig` 是一个冻结（frozen）数据类；`HIModelConfig.build(dict_or_config, h=...,
z=...)` 会把用户设置与快照头中的取值对照解析，其中 `h`、`z`
和盒子尺寸始终取自快照。

| 字段 | 默认值 | 含义 |
|-------|---------|---------|
| `ngrid` | 128 | 三维 FFT UV 网格的边长 |
| `f_esc` | 0.1 | UV 逃逸分数（Diemer+2018 §4.2） |
| `box_expand` | 1.1 | 自适应 UV 盒子 = 该值 × 气体坐标范围 |
| `ngrid_map` | 128 | 二维投影网格的边长 |
| `vol_fine_grid` | 0.0 | 第二层更高分辨率的 UV 盒子（论文附录 A.2）；默认关闭 |
| `r_soft_map_mode` | `'dx'` | 二维核的软化长度取自何处 |
| `map_clip_outside` | `False` | 保留投影框外的气体（质量守恒） |
| `k13_variant` | `'2p'` | 投影 K13 模型中由哪种方案给出 `n_CNM` |
| `rho_dm_mode` | `'none'` | 投影 K13 所需的暗物质中面密度 |

完整列表以及每个默认值背后的理由见该类的文档字符串。

## `uv_mode`

体密度模型的 UV 场有两个来源：对恒星形成分布的逐星系 FFT 卷积
（`uv_mode='adaptive_fft'`，默认值），
或者均匀的 Faucher-Giguère+2011 背景（`uv_mode='global'`）。两者都把
背景作为下限，因此 `'global'` 就是没有恒星形成的星系在
`'adaptive_fft'` 下的极限。

投影模型的 UV 场由 SFR 图导出 — 这是方法定义的一部分 —
因此 `uv_mode` 对它们不适用，把
`'global'` 与 `*_map` 模型一起传入会抛出异常。

> **v0.3.0 中的变更。** 默认值原本是 `'global'`，而 `'adaptive_fft'` 与它
> 完全等价：构建插值器时把某一个轴的边界数组
> 复用到了全部三个轴上，于是几乎每个粒子都落在其范围之外，
> `fill_value=0.0` 返回零，只剩下背景。现在的默认值
> 是 `'adaptive_fft'`，名副其实。

## 验证

物理表达式是从 `TNGHIcalc` 逐字移植过来的，后者已
与 Diemer+2018 的参考实现逐位（bit-for-bit）验证一致
（`TNGHI_GD14_model/gd14_model.py`）— 证据链和复现命令
见 [验证记录](../verification.md)。
