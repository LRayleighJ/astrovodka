# astrovodka

**面向 IllustrisTNG 的 HI（21 cm）模拟观测数据立方体。**

> **astrovodka 基于 [MARTINI](https://github.com/kyleaoman/martini)
> （PyPI 上的 `astromartini`）改动而来，不重复实现它的观测流程。**
> 数据立方体、波束、谱模型、SPH 核、噪声、`Martini` 驱动和 `SPHSource` 全部来自
> `astromartini` 包 —— 它是本项目的必需依赖，astrovodka 只导入并继承它，不复制
> 其代码。使用 astrovodka 时请同时引用 MARTINI。

astrovodka 把 SPH 粒子数据变成合成的射电数据立方体，在 MARTINI 之上加入两处
关键改动：距离改用 `astropy.cosmology` 的宇宙学计算，而不是线性哈勃定律；
中性气体质量改用有物理依据的 HI/H₂ 模型，而不是基于气压的分区。

[English](README.md) · **简体中文**

[![CI](https://github.com/LRayleighJ/astrovodka/actions/workflows/ci.yml/badge.svg)](https://github.com/LRayleighJ/astrovodka/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: GPLv3](https://img.shields.io/badge/license-GPLv3-green.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.11198185-blue.svg)](https://doi.org/10.5281/zenodo.11198185)

<img src="docs/images/sub113_projections.png" alt="TNG-100 子晕 113 的暗物质、HI 与 H2 密度投影，三个正交视角" width="100%">

*TNG-100 快照 99 子晕 113 的三个正交视角，由 `SubhaloDataProducts` 在默认的
`K13_vol` 模型下生成 —— 自上而下依次为暗物质、HI 和 H₂（后两行同一空间尺度）。
分子气体集中在内侧的密致盘上，原子气体则延伸到远得多的地方。*

---

## 它能做什么

| | |
|---|---|
| **宇宙学距离** | 光度距离 `D_L` 负责流量归一化，角直径距离 `D_A` 负责所有角尺度，两者都取自 `Planck15`。`z ≈ 0` 的快照里的星系，可以放到任意观测红移上。 |
| **四种可选 HI/H₂ 模型** | `GD14_vol`、`GD14_map`、`K13_vol`（默认）、`K13_map` —— Diemer+2018 对 Gnedin & Draine 2014 与 Krumholz+ 2013 的两两组合，各有体密度与投影两个版本，一个关键字切换。 |
| **自适应 FFT 紫外场** | 由星系自身的恒星形成率分布做三维 FFT 卷积，逐星系算出星际辐射场，并以 Faucher-Giguère+2011 背景为下限。 |
| **本地磁盘读 TNG** | `TNGSource` 直接从文件系统读快照，不走 Web API，不需要网络，也不需要 JupyterLab。 |
| **三维数据产品** | `SubhaloDataProducts` 输出 UV、DM、HI、H₂ 的网格化密度，实空间与红移空间都有。 |
| **经过验证** | 四种模型在真实 TNG-100 子晕上与参考实现**逐位一致**，见[验证](#验证)。 |

---

## 安装

```bash
git clone https://github.com/LRayleighJ/astrovodka.git
cd astrovodka

# 推荐用 conda，可复现出验证过的环境
conda env create -f environment.yml
conda activate astrovodka
```

或者用 pip，Python ≥ 3.9 即可：

```bash
python -m venv .venv && source .venv/bin/activate
python -m pip install -e ".[dev]"
```

| 依赖 | 版本 |
|---|---|
| Python | ≥ 3.9 |
| **astromartini** | 必需 —— 观测引擎，astrovodka 调用它 |
| numpy、scipy、astropy、tqdm | 自动安装 |
| **h5py** | 可选依赖，但**任何 TNG 工作流都必须有** —— 读快照靠它 |
| matplotlib、multiprocess | 可选：出图、多核并行 |

TNG 数据从本地磁盘读。`basePath` 要指向**包含** `groups_<snap>/` 与
`snapdir_<snap>/` 的那一级目录；TNG-100 单个快照的粒子数据约 1–2 TB，群目录约 13 GB。

**完整说明 —— 三种安装方式、可选依赖、离线/隔离网络集群、TNG 数据目录结构、
疑难排解表 —— 见 [`docs/zh/installation.md`](docs/zh/installation.md)。**

---

## 快速上手

### 解析星系，不需要模拟数据

```python
import numpy as np
import astropy.units as U
from astrovodka import Martini, DataCube
from astrovodka.sources import SPHSource
from astrovodka.beams import GaussianBeam
from astrovodka.spectral_models import GaussianSpectrum
from astrovodka.sph_kernels import CubicSplineKernel

n = 1000
rng = np.random.default_rng(42)
phi, r = rng.random(n) * 2 * np.pi, rng.random(n) * 10 * U.kpc

source = SPHSource(
    redshift=0.01,                      # 宇宙学模式：由 Planck15 给出 D_L 与 D_A
    ra=30 * U.deg, dec=45 * U.deg,
    T_g=np.full(n, 1e4) * U.K,
    mHI_g=np.full(n, 1e8) * U.Msun,
    xyz_g=U.Quantity([np.zeros(n) * U.kpc, r * np.cos(phi), r * np.sin(phi)]),
    vxyz_g=U.Quantity([np.zeros(n) * U.km / U.s,
                       -100 * np.sin(phi) * U.km / U.s,
                       100 * np.cos(phi) * U.km / U.s]),
    hsm_g=np.full(n, 1.0) * U.kpc,      # 平滑长度，按 FWHM 给
)

# 谱轴单位默认跟随 channel_width，除非显式指定 spectral_unit
datacube = DataCube(
    n_px_x=128, n_px_y=128, n_channels=64,
    px_size=2 * U.arcsec,
    channel_width=5 * U.km / U.s,
    spectral_centre=source.vsys,
    ra=source.ra, dec=source.dec,
)

M = Martini(source=source, datacube=datacube,
            beam=GaussianBeam(bmaj=8 * U.arcsec),
            spectral_model=GaussianSpectrum(sigma="thermal"),
            sph_kernel=CubicSplineKernel())
M.insert_source_in_cube()
M.convolve_beam()
M.write_fits("cube.fits")
```

### 真实 TNG 子晕

```python
from astrovodka.sources import TNGSource

source = TNGSource(
    basePath="/path/to/TNG/TNG-100/",
    snapNum=99,
    subID=113,
    hi_model="K13_vol",       # 或 K13_map / GD14_vol / GD14_map
    observed_redshift=0.02,   # 观测距离，与快照自身红移解耦
    load_subhalo_only=True,   # 只读子晕粒子，比读整个 FOF 群快得多
)

print(source.m_HI_total, source.m_H2_g.sum(), source.f_mol_g.mean())
print(source.D_L, source.D_A, source.vsys)
```

> **子晕要挑着选。** TNG 目录按质量降序排列，`subID=0` 就是团核心 ——
> 约 3800 万个气体单元，读一次超过 10 GB。100 多号的子晕是气体丰富的卫星星系，
> 几秒就能读完。

两个例子都有可直接运行的版本，另有一个批量脚本：

```bash
python examples/zh/astrovodka_demo.py                            # 解析星系
python examples/zh/astrovodka_demo.py --subID 113 --snapNum 99   # 需要 $TNG_BASE_PATH
python examples/zh/batch_subhalos.py --limit 20                  # 遍历目录
```

---

## 四种 HI/H₂ 模型

凡是构造源或计算质量的地方，都用关键字选择。

| `hi_model` | 处方 | 投影方式 | 需要恒星粒子 | 开销 |
|---|---|---|---|---|
| `K13_vol`（默认） | Krumholz+ 2013 | 体密度（三维） | 否 | 中 |
| `K13_map` | Krumholz+ 2013 | 投影 | 是 | 较高 |
| `GD14_vol` | Gnedin & Draine 2014 | 体密度（三维） | 否 | 中 |
| `GD14_map` | Gnedin & Draine 2014 | 投影 | 是 | 较高 |

**体密度**模型由局部的气体体密度、金属丰度和紫外场算出分子比例。**投影**模型改用
面密度，更接近观测的实际做法，但需要先把星系投影到网格上 —— 于是就必须对落在网格
之外的气体做出决定。这部分气体的质量占比会以 `f_gas_in_map` 报告，`TNGSource`
在它低于 0.9 时给出警告。

同一套物理也可以脱离观测流程单独使用：

```python
from astrovodka.himodel import METHODS, HIModelConfig, compute_hi_mass

METHODS    # ('GD14_vol', 'GD14_map', 'K13_vol', 'K13_map')

cfg = HIModelConfig(h=0.6774, z=0.0, boxsize_ckpc_per_h=75000.0, ngrid=64)
result = compute_hi_mass(gas, model="K13_map", cfg=cfg, stars=stars, sub=sub)
result.m_hi_msun, result.m_h2_msun
```

各模型假设了什么、紫外场怎么构造、数值参数怎么调，见
[`docs/zh/hi-models.md`](docs/zh/hi-models.md)。

---

## 与 MARTINI 的差别

astrovodka 并不是 fork，而是构建在 MARTINI 之上：观测流程（数据立方体、波束、谱模型、
SPH 核、噪声、`Martini` 驱动和 `SPHSource`）就是 `astromartini` 包本身。astrovodka
导入它，只在下表所列行为不同的地方继承并覆盖，其余原样转发；包内没有任何 MARTINI
代码的副本。不同的地方有：

| | MARTINI | astrovodka |
|---|---|---|
| **距离** | 线性哈勃定律 `v = H₀D` | `Planck15`：`D_L` 管流量，`D_A` 管角尺度 |
| **系统速度** | `h·100·D` | 射电约定 `c·z/(1+z)`，与谱轴一致 |
| **HI 质量** | 基于气压（Leroy+2008） | 四种可选 Diemer+2018 模型 |
| **H₂ 质量** | 未提供 | `source.m_H2_g`，逐粒子 |
| **紫外场** | — | 逐星系自适应 FFT 场，以背景为下限 |
| **红移覆盖** | — | `observed_redshift`：把 `z ≈ 0` 的星系放到任意距离 |
| **TNG 读取** | Web API + JupyterLab | 本地磁盘 |
| **子晕加载** | 整个 FOF 群 | 仅子晕（默认）或整个群 |
| **中心定位** | 目录里的 `SubhaloPos` | 暗物质质心 |
| **三维产品** | — | `SubhaloDataProducts`：UV、DM、HI、H₂，实空间与红移空间 |
| **已实现的源** | EAGLE、Simba、FIRE、Magneticum、Colibre、SWIFT…… | 仅 TNG |

依赖是**单向且显式**的：astrovodka 调用 MARTINI，反过来不成立；
`astrovodka.himodel` 是纯物理，既不导入观测层，也不导入 MARTINI，因此质量模型可以
单独使用、单独测试，完全不需要数据立方体。`tests/test_upstream_dependency.py` 强制
检查：`astromartini` 被声明为依赖且确实被导入；观测类要么就是上游对象、要么是它们的
子类，而不是重新实现；`himodel` 保持独立。

---

## 验证

四种模型是从一份参考实现逐表达式移植的（那份参考实现本身已对照过作者的原始代码），
并在真实数据上做了验证：

- **逐位一致**：在随机气体单元集合上、以及在真实 TNG-100 子晕上（4448、120、86902
  个气体单元），与参考实现逐位相同。
- **质量守恒**：流程会打印源质量经过粒子裁剪后存活的比例，教程断言它必须是 100%。
- **端到端 FITS 反算**：把写出的立方体积分回来的质量，与模型自己的总质量一致。
- **125 个单元测试**，不需要 TNG 数据。

完整记录 —— 包括上一版紫外插值中一个被实测出来的 bug（它让自适应场静默退化成了
均匀背景）—— 见 [`docs/zh/verification.md`](docs/zh/verification.md)。

---

## 文档

| | |
|---|---|
| [安装与环境配置](docs/zh/installation.md) | 三种安装方式、可选依赖、离线集群、TNG 数据目录、疑难排解 |
| [用户指南](docs/zh/guide.md) | 核心概念，以及一个完整的例子 |
| [HI/H₂ 模型](docs/zh/hi-models.md) | 四种模型各自假设什么、怎么调参 |
| [API 参考](docs/zh/api/index.md) | 全部公开的类与函数 |
| [验证记录](docs/zh/verification.md) | 模型是怎么验证的 |
| [更新日志](CHANGELOG.md) | 版本历史 |
| English docs | [`docs/en/`](docs/en/index.md) |

用 `python docs/build_html.py` 可以把 Markdown 源文件生成单页 HTML，输出在
`docs/html/`。

---

## 示例

| | |
|---|---|
| [`examples/zh/tutorial_TNG100-HI-mock-observations.ipynb`](examples/zh/tutorial_TNG100-HI-mock-observations.ipynb) | **从这里开始。** 在真实 TNG-100 子晕上的完整教程：怎么选子晕、四种模型怎么切换、怎么建数据立方体以及其中的坑、怎么验证结果、三维产物、可视化、疑难排解。 |
| [`examples/zh/astrovodka_demo.py`](examples/zh/astrovodka_demo.py) | 可直接运行的脚本：默认跑解析星系，加 `--basePath` 就跑真实 TNG 子晕。 |
| [`examples/zh/batch_subhalos.py`](examples/zh/batch_subhalos.py) | 遍历子晕目录，每个目标导出一个立方体。 |
| [`examples/en/`](examples/en/) | 同样三个示例的英文版。 |

Notebook 不保存输出，请自己跑。先指向你的数据：

```bash
export TNG_BASE_PATH=/path/to/TNG/TNG-100/
```

---

## 致谢与引用

### 引用 astrovodka

如果 astrovodka 对你的发表有帮助，请把 [`NOTICE`](NOTICE) 中列出的模型论文与下面的
MARTINI 论文一并引用 —— HI/H₂ 物理来自 Gnedin & Draine (2014)、Krumholz+ (2013) 和
Diemer+ (2018)，不是本包发明的。

### MARTINI

astrovodka 建立于 [MARTINI](https://github.com/kyleaoman/martini) 之上，并包含其
修改后的代码。请引用：

- Oman, K. A. 2024, *J. Open Source Softw.*, 9, 6860.
  [doi:10.21105/joss.06860](https://doi.org/10.21105/joss.06860)
- Oman, K. A., et al. 2019, *MNRAS*, 482, 821.
  [doi:10.1093/mnras/sty2687](https://doi.org/10.1093/mnras/sty2687)
- Oman, K. A. 2019, *ASCL*, ascl:1911.005.
  [ascl.net/1911.005](https://ascl.net/1911.005)

机器可读的 [`CITATION.cff`](CITATION.cff) 也已提供。

---

## 许可证

GNU 通用公共许可证第 3 版或更新版本 —— 见 [`LICENSE`](LICENSE)。

astrovodka 是 MARTINI 的衍生作品，MARTINI 采用同一许可证；衍生代码继承该许可证。
[`NOTICE`](NOTICE) 按 GPLv3 第 5(a) 条的要求，记录了各部分代码的来源。

---

## 联系方式

**Zerui Liu** — liuzr@bao.ac.cn

问题与 PR 欢迎提交到 <https://github.com/LRayleighJ/astrovodka/issues>。
开发环境的搭建见 [`CONTRIBUTING.md`](CONTRIBUTING.md)。
