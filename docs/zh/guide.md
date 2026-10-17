# 用户指南

astrovodka 如何工作，以及如何使用它。安装见
[installation.md](installation.md)；各类的签名见
[API 参考](api/index.md)；HI/H₂ 模型的深入讨论见
[hi-models.md](hi-models.md)。

- [流水线](#流水线)
- [源](#源)
- [两种距离](#两种距离)
- [数据立方体](#数据立方体)
- [子模块](#子模块)
- [完整示例](#完整示例)
- [常见任务](#常见任务)
- [陷阱](#陷阱)
- [与 MARTINI 的差异](#与-martini-的差异)

---

## 流水线

astrovodka 是面向对象的：`Martini` 持有一个**源**、一个**数据立方体**和
三个子模块，并按顺序运行它们。

```
    source ────────┐
    datacube ──────┤
    beam ──────────┤     ┌──────────────────────┐
    noise ─────────┼────▶│       Martini        │
    spectral_model ┤     └──────────────────────┘
    sph_kernel ────┘        │      │       │
                            ▼      ▼       ▼
                   insert_source_in_cube()  →  粒子落到像素上
                            │
                      add_noise()           →  （可选）每个像素的噪声
                            │
                    convolve_beam()         →  （可选）望远镜波束
                            │
                      write_fits()          →  一个能在 DS9 中打开的数据立方体
```

每一步都是独立的方法，因此你可以在其中任意一步之后停下来。停在
`insert_source_in_cube()` 之后得到的是一个无噪声、未卷积的数据立方体 ——
检查物理时你想要的正是这个，因为这时没有噪声与信号竞争，
也没有波束把信号重新分配。

## 源

源抽象出“发出 21 cm 谱线辐射的那个东西”。已实现的有两个。

### `SPHSource` —— 任意 SPH 数据

```python
SPHSource(
    redshift=0.01, ra=30*U.deg, dec=45*U.deg,
    T_g=..., mHI_g=..., xyz_g=..., vxyz_g=..., hsm_g=...,
)
```

粒子数组就是普通的 `astropy` 量：

| 数组 | 单位 | 含义 |
|---|---|---|
| `xyz_g` | kpc | 位置。**`x` 是视线方向**；`y`-`z` 是天空平面 |
| `vxyz_g` | km/s | 速度，采用相同的轴约定 |
| `mHI_g` | M☉ | 每个粒子的 HI 质量 |
| `T_g` | K | 温度，用于热谱线宽度 |
| `hsm_g` | kpc | SPH 平滑长度，**以 FWHM 计**，而不是核半径 |

最后一点值得强调。各模拟对“平滑长度”到底是什么并无共识 —— 有的存的是
核降到零处的半径，有的存 FWHM，还有的存别的什么。MARTINI 选择 FWHM
作为自己的约定。如果你传给它的其实是核半径，每个粒子被平滑的区域
都会与预期的大小不同，而且不会有任何报错。

### `TNGSource` —— 从本地磁盘读取 IllustrisTNG

```python
TNGSource(
    basePath="/path/to/TNG/TNG-100/",
    snapNum=99, subID=113,
    hi_model="K13_vol",
    observed_redshift=0.02,
)
```

`TNGSource` 直接从文件系统读取快照，运行选定的 HI/H₂ 模型，并把结果
传给 `SPHSource`。它额外提供：

- **`hi_model`** —— 四种模型之一（见 [hi-models.md](hi-models.md)）。
  在任何磁盘 I/O *之前*就已校验，因此拼写错误会立即失败，
  而不会在读入几个 GB 之后才发现。
- **`observed_redshift`** —— 见[两种距离](#两种距离)。
- **`load_subhalo_only`**（默认 `True`）—— 只加载该子晕的粒子，
  而不是它所在的整个 FOF 群。快得多，通常也正是你想要的。
- **以暗物质质心为中心**，而不是星表中的 `SubhaloPos`。

构建完成之后，这些属性值得了解：

| 属性 | 含义 |
|---|---|
| `m_HI_total` | HI 总质量 [M☉] |
| `m_H2_g` | 每个粒子的 H₂ 质量 |
| `f_mol_g`, `f_neutral_g` | 逐单元的比例 |
| `u_mw_g` | 每个粒子的 UV 场，单位为 Draine 单位 |
| `hi_diagnostics` | 模型诊断信息，包括 `f_gas_in_map` |
| `D_L`, `D_A`, `vsys` | 宇宙学距离与系统速度 |
| `source_redshift` | 快照自身的红移 |

## 两种距离

这是 astrovodka 与 MARTINI 区别最大的地方，弄错的代价很高，
因此值得讲清楚。

**`D_L` 与 `D_A` 不能互换，二者相差 `(1+z)²`。**

| 距离 | 用于 | 弄错的代价 |
|---|---|---|
| `D_L`，光度距离 | 把 HI 质量换算成流量；对数据立方体做归一化 | `(1+z)⁴` |
| `D_A`，角直径距离 | 所有角尺度：像素尺度、源的延展、波束 | `(1+z)²` |

astrovodka 用 `astropy.cosmology.Planck15`（或你传入的宇宙学）计算
这两个距离，并把各自用在该用的地方。

### 两个红移

一个 TNG 快照只有一个红移，而它并不总是你想要观测的那个红移。
`z ≈ 0` 的快照角直径距离趋近于零，因此其中的星系张开的角几乎可以忽略，
根本无法分辨。

`observed_redshift` 把这两个角色分开：

| 量 | 来自 | 原因 |
|---|---|---|
| 单位换算、UV 背景、HI/H₂ 模型 | **快照的**红移 | 那是气体的物理状态 |
| `D_L`、`D_A`、`vsys`、角尺度、流量 | `observed_redshift` | 那是观测的几何 |

```python
source = TNGSource(basePath=..., snapNum=99, subID=113, observed_redshift=0.05)
print(source.source_redshift)   # ~0，即快照的红移
print(source.D_L, source.D_A)   # 按 z = 0.05 计算
```

### 系统速度约定

`vsys` 采用**射电约定**，即 `c·z/(1+z)`，这正是 WCS 光谱轴所期望的。
这比听上去更要紧：在 z = 0.5 时，另一种做法（`H(z)·D_C`）大了 74%，
在 z ≈ 0.79 以上会超过 `c`，并在高红移处给出负频率。

## 数据立方体

```python
DataCube(
    n_px_x=256, n_px_y=256, n_channels=64,
    px_size=5*U.arcsec,
    channel_width=20*U.km/U.s,
    spectral_centre=source.vsys,
)
```

关于光谱轴有两点值得了解。

**单位默认跟随 `channel_width`。** 给一个速度宽度就得到速度轴；
WCS 会记录其中的换算。传入 `spectral_unit="MHz"` 可以改为强制使用频率轴。

**宽度按微分方式换算。** 速度宽度换到频率轴时按 `Δf = f_rest·Δv/c`
换算，而不是把两个端点变换再相减 —— 两者并不相同，用错的那个会
悄无声息地丢掉几乎全部信号（见[陷阱](#陷阱)）。

光谱轴始终要覆盖谱线，而不只是气体。`Martini` 会裁剪掉谱线轮廓放不进
通道范围的粒子，而轮廓宽 `±4σ`，因此轴必须比粒子速度范围宽大约 `8σ`。
对于有恒星形成的子晕，最热的单元达到 σ ≈ 100 km/s，这意味着大约
800 km/s 的余量。弄错这一点会悄无声息地丢掉质量 ——
质量回收那一行会掉到 100% 以下，而不会有别的任何提示。

无论当初是用哪种单位构建的，都可以按任一单位把这条轴读回来：

```python
datacube.velocity_channel_edges   # km/s
datacube.frequency_channel_edges  # Hz
```

## 子模块

| 组件 | 选项 | 说明 |
|---|---|---|
| **波束** | `GaussianBeam(bmaj, bmin, bpa, truncate)` | `bmin` 默认为 `bmaj` |
| **噪声** | `GaussianNoise(rms, seed)` | 在卷积*之前*加入，因此波束会使其相关 |
| **谱线模型** | `GaussianSpectrum(sigma=...)`、`DiracDeltaSpectrum()` | `sigma="thermal"` 表示逐粒子使用 `sqrt(k_B T / m_p)` |
| **SPH 平滑核** | `CubicSplineKernel`、`WendlandC2Kernel`、`WendlandC6Kernel`、`QuarticSplineKernel`、`GaussianKernel`、`DiracDeltaKernel` | 选与你的模拟所用核相匹配的那一个 |

SPH 平滑核是**自适应的**：当像素与平滑长度之比使主要近似不再准确时，
每个核都会依次回退到 `DiracDeltaKernel`，再回退到
`GaussianKernel(truncate=6)`，并且会说明这一点。当该比值处于临界时，
校验会给出警告；在 `insert_source_in_cube()` 上设置
`skip_validation=True` 可以把它关掉。

## 完整示例

```python
import astropy.units as U
from astrovodka import Martini, DataCube, SubhaloDataProducts
from astrovodka.sources import TNGSource
from astrovodka.beams import GaussianBeam
from astrovodka.spectral_models import GaussianSpectrum
from astrovodka.sph_kernels import CubicSplineKernel

basePath, snapNum, subID = "/path/to/TNG/TNG-100/", 99, 113

# --- 三维数据产品：UV 场，以及实空间和红移空间中的 DM/HI/H2 密度
dp = SubhaloDataProducts.from_subhalo(
    basePath=basePath, snapNum=snapNum, subID=subID,
    ngrid=64, hi_model="K13_vol",
)
dp.save_hdf5(f"sub{subID}_products.hdf5")

# --- 源
source = TNGSource(
    basePath=basePath, snapNum=snapNum, subID=subID,
    hi_model="K13_vol", observed_redshift=0.05, load_subhalo_only=True,
)

# --- 数据立方体。光谱轴沿用 channel_width 的单位（此处为 km/s）。
datacube = DataCube(
    n_px_x=128, n_px_y=128, n_channels=64,
    px_size=5 * U.arcsec,
    channel_width=20 * U.km / U.s,
    spectral_centre=source.vsys,
)

# --- 运行
M = Martini(source=source, datacube=datacube,
            beam=GaussianBeam(bmaj=15 * U.arcsec),
            spectral_model=GaussianSpectrum(sigma="thermal"),
            sph_kernel=CubicSplineKernel())
M.insert_source_in_cube(progressbar=False)
M.convolve_beam()
M.write_fits(f"sub{subID}_cube.fits")

# --- 按任一单位读回这条轴
print(f"channel 0: {datacube.frequency_channel_edges[0].to(U.MHz):.4f} MHz")
print(f"channel 0: {datacube.velocity_channel_edges[0].to(U.km/U.s):.1f} km/s")
```

`insert_source_in_cube()` 会打印放进立方体的质量占源总量的比例。
**要读那一行。** 它应该显示 100%。

## 常见任务

### 改变角分辨率

用不同的波束和像素尺寸重建数据立方体。干涉阵能分辨结构，
但会漏掉延展辐射；单口径天线则相反。

```python
datacube = DataCube(..., px_size=beam_fwhm / 3, ...)   # 每个波束约 3 个像素
```

### 把星系放到多个红移

```python
for z in (0.0, 0.02, 0.05, 0.1):
    source = TNGSource(basePath=..., snapNum=99, subID=113, observed_redshift=z)
    ...
```

每次都新建一个 `TNGSource`。构建数据立方体**会就地裁剪源**，
因此一个源对象只适用于一次观测 —— 见[陷阱](#陷阱)。

### 批量处理子晕

见 `examples/en/batch_subhalos.py`。简而言之：读取 FOF 群表，
按气体单元数过滤，然后循环。

```python
import h5py
with h5py.File(basePath + "groups_099/fof_subhalo_tab_099.0.hdf5", "r") as f:
    n_gas = f["SubhaloLenType"][:, 0]
ids = [i for i, n in enumerate(n_gas) if 64 <= n <= 5000]
```

**不要**取星表的前 N 个 —— 它按质量降序排列，因此最前面的条目
是星系团核心。要先过滤。

### 匹配已有的观测

```python
datacube = DataCube.from_wcs(fits.getheader("observation.fits")[0])
```

这样数据立方体就继承了真实的像素网格、通道宽度和光谱参考系。

### 用一条谱代替数据立方体

`GlobalProfile` 完全跳过数据立方体 —— 不需要波束、不需要噪声、
不需要 SPH 平滑核：

```python
from astrovodka import GlobalProfile
gp = GlobalProfile(source=source, spectral_model=GaussianSpectrum(sigma=7*U.km/U.s),
                   n_channels=64, channel_width=10*U.km/U.s,
                   spectral_centre=source.vsys)
gp.plot_spectrum()
```

## 陷阱

这四个都耗费过实实在在的调试时间。它们全都继承自 MARTINI。

### 1. `Martini` 会就地裁剪源

`Martini.__init__` 调用 `_prune_particles`，后者通过 `source.apply_mask`
把落在立方体之外的粒子从**源对象本身**移除 —— 而不只是从那一次观测中
移除。先建一个小的临时立方体检查点什么，再用同一个源建一个大的科学
立方体，第二个拿到的就是一个被悄悄截断的源。唯一的症状是质量回收
那一行掉到 100% 以下。

**一个源，一次观测。** 每个数据立方体都重建源。（裁剪并不会连同
`mHI_g` 一起把 `m_H2_g` 也掩掉，因此被裁剪过的源在内部也是不一致的。）

### 2. `write_fits` 和 `write_hdf5` 会丢弃填充区

两者都会调用 `datacube.drop_pad()`，因为填充区是为卷积而存在的。
在流水线中途写出一个中间产物，下一次 `convolve_beam()` 就会抛出
`ValueError: datacube padding insufficient for beam convolution`。把它加回来：

```python
M.datacube.add_pad(M.beam.needs_pad())
```

### 3. 谱线模型不能在形状不同的数据立方体之间共享

`GaussianSpectrum` 会缓存自己的 `spectra` 数组，而 `Martini` 只在该数组
为 `None` 时才重新初始化。同一个实例先用在一个 32 通道的数据立方体上，
再用在一个 42 通道的上面，就会抛出 `cannot reshape array of size ...`。

**每个 `Martini` 都新建一个谱线模型。**（如果你已经踩上了，
`M.init_spectra()` 是应急出口。）

### 4. 光谱轴必须覆盖谱线，而不只是气体

见[数据立方体](#数据立方体)一节：在粒子速度范围之外留出大约 `8σ`，
对有恒星形成的子晕大约就是 800 km/s。

## 与 MARTINI 的差异

观测流水线是 MARTINI 的：astrovodka 导入 `astromartini` 包，只在下表所列行为不同的
地方继承并覆盖，其余原样转发。不同之处：

| | MARTINI | astrovodka |
|---|---|---|
| **距离** | 线性哈勃定律 `v = H₀D` | `Planck15`：流量用 `D_L`，角尺度用 `D_A` |
| **系统速度** | `h·100·D` | 射电约定 `c·z/(1+z)` |
| **HI 质量** | 基于压强（Leroy+2008） | 四种可选的 Diemer+2018 模型 |
| **H₂ 质量** | 未暴露 | `source.m_H2_g`，逐粒子 |
| **UV 场** | — | 逐星系自适应 FFT 场，下限取 UV 背景 |
| **红移覆盖** | — | `observed_redshift`，把物理与几何解耦 |
| **TNG 访问** | Web API + JupyterLab | 本地磁盘 |
| **子晕加载** | 完整 FOF 群 | 仅子晕（默认）或完整群 |
| **居中方式** | 星表 `SubhaloPos` | 暗物质质心 |
| **三维产品** | — | `SubhaloDataProducts`：UV、DM、HI、H₂；实空间与红移空间 |
| **源** | EAGLE、Simba、FIRE、Magneticum、Colibre、SWIFT、… | 仅 TNG |

`SPHSource` 同时接受 `distance` + `h`（MARTINI 的接口）和
`redshift` + `cosmology`。在距离模式下 `D_L == D_A == distance`，
因此这条路径的行为与上游完全一致。

### 未实现的部分

- **只转发 TNG 源。** 上游的 EAGLE、Simba、FIRE、Magneticum、Colibre、SWIFT
  和 `CombinedSource` 类仍可从 `martini.sources` 使用，但 `astrovodka.sources`
  刻意不封装项目无法测试的源。
- **一个数据立方体中的多星系观测**，也就是 `CombinedSource` 的用途。
- **大子晕的内存分块** —— 三维网格是整体构建的。
