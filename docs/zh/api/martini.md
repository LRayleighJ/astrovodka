# astrovodka.martini — 核心流水线

`Martini` 类负责编排模拟观测流水线，把 source、DataCube、beam、noise、spectral model 与 SPH kernel 串联起来。

## Martini

### `__init__(self, *, source, datacube, beam=None, noise=None, sph_kernel, spectral_model, quiet=False)`

**参数：**

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `source` | `SPHSource` | 发射 HI 的源对象（必需） |
| `datacube` | `DataCube` | 数据立方体的网格定义（必需） |
| `beam` | `_BaseBeam \| None` | 望远镜波束。`None` 表示跳过卷积（默认：None） |
| `noise` | `_BaseNoise \| None` | 噪声模型。`None` 表示跳过加噪（默认：None） |
| `sph_kernel` | `_BaseSPHKernel` | SPH 平滑核（必需） |
| `spectral_model` | `_BaseSpectrum` | 21-cm 谱线模型（必需） |
| `quiet` | `bool` | 抑制标准输出（默认：False） |

**返回值：** None

> **说明：**
> - 所有参数均为仅关键字参数（keyword-only）。
> - 初始化过程中会裁剪粒子（按空间、按光谱、按质量），以跳过不贡献信号的粒子。
> - `DataCube._array` 的单位会随流程变化：插入源之后为 `Jy/pix²` → 卷积之后为 `Jy/beam`。

---

### `insert_source_in_cube(self, skip_validation=False, progressbar=None, ncpu=1)`

**参数：**

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `skip_validation` | `bool` | 跳过 SPH 核精度校验（默认：False） |
| `progressbar` | `bool \| None` | 是否显示 tqdm 进度条。`None` 表示自动判断 |
| `ncpu` | `int` | 并行进程数。需要 `multiprocess`（默认：1） |

**返回值：** None

> **说明：**
> - 这是开销最大的一步 — 每个像素都要遍历附近的粒子。
> - 并行使用 `multiprocess`（dill），而不是 `multiprocessing`（pickle）。
> - 并行粒度是逐像素的：每个 worker 独立处理一个像素。
> - 会报告质量守恒情况：由流量累加得到的质量与输入 HI 质量之差应在约 1% 以内。

---

### `convolve_beam(self)`

**返回值：** None

> **说明：**
> - 要求已设置 `beam`。若为 None 则给出警告并跳过。
> - 使用 `scipy.signal.fftconvolve` 做基于 FFT 的卷积。
> - 卷积后会移除填充区；数组单位变为 `Jy/beam`。

---

### `add_noise(self)`

**返回值：** None

> **说明：**
> - 要求同时设置 `beam` 和 `noise`。否则给出警告并跳过。
> - 噪声在卷积**之前**加入。RMS 会自动缩放：`rms_pre = rms_post × 2.19568 × √(π × σ_maj × σ_min)`。

---

### `write_fits(self, filename, overwrite=True, obj_name='MOCK', channels=None, dtype=None)`

**参数：**

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `filename` | `str` | 输出路径。若未带 `.fits` 后缀则自动追加 |
| `overwrite` | `bool` | 允许覆盖（默认：True） |
| `obj_name` | `str` | OBJECT 头（最多 16 个字符） |
| `dtype` | `str \| np.dtype \| None` | 输出类型。`None` 表示保持原生类型 |

**返回值：** None

> **说明：** 仪器头为 `ASTROVODKA`。轴顺序转置：`(RA,Dec,Chan)` → FITS `(Chan,Dec,RA)`。

---

### `write_hdf5(self, filename, overwrite=True, memmap=False, compact=False, channels=None)`

**参数：**

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `memmap` | `bool` | 在内存中创建文件并返回其句柄 |
| `compact` | `bool` | 省略坐标数组以节省空间 |

**返回值：** `h5py.File | None`

---

### `write_beam_fits(self, filename, overwrite=True, channels=None)`

**返回值：** None。若 beam 为 None 则抛出 `ValueError`。

---

### `preview(self, max_points=5000, fig=1, lim=None, vlim=None, point_scaling='auto', title='', save=None)`

**参数：**

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `max_points` | `int` | 每个面板最多绘制的粒子数（随机抽样） |
| `lim` | `str \| Quantity[deg] \| None` | 空间范围限制。`'datacube'` 表示取数据立方体的边界 |
| `vlim` | `str \| Quantity[km/s] \| None` | 速度范围限制 |
| `point_scaling` | `str` | `'auto'`（按质量 + 平滑长度）或 `'fixed'` |
| `save` | `str \| None` | 把图保存到文件 |

**返回值：** `matplotlib.figure.Figure`

> 三联图：天球投影 + 2 幅 PV 图。红框表示 DataCube 的范围。

---

### `reset(self)`

**返回值：** None。把 DataCube 以零重新初始化，并保留源与各子模块的配置。

---

## GlobalProfile

仅用于空间积分谱的简化类。不需要 beam、noise 或 SPH kernel。

### `__init__(self, *, source, spectral_model, n_channels, channel_width, spectral_centre=0 km/s, quiet=False)`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `source` | `SPHSource` | HI 源 |
| `spectral_model` | `_BaseSpectrum` | 谱线模型 |
| `n_channels` | `int` | 通道数 |
| `channel_width` | `Quantity[km/s \| Hz]` | 通道宽度 |
| `spectral_centre` | `Quantity` | 光谱中心（默认：0 km/s） |
| `quiet` | `bool` | 抑制输出 |

> 使用 1×1 像素的 DataCube 配合 DiracDeltaKernel。所有粒子都会贡献信号。

### `spectrum`（属性）

**返回值：** `Quantity[Jy]` — 一维流量密度数组。惰性求值。

### `plot_spectrum(self, fig=1, title='', channels='velocity', show_vsys=True, save=None)`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `channels` | `str` | `'velocity'` 或 `'frequency'` |
| `show_vsys` | `bool` | 是否在系统速度处画线 |
| `save` | `str \| None` | 把图保存到文件 |

**返回值：** `Figure`
