# astrovodka.datacube — DataCube

定义输出数据立方体的网格：像素、通道与 WCS。

## 模块常量

- **`HIfreq`** = `1.420405751e9 * U.Hz` — HI 21-cm 静止频率。

## DataCube

### `__init__(self, *, n_px_x, n_px_y, n_channels, px_size, channel_width, spectral_centre=0 km/s, ra=0 deg, dec=0 deg, stokes_axis=False, coordinate_frame=ICRS(), specsys='icrs', spectral_unit=None)`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `n_px_x` | `int` | 沿 RA 方向的像素数。建议取偶数 |
| `n_px_y` | `int` | 沿 Dec 方向的像素数 |
| `n_channels` | `int` | 光谱通道数 |
| `px_size` | `Quantity[arcsec]` | 像素的角尺寸 |
| `channel_width` | `Quantity[km/s \| Hz]` | 通道宽度 |
| `spectral_centre` | `Quantity` | 光谱轴的中心（默认：0 km/s） |
| `ra` | `Quantity[deg]` | 数据立方体质心的 RA（默认：0） |
| `dec` | `Quantity[deg]` | 数据立方体质心的 Dec（默认：0） |
| `stokes_axis` | `bool` | 是否包含 Stokes 轴 |
| `coordinate_frame` | `BaseRADecFrame` | 坐标系（默认：ICRS） |
| `specsys` | `str` | 光谱参考系：`'icrs'`、`'gcrs'`、`'lsrk'` |
| `spectral_unit` | `str \| None` | `None`（默认）沿用 `channel_width` 的单位；也可强制为 `'MHz'` / `'km/s'` |

> **说明：**
> - 所有参数均为仅关键字参数。
> - `spectral_unit=None`（默认）与 MARTINI 行为一致：`channel_width` 为速度单位时得到
>   速度轴（WCS CTYPE3=`VRAD`），为频率单位时得到 `FREQ`。
> - `spectral_unit='MHz'` 会把速度单位的 `channel_width` 换算为频率。该换算是*微分*意义上的
>   （`Δf = f_rest·Δv/c`，即在 HI 线附近 20 km/s → 94.76 kHz），而不是该速度处的
>   绝对频率。
> - `channel_width` 必须带频率或速度单位，否则抛出 `ValueError`。
> - 解析出的单位保存在 `_spectral_unit` 中，并被 `copy()`、
>   `save_state()`/`load_state()`、`from_wcs()` 和 `Martini.reset()` 复用，
>   因此数据立方体的光谱轴不会悄悄改变类型。
> - 内部 `_array`：形状 `(n_px_x, n_px_y, n_channels)`，单位 `Jy/pix²`。

---

### `from_wcs(cls, input_wcs, specsys=None)`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `input_wcs` | `astropy.wcs.WCS` | 已有的 WCS（例如来自 FITS 头） |
| `specsys` | `str \| None` | 覆盖光谱参考系 |

**返回值：** `DataCube`

> 从 WCS 中提取全部维度。要求像素为正方形。缺少 RESTFRQ 时给出警告。

---

## 属性

| 属性 | 返回 | 说明 |
|----------|---------|-------------|
| `wcs` | `WCS` | 惰性构建的 WCS |
| `channel_mids` | `SpectralCoord` | 通道中心（原生单位） |
| `channel_edges` | `SpectralCoord` | 通道边界（n_channels+1 个值） |
| `velocity_channel_mids` | `Quantity[km/s]` | 以 km/s 表示的中心 |
| `velocity_channel_edges` | `Quantity[km/s]` | 以 km/s 表示的边界 |
| `frequency_channel_mids` | `Quantity[Hz]` | 以 Hz 表示的中心 |
| `frequency_channel_edges` | `Quantity[Hz]` | 以 Hz 表示的边界 |
| `spatial_slices` | `Iterator[Quantity]` | 迭代 (ny,nx) 切片 |
| `spectra` | `Iterator[Quantity]` | 按像素迭代 (nchan,) 谱 |

## 方法

### `save_state(self, filename, overwrite=False)`

把数组与 WCS 保存到 HDF5。需要 `h5py`。

### `load_state(cls, filename)`

从 `save_state()` 恢复。返回 `DataCube`。

### `copy(self)` → `DataCube`

深拷贝。

### `add_pad(self, pad: tuple[int,int])` / `drop_pad(self)`

波束卷积所需的填充区管理。由内部调用。
