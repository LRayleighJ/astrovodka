# astrovodka.noise — 噪声模型

生成噪声数据立方体。内置：`GaussianNoise`。

## GaussianNoise

### `__init__(self, rms=1 Jy/beam, seed=0)`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `rms` | `Quantity[Jy/beam]` | 卷积后的 RMS（默认：1） |
| `seed` | `int` | 随机种子（默认：0） |

### `generate(self, datacube, beam)` → `Quantity[Jy/arcsec²]`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `datacube` | `DataCube` | 数组形状与像素尺度 |
| `beam` | `_BaseBeam` | 用于 RMS 缩放的 FWHM |

> 噪声在卷积**之前**加入。RMS 按 `rms_pre = rms_post × 2.19568 × √(π × σ_maj_px × σ_min_px)` 缩放。

### `reset_rng(self)`

用原始种子重新播种，以便重复得到同一实现。

---

> **自定义噪声：** 从 `_BaseNoise` 继承，实现 `generate(datacube, beam)`。
