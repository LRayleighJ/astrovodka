# astrovodka.beams — 望远镜波束模型

对望远镜波束进行建模。内置：`GaussianBeam`（椭圆二维高斯）。

实现来自 MARTINI：从 `martini.beams` 导入并在此转发。

## GaussianBeam

### `__init__(self, bmaj=15 arcsec, bmin=None, bpa=0 deg, truncate=4.0)`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `bmaj` | `Quantity[arcsec]` | 长轴 FWHM（默认：15"） |
| `bmin` | `Quantity[arcsec] \| None` | 短轴 FWHM。`None` 表示圆形波束 |
| `bpa` | `Quantity[deg]` | 位置角，从北向东量（默认：0°） |
| `truncate` | `float` | 截断半径，以 FWHM 为单位（默认：4.0） |

> 面积 = π × bmaj × bmin / (4 × ln2)。波束在每个像素上使用 RectBivariateSpline 做积分。

### `needs_pad(self)` → `tuple[int, int]`

返回 (pad_x, pad_y) = 核尺寸的一半。在卷积前使用。

---

> **自定义波束：** 从 `_BaseBeam` 继承，实现 `f_kernel()`（返回可调用对象）、`kernel_size_px()`（返回半尺寸），并可选实现 `init_beam_header()`。
