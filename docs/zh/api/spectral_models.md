# astrovodka.spectral_models — 21-cm 谱线模型

对每个 SPH 粒子的 21-cm 谱线轮廓建模。

这些类继承 MARTINI 的 `GaussianSpectrum` 与 `DiracDeltaSpectrum`；唯一新增的是模块文档中说明的光度距离流量标度。

## GaussianSpectrum

### `__init__(self, sigma=7 km/s, ncpu=None, spec_dtype=np.float64)`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `sigma` | `Quantity[km/s] \| str` | 谱线宽度。`'thermal'` 表示逐粒子取 √(k<sub>B</sub>T/m<sub>p</sub>) |
| `ncpu` | `int \| None` | 并行 CPU 数（默认：1） |
| `spec_dtype` | `type` | 数据类型（默认：float64） |

> 按 `0.5 × [erf((b-v)/(σ√2)) - erf((a-v)/(σ√2))]` 求值，其中 (a,b) 为通道边界。

---

## DiracDeltaSpectrum

### `__init__(self, ncpu=None, spec_dtype=np.float64)`

δ 函数模型。全部流量都进入包含该粒子速度的那一个通道。

---

> **流量如何计算：** 每个粒子的谱会在每个通道上做积分。无量纲的结果乘以 `mHI / D_L² / dv`，再按 M<sub>HI</sub>/M<sub>sun</sub> = 2.36×10⁵ × (D<sub>L</sub>/Mpc)² × (S<sub>21</sub>/Jy km/s) 换算为 Jy。
>
> **自定义谱模型：** 从 `_BaseSpectrum` 继承。实现 `half_width(source)` 与 `spectral_function(a, b, vmids, extra_data)`。
