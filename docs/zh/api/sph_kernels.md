# astrovodka.sph_kernels — SPH 平滑核

把 SPH 粒子平滑到像素网格上。所有公开核都带**自适应回退**：当像素尺寸与平滑长度之比过大时，会自动回退到 `DiracDeltaKernel`，再回退到 `GaussianKernel(truncate=6)`。

实现来自 MARTINI：从 `martini.sph_kernels` 导入并在此转发。

## 公开核类

| 核类 | min_valid_size | 典型模拟 |
|-------|---------------|-------------------|
| `CubicSplineKernel()` | 1.16 px | Gadget/Gadget2（经典 M4） |
| `WendlandC2Kernel()` | 1.51 px | EAGLE |
| `WendlandC6Kernel()` | 1.29 px | Magneticum |
| `QuarticSplineKernel()` | 1.24 px | Colibre（SPHENIX） |
| `GaussianKernel(truncate=3.0)` | 2.34 px | 通用 |
| `DiracDeltaKernel()` | ≤0.5 px | 点状粒子 |

> **min_valid_size**：使核积分精度达到 1% 所需的最小平滑长度（以像素计）。低于该值时核会回退。
>
> 若需要不带回退的原始核，请使用**下划线前缀**版本（例如 `_CubicSplineKernel`）。
>
> 设置 `Martini.insert_source_in_cube(skip_validation=True)` 可跳过精度校验。

### `find_fwhm(f: Callable) -> float`

用 `scipy.optimize.fsolve` 求解核函数的 FWHM。假定极大值出现在 f(0) 处且函数对称。
