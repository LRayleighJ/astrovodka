# astrovodka.sources.sph_source — SPHSource

从 SPH 粒子数组抽象出发射 HI 的天体。同时支持旧版距离模式和新的宇宙学红移模式。

本类继承 `martini.sources.SPHSource`。粒子机制来自 MARTINI，只有红移距离模式是 astrovodka 的。

## SPHSource

### `__init__(self, *, distance=None, redshift=None, cosmology=None, vpeculiar=0 km/s, rotation=None, L_coords=None, ra=0 deg, dec=0 deg, h=0.7, T_g=None, mHI_g=None, xyz_g=None, vxyz_g=None, hsm_g=None, coordinate_axis=None, coordinate_frame=ICRS())`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `distance` | `Quantity[Mpc] \| None` | 旧版模式：物理距离（与 `h` 配合使用） |
| `redshift` | `float \| None` | 宇宙学红移模式（与 `cosmology` 配合使用） |
| `cosmology` | `Cosmology \| None` | Astropy 宇宙学（默认：Planck15） |
| `vpeculiar` | `Quantity[km/s]` | 本动速度（默认：0） |
| `rotation` | `Rotation \| None` | scipy 的 Rotation |
| `L_coords` | `L_coords \| None` | 角动量对齐 |
| `ra` / `dec` | `Quantity[deg]` | 源的中心（默认：0） |
| `h` | `float` | 哈勃常数（旧版模式，默认：0.7） |
| `T_g` | `Quantity[K] \| None` | 粒子温度 |
| `mHI_g` | `Quantity[M_sun]` | 每个粒子的 HI 质量 |
| `xyz_g` | `Quantity[kpc]` | 位置。形状为 (3,N) 或 (N,3)。x 为视线方向 |
| `vxyz_g` | `Quantity[km/s]` | 速度。必须与 xyz_g 对应 |
| `hsm_g` | `Quantity[kpc] \| None` | SPH 平滑长度（**FWHM**） |
| `coordinate_axis` | `int \| None` | 自动检测 |
| `coordinate_frame` | `BaseRADecFrame` | 坐标系（默认：ICRS） |

> **模式：**
> - **红移模式**（新增）：设置 `redshift`，可选设置 `cosmology`。粒子按**角直径距离**
>   D<sub>A</sub> 放置到天球上（因此在任意红移下角尺寸都正确）；流量定标使用
>   **光度距离** D<sub>L</sub>；并且 `vsys = c·z/(1+z)`，即与 WCS 光谱轴一致的
>   射电约定。
> - **距离模式**（旧版）：设置 `distance` 与 `h`。计算 vsys=h×100×D + vpeculiar。
>   此时 D<sub>A</sub> == D<sub>L</sub> == `distance`，也就是未作改动的 MARTINI 行为。
> - 两者不能同时设置（会抛出 ValueError）。
> - 平滑长度必须是 FWHM（而不是核半径）。

---

## 关键属性

| 属性 | 类型 | 说明 |
|-----------|------|-------------|
| `vsys` | `Quantity[km/s]` | 系统速度（射电约定）。可用作 `DataCube.spectral_centre` |
| `distance` | `Quantity[Mpc]` | 光度距离 D<sub>L</sub> |
| `luminosity_distance` | `Quantity[Mpc]` | `distance` 的别名；流量归一实际使用的就是它 |
| `angular_diameter_distance` | `Quantity[Mpc]` | D<sub>A</sub>；天球放置实际使用的就是它 |
| `vhubble` | `Quantity[km/s]` | 系统哈勃退行速度（宇宙学模式下为 `c·z/(1+z)`） |
| `npart` | `int` | 裁剪后的粒子数 |
| `input_mass` | `Quantity[M_sun]` | 裁剪前的 HI 总质量 |

## 方法

### `rotate(self, rotation=None, *, L_coords=None)`

施加旋转。返回 `Rotation`。

### `translate(self, translation_vector: Quantity[kpc])`

平移位置。

### `boost(self, boost_vector: Quantity[km/s])`

叠加速度偏移。

### `preview(self, max_points=5000, fig=1, lim=None, vlim=None, point_scaling='auto', title='', save=None)`

三联图：天球投影 + 2 幅 PV 图。返回 `Figure`。

### `save_current_rotation(self, fname: str)`

把 3×3 旋转矩阵保存为文本。

### `save_current_affine_transformations(self, fname: str)`

把 (2,4,4) 仿射变换保存为 `.npy`。

### `load_affine_transformations(self, fname: str)`

从 `.npy` 读入并施加仿射变换。
