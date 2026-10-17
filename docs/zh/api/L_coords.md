# astrovodka.L_coords — 角动量坐标

一个 `NamedTuple`，用角动量来指定源的取向。

MARTINI 2.1.16 及以后提供同名的 `martini.L_coords`；在 2.1.9–2.1.15 上由 astrovodka 本地定义。

## `L_coords(incl=0 deg, pa=270 deg, az_rot=0 deg)`

| 字段 | 类型 | 默认值 | 说明 |
|-------|------|---------|-------------|
| `incl` | `Quantity[deg]` | 0° | 倾角（0° = 面向观测者，90° = 侧向） |
| `pa` | `Quantity[deg]` | 270° | 天球上的位置角（从北向东量，指退行侧） |
| `az_rot` | `Quantity[deg]` | 0° | 绕 L 轴的方位角旋转 |

**算法：** 先把中心 1/3 HI 粒子的角动量矢量对齐到 x 轴（面向观测者）。随后依次：绕 L 轴做 `az_rot` → 绕 y 轴做 `incl` 倾斜 → 绕 x 轴做 `pa` 天球旋转。所有旋转均为右手、外旋（extrinsic）。

**用法：** `source.rotate(L_coords=L_coords(incl=60*U.deg, pa=270*U.deg))`
