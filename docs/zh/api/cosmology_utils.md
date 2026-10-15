# astrovodka.cosmology_utils — 宇宙学工具

宇宙学距离（默认 Planck15）与 TNG 代码单位换算的辅助函数。

## 距离函数

### `get_cosmology(cosmo=None)` → `Cosmology`

返回宇宙学实例（默认 Planck15）。参数非法时抛出 TypeError。

### `luminosity_distance(redshift, cosmo=None)` → `Quantity[Mpc]`

D<sub>L</sub> — 用于流量：flux ∝ M<sub>HI</sub> / D<sub>L</sub>²。当 z < 1e-6 时改用哈勃定律下限 cz/H₀。

### `angular_diameter_distance(redshift, cosmo=None)` → `Quantity[Mpc]`

D<sub>A</sub> — 用于角尺度：θ = arctan(size / D<sub>A</sub>)。D<sub>L</sub>/D<sub>A</sub> = (1+z)²。

### `comoving_distance(redshift, cosmo=None)` → `Quantity[Mpc]`

视线方向的共动距离。

### `hubble_parameter(redshift, cosmo=None)` → `Quantity[km/s/Mpc]`

H(z) — 红移 z 处的哈勃参数。

### `systemic_velocity_from_redshift(redshift)` → `Quantity[km/s]`

狭义相对论（运动学/光学）多普勒公式：v/c = ((1+z)²-1)/((1+z)²+1)。

> **注意：** 这不是本流水线采用的约定。`SPHSource.vsys` 来自
> `radio_velocity_from_redshift`，后者与 WCS 光谱轴一致。

### `radio_velocity_from_redshift(redshift, rest_freq=None)` → `Quantity[km/s]`

射电约定的退行速度，v = c(1 − f<sub>obs</sub>/f<sub>rest</sub>) = c·z/(1+z)。
这是 `SpectralCoord(doppler_convention="radio")`、FITS `VRAD` 轴和 `doppler_radio` 所采用的
约定，因此系统速度必须用它才能与数据立方体的光谱轴对齐。与 `H(z)·D_C` 不同，它永远不会
超过 c。默认 rest_freq = HI 1.42 GHz。

### `observed_frequency(redshift, rest_freq=None)` → `Quantity[Hz]`

f<sub>obs</sub> = f<sub>rest</sub>/(1+z)。默认 rest_freq = HI 1.42 GHz。

---

## 单位换算函数

### `rho_code_to_physical(rho_code, redshift, hubble_h)` → `ndarray [g/cm³]`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `rho_code` | `ndarray` | TNG 密度 [1e10 M_sun/kpc³ × h²] |
| `redshift` | `float` | 红移 |
| `hubble_h` | `float` | h = H₀/100 |

ρ<sub>phys</sub> = ρ<sub>com</sub> × (1+z)³。

### `rho_code_to_nH(rho_code, redshift, hubble_h, X_H=0.76)` → `ndarray [cm⁻³]`

| 名称 | 类型 | 说明 |
|------|------|-------------|
| `X_H` | `float` | 氢的质量分数（默认：0.76） |

n<sub>H</sub> = X<sub>H</sub> × ρ<sub>phys</sub> / m<sub>p</sub>。

### `calc_Jeans_length_physical(u_code, rho_code, redshift, hubble_h)` → `ndarray [kpc]`

λ<sub>J</sub> = √(γ(γ-1)u / (Gρ))。

### `calc_Sigma_neutral_physical(u_code, rho_code, f_neutral_XH, redshift, hubble_h)` → `ndarray [M_sun/kpc²]`

Σ = ρ<sub>phys</sub> × λ<sub>J</sub> × f<sub>neutral</sub> × X<sub>H</sub>。体密度 HI 模型（`GD14_vol`、`K13_vol`）会用到它。
