# astrovodka.cosmology_utils — Cosmology Utilities

Helper functions for cosmological distances (Planck15 default) and TNG code unit conversions.

## Distance Functions

### `get_cosmology(cosmo=None)` → `Cosmology`

Returns cosmology instance (default Planck15). Raises TypeError if invalid.

### `luminosity_distance(redshift, cosmo=None)` → `Quantity[Mpc]`

D<sub>L</sub> — used for flux: flux ∝ M<sub>HI</sub> / D<sub>L</sub>². For z < 1e-6, applies Hubble-law floor cz/H₀.

### `angular_diameter_distance(redshift, cosmo=None)` → `Quantity[Mpc]`

D<sub>A</sub> — used for angular scale: θ = arctan(size / D<sub>A</sub>). D<sub>L</sub>/D<sub>A</sub> = (1+z)².

### `comoving_distance(redshift, cosmo=None)` → `Quantity[Mpc]`

Line-of-sight comoving distance.

### `hubble_parameter(redshift, cosmo=None)` → `Quantity[km/s/Mpc]`

H(z) — Hubble parameter at redshift z.

### `systemic_velocity_from_redshift(redshift)` → `Quantity[km/s]`

Relativistic Doppler: v/c = ((1+z)²-1)/((1+z)²+1).

> **Note:** SPHSource uses H(z)×D<sub>L</sub>/(1+z) for `vsys` consistency, NOT this function.

### `observed_frequency(redshift, rest_freq=None)` → `Quantity[Hz]`

f<sub>obs</sub> = f<sub>rest</sub>/(1+z). Default rest_freq = HI 1.42 GHz.

---

## Unit Conversion Functions

### `rho_code_to_physical(rho_code, redshift, hubble_h)` → `ndarray [g/cm³]`

| Name | Type | Description |
|------|------|-------------|
| `rho_code` | `ndarray` | TNG density [1e10 M_sun/kpc³ × h²] |
| `redshift` | `float` | Redshift |
| `hubble_h` | `float` | h = H₀/100 |

ρ<sub>phys</sub> = ρ<sub>com</sub> × (1+z)³.

### `rho_code_to_nH(rho_code, redshift, hubble_h, X_H=0.76)` → `ndarray [cm⁻³]`

| Name | Type | Description |
|------|------|-------------|
| `X_H` | `float` | Hydrogen mass fraction (default: 0.76) |

n<sub>H</sub> = X<sub>H</sub> × ρ<sub>phys</sub> / m<sub>p</sub>.

### `calc_Jeans_length_physical(u_code, rho_code, redshift, hubble_h)` → `ndarray [kpc]`

λ<sub>J</sub> = √(γ(γ-1)u / (Gρ)).

### `calc_Sigma_neutral_physical(u_code, rho_code, f_neutral_XH, redshift, hubble_h)` → `ndarray [M_sun/kpc²]`

Σ = ρ<sub>phys</sub> × λ<sub>J</sub> × f<sub>neutral</sub> × X<sub>H</sub>. Used in K13 model.
