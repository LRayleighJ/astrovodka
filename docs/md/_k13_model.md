# astrovodka.sources._k13_model — K13 HI/H2 Model

Implements the Krumholz 2013 (Diemer+2018) volumetric HI/H2 model.

### `calc_UMW_global_Faucher11(redshift)` → `float [Draine units]`

Global UV background from Faucher-Giguère+2011. 4th-order polynomial:
log<sub>10</sub>U<sub>MW</sub> = -0.004292z⁴ + 0.061446z³ - 0.358413z² + 1.204787z - 2.846354.

At z=0, U<sub>MW</sub> ≈ 1.4×10⁻³.

---

### `get_f_neutral_H(nH_code, sfr_code)` → `ndarray`

| Name | Type | Description |
|------|------|-------------|
| `nH_code` | `ndarray` | TNG NeutralHydrogenAbundance (f_HI only) |
| `sfr_code` | `ndarray` | TNG StarFormationRate |

Estimates f<sub>neutral</sub> = f<sub>HI</sub> + f<sub>H2</sub>:
- SF particles (SFR>0): f<sub>neutral</sub> = 1
- Non-SF: f<sub>neutral</sub> ≈ f<sub>HI</sub>

---

### `f_mol_K13_vol(rho_code, u_code, Z, nH_code, sfr_code, X_H, U_MW, redshift=0.0, hubble_h=0.6774, mgas=None)` → `(f_mol, f_HI, f_H2)`

| Name | Type | Description |
|------|------|-------------|
| `rho_code` | `ndarray` | TNG comoving density |
| `u_code` | `ndarray` | TNG internal energy [(km/s)²] |
| `Z` | `ndarray` | Metallicity (mass fraction) |
| `nH_code` | `ndarray` | NeutralHydrogenAbundance |
| `sfr_code` | `ndarray` | Star formation rate |
| `X_H` | `ndarray` | Hydrogen mass fraction |
| `U_MW` | `float \| ndarray` | UV radiation field (scalar or per-particle) |
| `redshift` | `float` | For unit conversions |
| `hubble_h` | `float` | h = H₀/100 |
| `mgas` | `ndarray \| None` | Gas masses (unused, interface compat) |

**Algorithm:**
1. f<sub>neutral</sub> from `get_f_neutral_H`
2. n<sub>CNM</sub> = max(n<sub>H</sub>×f<sub>neutral</sub>, n<sub>CNM,2p</sub>) — two-phase equilibrium (K13 Eq. C14)
3. χ = 7.2 × U<sub>MW</sub> × (n<sub>CNM</sub>/10)⁻¹ — radiation parameter
4. Σ<sub>neutral</sub> = ρ × λ<sub>J</sub> × f<sub>neutral</sub> × X<sub>H</sub> — surface density
5. τ<sub>c</sub> = 0.066 × f<sub>c</sub> × D<sub>MW</sub> × Σ — dust optical depth (f<sub>c</sub>=5)
6. s = ln(1+0.6χ+0.01χ²) / (0.6×τ<sub>c</sub>) — K13 Eq. 6
7. f<sub>mol</sub> = 1-3s/(4+s) if s<2 else 0 — K13 Eq. 2

---

### `compute_HI_mass_particle(gas_data, U_MW, redshift=0.0, hubble_h=0.6774)` → `(m_HI, m_H2) [code units: 1e10 M_sun/h]`

| Name | Type | Description |
|------|------|-------------|
| `gas_data` | `dict` | Keys: Masses, Density, InternalEnergy, SFR, NeutralHydrogenAbundance, GFM_Metals, GFM_Metallicity |

Convenience wrapper. Multiply by 1e10/h to get M<sub>sun</sub>. GFM_Metals handled for both 1D and 2D formats.

---

### `compute_galaxy_fft_uv(gas_data, redshift=0.0, hubble_h=0.6774, ngrid=128, f_esc=0.1, box_expand=1.1, return_grid=False)` → `ndarray | (ndarray, ndarray, list)`

| Name | Type | Description |
|------|------|-------------|
| `gas_data` | `dict` | Keys: Coordinates, StarFormationRate |
| `ngrid` | `int` | FFT grid per dimension (default: 128) |
| `f_esc` | `float` | UV escape fraction (default: 0.1) |
| `box_expand` | `float` | Adaptive box multiplier (default: 1.1) |
| `return_grid` | `bool` | If True, return `(per_particle, grid, edges)` |

Builds SFR grid → FFT-convolves with 1/r² kernel → interpolates to particles.
Always applies UVB floor. SFR→FUV: 3.3×10⁻⁶ photons/s/cm²/Hz per M<sub>sun</sub>/yr.
Draine normalization: 3.43×10⁻⁸ erg/s/cm².
