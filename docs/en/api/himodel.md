# `astrovodka.himodel` — HI/H2 mass models

Four selectable models for the neutral and molecular hydrogen content of
IllustrisTNG gas cells, all as implemented by
[Diemer+2018 (arXiv:1806.02341)](https://arxiv.org/abs/1806.02341):

| Method | Surface density from | Reference |
|--------|----------------------|-----------|
| `GD14_vol` | per-cell Jeans approximation, `Sigma = lambda_J * rho` | Gnedin & Draine (2014), Appendix C.3 |
| `GD14_map` | face-on 128² projection | same + §2.3 |
| `K13_vol` | per-cell Jeans approximation | Krumholz, McKee & Tumlinson (2013), Appendix C.4 |
| `K13_map` | face-on 128² projection | same + §2.3 |

All four replace MARTINI's pressure-based partition, which is calibrated on
nearby galaxies and has no redshift dependence.

## Quick start

```python
from astrovodka.himodel import HIModelConfig, compute_hi_mass

cfg = HIModelConfig(h=0.6774, z=0.0)        # h, z from the snapshot header

res = compute_hi_mass(gas, model='K13_vol', cfg=cfg)
res.m_hi          # per gas cell, TNG code units (1e10 Msun/h)
res.m_h2
res.m_hi_msun     # total, Msun

# the projected models need the stars and the subhalo record as well
res = compute_hi_mass(gas, model='K13_map', cfg=cfg, stars=stars, sub=sub)
```

## Selecting a model by keyword

`model` is case-insensitive. The bare model names `'K13'` and `'GD14'` mean the
volume variant, and the two-argument spelling `('GD14', 'map')` is accepted.

```python
from astrovodka.himodel import METHODS, resolve, is_known

METHODS           # ('GD14_vol', 'GD14_map', 'K13_vol', 'K13_map')
is_known('k13_map')          # True
resolve('K13')['label']      # 'K13_vol'
```

The same keyword is available on the source class:

```python
from astrovodka.sources import TNGSource

source = TNGSource(basePath=..., snapNum=99, subID=113, hi_model='GD14_map')
source.m_HI_total     # Msun
source.m_H2_g         # per-particle H2 mass
source.f_mol_g        # per-cell molecular fraction
```

## Functions

| Function | Description |
|----------|-------------|
| `compute_hi_mass(gas, model, cfg, *, stars, sub, uv_mode, uv)` | Per-cell HI and H2 masses under any of the four models. The main entry point. |
| `compute_galaxy_hi_h2_map(gas, stars, sub, cfg, model, ...)` | The projected path on its own: builds the face-on maps and returns the galaxy totals plus alignment diagnostics. |
| `calc_UMW_global_Faucher11(redshift)` | Faucher-Giguère+2011 UV background in Draine units. |
| `resolve(method)` / `is_known(method)` | Method-name lookup. |

Lower-level pieces, if you need them directly:

| Module | Contents |
|--------|----------|
| `himodel.physics` | Unit conversions, `f_mol_GD14_vol`, `f_mol_K13_vol`, `get_f_neutral_H` |
| `himodel.models` | `f_mol_GD14_map`, `f_mol_K13_map`, `calc_rho_star_map`, `calc_rho_dm_nfw_map`, `MODEL_SPEC` |
| `himodel.uvfield` | `compute_galaxy_fft_uv`, `compute_galaxy_fft_uv_map`, `compute_align_vector`, `deposit_gaussian_maps` |
| `himodel.constants` | Physical constants and model defaults (no `h`) |
| `himodel.config` | `HIModelConfig` |

## Units

| Quantity | Unit |
|----------|------|
| `HiMassResult.m_hi`, `.m_h2` | TNG code units, `1e10 Msun / h` — the same unit as the `Masses` field |
| `HiMassResult.m_hi_msun`, `.m_h2_msun`, `.m_neutral_msun` | Msun |
| `HiMassResult.u_mw` | Draine field units (`U_MW = 1` is the Milky Way) |

Densities inside the model are converted with `HIModelConfig.h` and `.z`;
`h` is a per-snapshot quantity and deliberately does not appear as a module
constant anywhere.

## Results

`compute_hi_mass` returns a `HiMassResult`:

| Attribute | Description |
|-----------|-------------|
| `m_hi`, `m_h2` | Per-cell masses, **full length** — cells the model excludes get zero, never a shorter array |
| `m_neutral` | `m_hi + m_h2`, exactly the neutral hydrogen mass |
| `f_mol`, `f_neutral` | Per-cell fractions |
| `u_mw` | Per-cell UV field (volume models) |
| `diagnostics` | UV mode, projection geometry, alignment, `f_gas_in_map` |
| `cells` | The underlying `Cells` table |

## The projected models and gas outside the box

The projected models build a face-on map sized
`L = 2 * max(2 R_half,gas, 6 R_half,*, 30 kpc)` and read each cell's molecular
fraction off it. Cells outside that box get **zero** HI and H2: the model's
galaxy *is* the projection box, and the catalog total for the method is the
in-box total.

The fraction of gas mass that falls inside is reported as `f_gas_in_map`;
`TNGSource` prints a warning when it drops below 0.9. To keep the outside gas
instead, clamp it onto the edge pixel with

```python
cfg = HIModelConfig(h=0.6774, z=0.0, map_clip_outside=True)
```

which conserves mass at the cost of piling it on the boundary.

## Configuration

`HIModelConfig` is a frozen dataclass; `HIModelConfig.build(dict_or_config, h=...,
z=...)` resolves user settings against a snapshot's header values, with `h`, `z`
and the box size always taken from the snapshot.

| Field | Default | Meaning |
|-------|---------|---------|
| `ngrid` | 128 | 3D FFT UV grid edge length |
| `f_esc` | 0.1 | UV escape fraction (Diemer+2018 §4.2) |
| `box_expand` | 1.1 | adaptive UV box = this × the gas coordinate extent |
| `ngrid_map` | 128 | 2D projection grid edge length |
| `vol_fine_grid` | 0.0 | second, higher-resolution UV box (paper App. A.2); off by default |
| `r_soft_map_mode` | `'dx'` | where the 2D kernel's softening length comes from |
| `map_clip_outside` | `False` | keep gas outside the projection box (mass conserving) |
| `k13_variant` | `'2p'` | which prescription sets `n_CNM` in the projected K13 model |
| `rho_dm_mode` | `'none'` | dark-matter mid-plane density for projected K13 |

See the class docstring for the full list and the reasoning behind each default.

## `uv_mode`

The volume models take their UV field either from a per-galaxy FFT convolution
of the star-formation distribution (`uv_mode='adaptive_fft'`, the default) or
from the uniform Faucher-Giguère+2011 background (`uv_mode='global'`). Both
apply the background as a floor, so `'global'` is the `'adaptive_fft'` limit of
a galaxy with no star formation.

The projected models derive their UV field from the SFR map — that is part of
the method definition — so `uv_mode` does not apply to them and passing
`'global'` alongside a `*_map` model raises.

> **Changed in v0.3.0.** The default was `'global'`, and `'adaptive_fft'` was
> silently identical to it: the interpolator was built with one axis' edge array
> reused for all three, so almost every particle fell outside its bounds and
> `fill_value=0.0` returned zero, leaving only the background. The default is
> now `'adaptive_fft'`, which does what its name says.

## Verification

The physical expressions are ported verbatim from `TNGHIcalc`, which is
verified bit-for-bit against Diemer+2018's reference implementation
(`TNGHI_GD14_model/gd14_model.py`) — see [the verification record](../verification.md)
for the chain of evidence and the reproduction commands.
