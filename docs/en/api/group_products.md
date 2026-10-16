# astrovodka.group_products — Group/Cluster-Scale 3D Products

Grids a whole **FOF group** — every subhalo it contains plus the group's unbound
fuzz — and classifies the gas into ISM, CGM and IGM. Where
[data_products](data_products.md) shows one galaxy, this shows the group that
galaxy belongs to.

| Name | Value | Meaning |
|------|-------|---------|
| `ISM` | `0` | Bound to a subhalo and within `k x R_ISM` of its centre |
| `CGM` | `1` | Bound to a subhalo but outside that radius (no R<sub>200c</sub> cut) |
| `IGM` | `2` | In the group but bound to no subhalo (the fuzz) |

`CLASS_LABELS` maps those integers back to `"ISM"`, `"CGM"`, `"IGM"`.

---

## Membership

Cell ownership comes from the group catalog offsets, not from a spatial guess:

```
gas cells of subhalo s  =  [Subhalo/SnapByType[s, 0], + SubhaloLenType[s, 0])
cells of the group      =  [Group/SnapByType[g, 0],  + GroupLenType[g, 0])
fuzz                    =  the group's cells minus every subhalo's
```

> **Lengths must come from `SubhaloLenType`, never from differencing the offsets
> array.** `SnapByType` is a cumulative start-offset table, and for the *last*
> subhalo of a group its next row already belongs to the next group, so the
> difference silently swallows that group's fuzz — on TNG-100 snapshot 99 group
> 59 the last subhalo has `SubhaloLenType == 0` yet a 319,249-cell step.
> `from_group` sums the mass of every block and compares it with
> `SubhaloMassType`, refusing to continue if they disagree; that check is what
> makes the membership trustworthy rather than assumed.

The fuzz is a real component: 1.5% of the gas mass in the most massive TNG-100
group, ~10% at the mass of group 59, ~24% in smaller groups.

---

## `ism_radius(rad_gas, rad_star, k_gas=1.0, k_star=2.0)` → `ndarray`

`min(k_gas x R_half,gas, k_star x R_half,star)`, from
`SubhaloHalfmassRadType` columns 0 and 4, in the same unit as the inputs
(ckpc/h).

| Name | Type | Description |
|------|------|-------------|
| `rad_gas`, `rad_star` | `array-like` | Gas and stellar half-mass radii [ckpc/h] |
| `k_gas`, `k_star` | `float` | Multipliers (defaults 1 and 2) |

> The `min` is what keeps the definition usable at group scale. A group's
> central subhalo has an `R_half,gas` of hundreds of ckpc/h because it contains
> the whole hot halo — 263 ckpc/h for TNG-100 group 59, against 19 ckpc/h for
> `R_half,star` — so a gas-only radius would swallow most of the group.
> A non-positive or non-finite radius counts as "no constraint"; if neither is
> usable the result is 0 and that subhalo contributes no ISM.

---

## `classify_ism_cgm_igm(owner, r_gal, r_ism, nH=None, T=None, nH_max=None, T_max=None)` → `ndarray[int8]`

| Name | Type | Description |
|------|------|-------------|
| `owner` | `ndarray of int` | Owning subhalo per cell, `-1` for the fuzz |
| `r_gal` | `ndarray` | Distance to the owning subhalo's centre [ckpc/h]; ignored where `owner < 0` |
| `r_ism` | `ndarray or float` | ISM radius per cell (already indexed by `owner`), or a scalar |
| `nH` | `ndarray` | H number density [cm⁻³]; required with `nH_max` |
| `T` | `ndarray` | Temperature [K]; required with `T_max` |
| `nH_max`, `T_max` | `float` | Optional phase caps: a cell failing either is demoted out of the ISM |

> The caps are **off by default**, and the measurement is the reason: the
> population they target — dense *and* hot shocked gas, `n_H > 0.1 cm^-3` with
> `T > 1e6 K` — is 0.000 of the ISM mass in group 59, while a plain `n_H < 100`
> cap removes 0.6%.

---

## `load_group_members(basePath, snapNum, groupID)` → `dict`

The catalog records of every subhalo in the group, from the contiguous rows
`[GroupFirstSub, + GroupNsubs)` of the subhalo table.

Returns `group` (the FOF group's own record), `first`, `pos` (n, 3) [ckpc/h],
`len_type` (n, 6), `mass_type` (n, 6) [1e10 M<sub>sun</sub>/h], `mass`, `rhalf_type`
(n, 6) [ckpc/h], `sfr`, `flag`, and `start` (n, 6) — each subhalo's start offset
in the snapshot array, per particle type.

---

## `GroupDataProducts`

### `__init__(self, ngrid=256)`

### `from_group(cls, basePath, snapNum, groupID, *, ngrid=256, hi_model='K13_vol', hi_config=None, uv_mode='adaptive_fft', compute_hi=True, compute_rsd=False, k_gas=1.0, k_star=2.0, nH_max=None, T_max=None, igm_uv='group', ncpu=1, cosmology=None, box_halfsize_ckpc=None)` → `GroupDataProducts`

| Name | Type | Description |
|------|------|-------------|
| `basePath` | `str` | TNG output directory (with `groups_<snapNum>/` and the `postprocessing/offsets` files) |
| `snapNum` | `int` | Snapshot number |
| `groupID` | `int` | FOF group ID (`haloID`) |
| `ngrid` | `int` | Grid edge length (default 256; 128 for the most massive groups) |
| `hi_model` | `str` | `'K13_vol'` (default) or `'GD14_vol'`; the projected models are rejected |
| `hi_config` | `HIModelConfig \| dict \| None` | Model settings; `h`, `z` and box size come from the snapshot |
| `uv_mode` | `str` | `'adaptive_fft'` (per subhalo) or `'global'` |
| `compute_hi` | `bool` | `False` skips the HI/H2 model entirely |
| `compute_rsd` | `bool` | Also build the redshift-space grids (off by default, see below) |
| `k_gas`, `k_star` | `float` | ISM radius multipliers |
| `nH_max`, `T_max` | `float \| None` | Optional ISM phase caps |
| `igm_uv` | `str` | UV field for the fuzz: `'group'` (default) or `'uvb'` |
| `ncpu` | `int` | Processes for the per-subhalo HI loop (needs `fork`) |
| `cosmology` | `Cosmology \| None` | For the RSD shift; default Planck15 |
| `box_halfsize_ckpc` | `float \| None` | Pin the box; otherwise 1.05x the particle extent |

> **HI is computed per subhalo**, each in its own adaptive FFT box — the
> prescription the volumetric models were calibrated under. Running one FFT for
> the whole group would smear the field over ~30 kpc cells at cluster scale.
> Only gas-bearing subhalos are modelled (68 of 1710 members in group 59), so
> the loop is far cheaper than the subhalo count suggests. The fuzz is modelled
> separately with the group-wide UV field.
>
> **Redshift space is opt-in.** At z=0 one km/s of line-of-sight velocity is
> 10 ckpc/h, so a group whose members span a few hundred km/s occupies several
> Mpc along the line of sight: group 59 needs a 17 cMpc axis against 1.9 cMpc
> in the sky plane. The shift is applied after subtracting the group's bulk
> velocity, which is the frame the mock observations are made in; without that
> subtraction the map comes out empty (and `SubhaloDataProducts`' redshift-space
> grids do lose ~99.5% of their mass for exactly that reason).

### Real space first

**Classification happens in real space and is then carried into redshift space;
it is never computed there.** `classify_ism_cgm_igm` runs on the real-space
coordinates, and the redshift-space grids are built by moving those same cells
along the line of sight with their real-space `class_id` as the mask. The
consequences, all verified on TNG-100 group 59:

- `cells`, the class masses and every real-space grid are **bit-identical
  whether or not `compute_rsd` is set** — the real-space products are canonical.
- Each class keeps its mass exactly through the transform: real-space and
  redshift-space grid masses agree at a ratio of 1.000000 for ISM, CGM and IGM,
  and both match the metadata totals.
- Redshift space gets its **own grid box**, longer along the line of sight
  (`grid_edges_rsd`, plus `box_halfsize_ckpc_rsd` and `box_bounds_ckpc_h_rsd` in
  the metadata), so RSD never degrades the real-space resolution. The two boxes
  have different cell volumes, which `f_mass_in_rsd_box_gas` accounts for.

---

### `save_hdf5(self, filename)` / `load_hdf5(cls, filename)`

Grids (LZF-compressed, chunked), the per-cell arrays, the subhalo table and the
metadata. `None` metadata values are stored as the string `"none"`, and scalar
subhalo entries come back as length-1 arrays, since HDF5 scalar datasets cannot
carry the compression the rest of the file uses.

---

## Grid attributes

All grids are `(ngrid, ngrid, ngrid)` float32 on one shared box, with `x` along
the line of sight, in M<sub>sun</sub>/(ckpc/h)³:

| Attribute | Description |
|-----------|-------------|
| `dm_real`, `dm_rsd` | Total dark matter |
| `dm_igm_real` | Dark matter bound to no subhalo |
| `gas_real`, `gas_rsd` | Total gas |
| `gas_ism_real`, `gas_cgm_real`, `gas_igm_real` | Gas by class |
| `hi_real`, `hi_rsd` | Neutral hydrogen |
| `hi_ism_real`, `hi_ism_rsd`, `hi_cgm_real`, `hi_cgm_rsd`, `hi_igm_real`, `hi_igm_rsd` | HI by class |
| `h2_real` | Molecular hydrogen |
| `t_mw`, `nh_mw` | Gas-mass-weighted temperature [K] and n<sub>H</sub> [cm⁻³] |
| `grid_edges` | Tuple of 3 edge arrays [ckpc/h] for the real-space box |
| `grid_edges_rsd` | Same, for the redshift-space box (longer along the line of sight); `None` when `compute_rsd` is off |

## `cells` — per-cell arrays

`coords_ckpc_h`, `m_gas_msun`, `m_hi_msun`, `m_h2_msun`, `nH_cm3`, `T_K`,
`r_gal_ckpc_h`, `r_grp_over_r200c`, `owner`, `class_id`.

> These exist so that a re-classification, a radial profile or a zoomed re-grid
> does not need the snapshot again — a group-wide grid cannot resolve a galaxy,
> so the zoomed figures are built from these by re-histogramming a sub-volume.

## `metadata`

Group identity (`groupID`, `snapNum`, `redshift`, `hubble_h`, `M200c_msun`,
`R200c_ckpc_h`, `R200c_kpc`, `group_pos_ckpc_h`, `Nsubs`, `N_sub_with_gas`),
counts (`n_gas`, `n_dm`, `n_cells_ism/cgm/igm`), the mass budget
(`m_gas_*_msun`, `m_hi_*_msun`, `m_h2_total_msun`, `m_dm_total_msun`,
`m_dm_igm_total_msun`), the classification settings (`ism_k_gas`, `ism_k_star`,
`ism_nH_max`, `ism_T_max`), the box (`box_halfsize_ckpc`,
`box_bounds_ckpc_h`), the model (`hi_model`, `uv_mode`, `igm_uv`,
`vol_fine_grid`) and, when RSD is on, `f_mass_in_rsd_box_gas`.

## `subhalos`

The `load_group_members` table, so that per-subhalo quantities can be joined to
the grids without re-reading the catalog.
