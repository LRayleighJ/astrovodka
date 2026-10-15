# The HI/H₂ models

astrovodka assigns a neutral and molecular hydrogen mass to every TNG gas cell
using one of four published prescriptions. This page explains what they compute,
how they differ, how to choose between them, and where their assumptions bite.

Signatures, defaults and the full configuration table live in the
[API reference](api/himodel.md); this page is about the physics and the
practical consequences.

- [Why not MARTINI's partition](#why-not-martinis-partition)
- [What the models compute](#what-the-models-compute)
- [The 2×2: GD14 and K13, volume and projected](#the-22-gd14-and-k13-volume-and-projected)
- [The UV field](#the-uv-field)
- [Choosing a model](#choosing-a-model)
- [Tuning the configuration](#tuning-the-configuration)
- [What to watch out for](#what-to-watch-out-for)

---

## Why not MARTINI's partition

MARTINI splits gas into HI and H₂ using a pressure-based relation calibrated on
nearby galaxies. It has no redshift dependence, so at high redshift — where
galaxies are denser, more metal-poor and bathed in a stronger UV field — it has
no reason to be right.

The Diemer+2018 models replace it with a calculation from local conditions:
gas density, metallicity, and the interstellar radiation field. That is what
makes them usable across the redshift range astrovodka targets.

## What the models compute

Every model produces, per gas cell, a **molecular fraction**

```
f_mol = M_H2 / (M_HI + M_H2)
```

and then splits the cell's neutral gas accordingly:

```
M_H2      = f_mol       × M_neutral
M_HI      = (1 - f_mol) × M_neutral
M_neutral = the cell's total neutral hydrogen mass
```

`f_mol` rises with gas density and metallicity — denser gas shields itself from
dissociating radiation more effectively, and dust is what does the shielding, so
metallicity enters through the dust-to-gas ratio. It falls with the strength of
the UV field, which is what breaks H₂ apart.

Every cell also gets a **neutral fraction**, which accounts for ionised gas;
`compute_hi_mass` returns `m_neutral` as `m_hi + m_h2` by construction.

## The 2×2: GD14 and K13, volume and projected

The four models are the cross product of two prescription families and two ways
of estimating the gas surface density the prescriptions need.

| | **Volume** | **Projected** |
|---|---|---|
| **Gnedin & Draine 2014** | `GD14_vol` | `GD14_map` |
| **Krumholz+ 2013** | `K13_vol` *(default)* | `K13_map` |

### GD14 vs K13

Both families predict `f_mol` from a radiation field parameter and a
dust-to-gas-weighted surface density, but they were calibrated on different
simulations and treat the shielding differently. In practice, on TNG galaxies,
**K13 tends to produce much less H₂ than GD14** — on subhalo 122 the ratio is
about 0.18 (`K13_vol`) versus 0.43 (`GD14_vol`). That is a real physical
disagreement between the prescriptions, not a bug, and it is the main reason
for offering both.

### Volume vs projected

- The **volume** models estimate each cell's surface density from its own
  density and a Jeans-length approximation, `Σ = λ_J · ρ`. Purely local, no
  projection, no stellar particles needed.
- The **projected** models build a face-on map of the galaxy and read the
  surface density off it. This is closer to how the observations are actually
  made — a telescope measures a column, not a volume density — and it is what
  Diemer+2018 use for the projected variants.

The projected models therefore need the **stellar particles** as well as the gas
(to compute the stellar surface density) and cost more. They also have to make a
decision about gas that falls outside the projection box — see
[What to watch out for](#what-to-watch-out-for).

## The UV field

`f_mol` depends on how strong the dissociating radiation field is at each gas
cell, so a model is only as good as its UV field.

The **volume** models take theirs from one of two sources, chosen with
`uv_mode`:

| `uv_mode` | What it does |
|---|---|
| `'adaptive_fft'` *(default)* | Convolves the galaxy's star-formation-rate distribution with a 1/r² kernel on a 3D FFT grid, giving a per-cell field that responds to where the stars actually are |
| `'global'` | Uses the uniform Faucher-Giguère+2011 extragalactic background everywhere |

Both apply the background as a **floor**, so `'global'` is the limit of
`'adaptive_fft'` for a galaxy with no star formation. In a star-forming galaxy
the two differ substantially: on the tutorial's subhalo 113, switching from the
background to the adaptive field changes the HI mass by about **+14%**, because
the strong local UV field pushes gas out of the molecular phase.

The **projected** models derive their UV field from the SFR map as part of the
method definition, so `uv_mode` does not apply to them. Passing `'global'`
alongside a `*_map` model raises rather than silently ignoring the argument.

> **A bug worth knowing about.** Before 0.3.0 the adaptive field did not work:
> the 3D interpolator was built with one axis' coordinates reused for all three,
> so nearly every particle fell outside its bounds and got `fill_value=0.0`.
> The result was a field that was numerically identical to the uniform
> background, meaning `uv_mode='adaptive_fft'` and `uv_mode='global'` gave the
> same masses and molecular masses were wrong by up to 42%. The test that guards
> this is placement-based rather than a pinned number — see the
> [verification record](verification.md).

## Choosing a model

The default is `K13_vol`, and it is a reasonable default: it needs no stellar
particles, is the cheapest of the four, and is the volume variant that
Diemer+2018 use as their baseline.

| If you want to… | Use |
|---|---|
| Compare against the Diemer+2018 catalog | `K13_vol` |
| Match how observations are actually made (a column through the galaxy) | `K13_map` or `GD14_map` |
| See the effect of a different prescription | switch between `K13_*` and `GD14_*` |
| Run with the least input data | a `*_vol` model — no stellar particles required |

```python
from astrovodka.sources import TNGSource

source = TNGSource(basePath=..., snapNum=99, subID=113, hi_model="K13_map")
print(source.m_HI_total, source.m_H2_g.sum(), source.f_mol_g.mean())
```

Because the models disagree by design, **report which one you used**. A single
number labelled "the HI mass" is not reproducible without it.

## Tuning the configuration

`HIModelConfig` is a frozen dataclass holding the snapshot-dependent values
(`h`, `z`, box size) and the numerical settings. It can be passed as an object
or as a plain dict:

```python
from astrovodka.himodel import HIModelConfig

TNGSource(..., hi_config=HIModelConfig(h=0.6774, z=0.0, ngrid=64))
TNGSource(..., hi_config={"ngrid": 64})          # same thing, merged for you
```

`h` and `z` always come from the snapshot, never from the caller — passing them
explicitly to `TNGSource` has no effect, by design, because a mismatch between
the snapshot's cosmology and the model's would be silent.

The settings most likely to matter:

| Field | Default | When to change it |
|---|---|---|
| `ngrid` | 128 | Lower it for speed while exploring; the UV grid becomes coarse |
| `f_esc` | 0.1 | If you want to test a different UV escape fraction |
| `box_expand` | 1.1 | How far the adaptive UV box extends beyond the gas |
| `map_clip_outside` | `False` | See below |
| `k13_variant` | `'2p'` | Which prescription sets `n_CNM` in the projected K13 model |
| `vol_fine_grid` | `0.0` | Diemer+2018 Appendix A.2's second, finer UV grid — **off by default**, see below |

## What to watch out for

### Gas outside the projection box

The projected models define the galaxy *as* their projection box. Cells outside
it get **zero** HI and H₂ — that is the method, and the catalog total for the
method is the in-box total.

The fraction of gas mass that falls inside is reported as `f_gas_in_map`, and
`TNGSource` warns when it drops below 0.9. For the tutorial's subhalo 113 it is
0.955; for subhalo 122 it is 0.883, which is low enough to be worth noticing.

If you need mass conservation instead, clamp the outside gas onto the edge pixel:

```python
HIModelConfig(h=0.6774, z=0.0, map_clip_outside=True)
```

This conserves mass but piles it on the boundary, which is not obviously better
— the paper does not specify a convention here, so both options are defensible
and you should say which you used.

### The second UV grid is off by default

Diemer+2018 Appendix A.2 runs the volume calculation twice, the second time on a
grid spanning two gas half-mass radii, because a single coarse grid can
underestimate the UV field at the centre of a galaxy. `vol_fine_grid` implements
that switch and measurably moves `f_mol` by −3 to −5 percentage points.

It is **off by default** for two reasons: the paper does not say how to combine
the two grids, and turning it on breaks the bit-exactness with the reference
implementation that the [verification record](verification.md) documents.
Enabling it is an experiment, not a correction.

### Memory

The volume models build a full `ngrid³` grid and the projected models a
`ngrid_map²` map. At the default 128 these are modest, but subhalos with more
than ~10⁷ gas cells would benefit from chunking, which is not implemented.

### The physics layer stands alone

`astrovodka.himodel` imports nothing from the observation layer, so you can use
the models on a gas cell array from anywhere — another simulation, a
post-processing pipeline, an interactive session — with no data cube involved:

```python
from astrovodka.himodel import HIModelConfig, compute_hi_mass

result = compute_hi_mass(gas, model="K13_vol",
                         cfg=HIModelConfig(h=0.6774, z=0.0))
result.m_hi_msun, result.m_h2_msun
```
