# Verification record

This page records how the HI/H₂ models in astrovodka were validated. The claim
being supported is narrow and specific: **the physics agrees, bit for bit, with
the reference implementation**, and the observation pipeline built on top of it
conserves mass and puts objects at the right distance.

Nothing here needs to be taken on trust — every check is either a test in
`tests/` or a command you can re-run.

- [What is being verified](#what-is-being-verified)
- [The test suite](#the-test-suite)
- [Bit-for-bit agreement with the reference](#bit-for-bit-agreement-with-the-reference)
- [End-to-end checks on real data](#end-to-end-checks-on-real-data)
- [Choosing the demo subhalo](#choosing-the-demo-subhalo)
- [Reproducing these checks](#reproducing-these-checks)

---

## What is being verified

The four models are implementations of published prescriptions, not new
physics:

| `hi_model` | Prescription | Projection |
|---|---|---|
| `GD14_vol` | Gnedin & Draine (2014), Appendix C.3 | volume |
| `GD14_map` | Gnedin & Draine (2014), §2.3 | projected |
| `K13_vol` | Krumholz, McKee & Tumlinson (2013), Appendix C.4 | volume |
| `K13_map` | Krumholz+ (2013), §2.3 | projected |

all as assembled by [Diemer+2018 (arXiv:1806.02341)](https://arxiv.org/abs/1806.02341).

So the question is not "is the physics right" — that is settled by the papers —
but "is this code the same code as the reference". Hence the emphasis on exact
agreement rather than on plausible-looking numbers.

The chain is: Diemer+2018's reference implementation
(`TNGHI_GD14_model/gd14_model.py`) → the `TNGHIcalc` package → `astrovodka.himodel`.
Each arrow was checked by re-running both sides on identical inputs.

---

## The test suite

```bash
python -m pytest tests/ -q      # 125 tests, ~14 s
```

It needs **no TNG data, no network and no optional dependencies**.

| Module | Tests | What it covers |
|---|---|---|
| `test_hi_models.py` | 52 | Keyword selection across the four methods, config resolution, per-model invariants (fraction ranges, mass conservation, full-length output), the four methods genuinely differing from one another, `uv_mode` handling, the off-centre FFT UV regression, and the underlying physics |
| `test_pipeline.py` | 13 | Cosmology mode, spectral-axis/velocity consistency, angular scale equals `r/D_A`, mass recovery at z = 0.02–1.0 |
| `test_selfcontained.py` | 19 | No `import martini` anywhere (checked by AST), no `martini.*` in `sys.modules` after importing astrovodka, every vendored and own module importable, `himodel` never importing the observation layer |
| `test_cosmology.py` | 14 | Distances, H(z), both velocity conventions, observed frequency |
| `test_docs_sync.py` | 11 | The generated HTML matches the Markdown; the English and Chinese trees stay parallel; relative links resolve |
| `test_data_products.py` | 16 | `spectral_unit` handling including the width-conversion regression, UV grid exposure, HDF5 round-trip |

Every defect fixed in 0.2.0 and 0.3.0 has a regression test.

Two of those tests are worth describing, because the obvious way to write them
would not have caught the bug:

- **The 0.2.0 spectral-width tests assert the numeric value**, not the unit. The
  original test only checked that the axis carried *a* frequency unit, which is
  why a cube retaining 0.007% of its flux passed for two releases.
- **The 0.3.0 FFT UV regression is placement-based, not a pinned number.** A
  broken interpolator returns the background wherever the source happens to
  sit, so the property worth asserting is "the field responds to where the
  source is" — pinning a value would have let the same failure back in under a
  different arrangement of the grid.

---

## Bit-for-bit agreement with the reference

The physical expressions were ported verbatim — operand order included, since
that is what the agreement rests on. The comparisons below were re-run **after**
the port:

| Comparison | Result |
|---|---|
| Every `himodel.physics` function against `TNGHIcalc`, on 4000 random cells spanning physical ranges (h = 0.6774, z = 0.5) | 12/12 `np.array_equal` |
| `compute_hi_mass` against `TNGHIcalc.cells_vol` / `cells_map`, all four methods, synthetic galaxy placed **off-centre** | bit-exact in masses, `f_mol`, `U_MW` and `row_index` |
| All four methods on **real TNG-100 snap 99 subhalos** 122 (4448 cells), 7 (120 cells) and 4 (86902 cells) | bit-exact |

### Why the unit conversions deserve their own paragraph

`rho_code_to_physical` deliberately goes through **astropy's unit chain** rather
than multiplying by a hard-coded constant. The two differ: astropy's kpc is
3.0856775814913673 × 10²¹ cm, while the module's own `kpc_to_cm` is 3.08568 ×
10²¹ cm — a relative difference of about 8 × 10⁻⁷. Substituting one for the
other is a silent change at the 10⁻⁶ level in **every** density, which
propagates into every molecular fraction.

This is precisely the class of error that bit-exact comparison catches and
spot-checking does not, which is why it is asserted.

### A trap removed by the port

The reference's `f_mol_*` functions take `z = 0.0` as a **positional** argument
with a default. Calling one without passing the redshift therefore does not fail
— it silently overrides the configured redshift with zero. In `himodel` the
redshift comes only from `HIModelConfig`, and there is no positional default to
fall into.

---

## End-to-end checks on real data

### All four models through the full pipeline

TNG-100 snapshot 99, subhalo 122, built through `TNGSource`:

| `hi_model` | `M_HI` [M☉] | `M_H₂` [M☉] | Cube mass recovery |
|---|---|---|---|
| `GD14_vol` | 2.4169 × 10⁹ | 1.0397 × 10⁹ | 100% |
| `K13_vol` | 3.2710 × 10⁹ | 1.8564 × 10⁸ | 100% |
| `GD14_map` | 2.1310 × 10⁹ | 1.2287 × 10⁹ | 100% |
| `K13_map` | 3.0689 × 10⁹ | 2.9075 × 10⁸ | 100% |

The projected models carry 88.3% of the gas mass inside their projection box,
so 474 cells get zero mass. Those are pruned automatically before cube
insertion, and the cube's `M_HI` still matches `source.m_HI_total` exactly —
which is the point: the zero-mass cells never contribute flux, so the recovery
stays at 100%.

### Distances, velocities and the spectral axis

Subhalo 122 under the default model:

| z | `vsys` | maps to | expected `f_rest/(1+z)` | angular scale | recovered mass |
|---|--------|---------|--------------------------|---------------|----------------|
| 0.02 | 5 878 km/s | 1392.55 MHz | 1392.55 MHz | exact | 100% |
| 0.10 | 27 254 km/s | 1291.28 MHz | 1291.28 MHz | exact | 100% |
| 0.50 | 99 931 km/s | 946.94 MHz | 946.94 MHz | exact | 100% |
| 1.00 | 149 896 km/s | 710.20 MHz | 710.20 MHz | exact | 100% |
| 2.00 | 199 862 km/s | 473.47 MHz | 473.47 MHz | exact | 100% |

Angular separations match `r / D_A` at every redshift, and mass recovery is
exact. Together these confirm the `D_L`/`D_A` split and the flux normalisation —
the two places where getting the cosmology wrong costs a factor of `(1+z)²` or
`(1+z)⁴`.

---

## Choosing the demo subhalo

The tutorial uses **subhalo 113** of TNG-100 snapshot 99. It was picked by
measuring candidates rather than by taking whichever had the most HI.

The requirement is HI *and* star formation. The adaptive FFT UV field only does
visible work in star-forming gas, so a gas-rich but quiescent subhalo
demonstrates nothing about it — the previous demonstration subhalo was
particularly bad in this respect, with the model assigning it about 16 M☉ of HI
across 120 gas cells and no star formation at all, which made its cubes
noise-dominated for a purely physical reason.

| | **113** (chosen) | 139 | 122 |
|---|---|---|---|
| Gas cells | 8 169 | 7 246 | 4 448 |
| `M_HI` [M☉] | 4.3 × 10⁹ | **4.5 × 10⁹** | 3.3 × 10⁹ |
| `M_H₂` [M☉] | **1.8 × 10⁹** | 3.0 × 10⁸ | 1.9 × 10⁸ |
| SFR [M☉/yr] | 0.98 | 0.17 | 0.19 |
| UV field peak / background | **27.6 / 0.0014 (≈ 19 000×)** | 2.2 | 6.0 |
| HI mass, `adaptive_fft` vs `global` | **+13.9%** | +5.7% | +4.3% |
| `f_gas_in_map` | **0.955** | 0.951 | 0.883 |

Subhalo 139 has 5% more HI, but its UV field peaks only about 1550× above the
background, so the model whose bug 0.3.0 fixed would look almost identical to
the broken one there. Subhalo 113 also keeps 95.5% of its gas inside the
projection box, which makes the `*_map` models directly comparable.

Measured end to end at z = 0.02: 100% mass recovery in both cubes, moment-0 peak
S/N 5.5 and total-flux S/N 16.9 for the 10 arcsec beam, 6.5 for the 3 arcmin
beam.

---

## Reproducing these checks

```bash
# 1. The unit tests -- no TNG data needed
python -m pytest tests/ -q

# 2. The tutorial, which asserts on its own results as it runs
#    (point it at your TNG data first)
export TNG_BASE_PATH=/path/to/TNG/TNG-100/
#    ... then run examples/{en,zh}/tutorial_TNG100-HI-mock-observations.ipynb
```

Several sections of the tutorial carry `assert` statements rather than printing
numbers for a reader to eyeball: the mass recovery must be 100%, the integrated
FITS flux must reproduce the model's own total to within 1%, and the velocity-
and frequency-axis cubes must agree in flux. If a change breaks the physics,
that notebook fails rather than quietly reporting a different number.
