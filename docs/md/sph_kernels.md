# astrovodka.sph_kernels — SPH Smoothing Kernels

Smooth SPH particles onto the pixel grid. All public kernels include **adaptive fallback**: when the pixel-to-smoothing-length ratio is too large, they automatically fall back to `DiracDeltaKernel` then `GaussianKernel(truncate=6)`.

## Public Kernel Classes

| Class | min_valid_size | Typical Simulation |
|-------|---------------|-------------------|
| `CubicSplineKernel()` | 1.16 px | Gadget/Gadget2 (classic M4) |
| `WendlandC2Kernel()` | 1.51 px | EAGLE |
| `WendlandC6Kernel()` | 1.29 px | Magneticum |
| `QuarticSplineKernel()` | 1.24 px | Colibre (SPHENIX) |
| `GaussianKernel(truncate=3.0)` | 2.34 px | General purpose |
| `DiracDeltaKernel()` | ≤0.5 px | Point-like particles |

> **min_valid_size**: Minimum smoothing length (pixels) for 1% kernel integral accuracy. Below this, the kernel falls back.
>
> Use **underscore-prefixed** versions (e.g., `_CubicSplineKernel`) for raw kernels without fallback.
>
> Set `Martini.insert_source_in_cube(skip_validation=True)` to bypass accuracy checks.

### `find_fwhm(f: Callable) -> float`

Solves for FWHM of a kernel function via `scipy.optimize.fsolve`. Assumes max at f(0) and symmetry.
