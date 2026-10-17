# astrovodka.spectral_models — 21-cm Spectral Line Models

Models the 21-cm line profile per SPH particle.

These classes subclass MARTINI's `GaussianSpectrum` and `DiracDeltaSpectrum`; the only addition is the luminosity-distance flux scale described in the module docstring.

## GaussianSpectrum

### `__init__(self, sigma=7 km/s, ncpu=None, spec_dtype=np.float64)`

| Name | Type | Description |
|------|------|-------------|
| `sigma` | `Quantity[km/s] \| str` | Width. `'thermal'` = √(k<sub>B</sub>T/m<sub>p</sub>) per particle |
| `ncpu` | `int \| None` | Parallel CPUs (default: 1) |
| `spec_dtype` | `type` | Data type (default: float64) |

> Evaluated as `0.5 × [erf((b-v)/(σ√2)) - erf((a-v)/(σ√2))]` where (a,b) are channel edges.

---

## DiracDeltaSpectrum

### `__init__(self, ncpu=None, spec_dtype=np.float64)`

Delta-function model. All flux goes into the single channel containing the particle velocity.

---

> **How flux is computed:** Each particle's spectrum is integrated over each channel. The dimensionless result is multiplied by `mHI / D_L² / dv` and converted to Jy via: M<sub>HI</sub>/M<sub>sun</sub> = 2.36×10⁵ × (D<sub>L</sub>/Mpc)² × (S<sub>21</sub>/Jy km/s).
>
> **Custom spectra:** Inherit from `_BaseSpectrum`. Implement `half_width(source)` and `spectral_function(a, b, vmids, extra_data)`.
