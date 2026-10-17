"""Re-export MARTINI's SPH smoothing kernels.

The implementation lives in :mod:`martini.sph_kernels`; astrovodka imports it
rather than carrying a copy. This module exists so that the names astrovodka's
documentation and examples use (``CubicSplineKernel`` and friends) resolve
inside this package.
"""

from martini.sph_kernels import (
    AdaptiveKernel,
    CubicSplineKernel,
    DiracDeltaKernel,
    GaussianKernel,
    QuarticSplineKernel,
    WendlandC2Kernel,
    WendlandC6Kernel,
    _AdaptiveKernel,
    _BaseSPHKernel,
    _CubicSplineKernel,
    _GaussianKernel,
    _QuarticSplineKernel,
    _WendlandC2Kernel,
    _WendlandC6Kernel,
    find_fwhm,
)

__all__ = [
    "AdaptiveKernel",
    "CubicSplineKernel",
    "DiracDeltaKernel",
    "GaussianKernel",
    "QuarticSplineKernel",
    "WendlandC2Kernel",
    "WendlandC6Kernel",
    "_AdaptiveKernel",
    "_BaseSPHKernel",
    "_CubicSplineKernel",
    "_GaussianKernel",
    "_QuarticSplineKernel",
    "_WendlandC2Kernel",
    "_WendlandC6Kernel",
    "find_fwhm",
]
