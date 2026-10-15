"""
astrovodka source classes for IllustrisTNG synthetic HI observations.

Provides cosmology-aware source classes. The HI/H2 masses that feed them come
from `astrovodka.himodel`, which offers four selectable models; `TNGSource`
exposes them through its ``hi_model`` keyword.

Other simulation-specific sources from upstream MARTINI (EAGLE, Simba, FIRE,
Magneticum, Colibre, SWIFT, simobj) are not included: astrovodka is built
around TNG data on local disk, and vendoring sources it cannot test would mean
shipping code nobody here can check.
"""

from .sph_source import SPHSource
from .tng_source import TNGSource

__all__ = [
    "SPHSource",
    "TNGSource",
]
