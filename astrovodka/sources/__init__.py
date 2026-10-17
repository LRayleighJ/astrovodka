"""
astrovodka source classes for IllustrisTNG synthetic HI observations.

:class:`SPHSource` extends MARTINI's source class with cosmological distances;
:class:`TNGSource` loads TNG data from local disk and computes HI/H2 masses
with one of the four ``himodel`` models.

The other simulation-specific sources shipped by upstream MARTINI (EAGLE,
Simba, FIRE, Magneticum, Colibre, SWIFT, simobj, and ``CombinedSource``) are
still available from :mod:`martini.sources`; astrovodka simply does not
re-export them, because it is built around TNG data on local disk.
"""

from .sph_source import SPHSource
from .tng_source import TNGSource

__all__ = [
    "SPHSource",
    "TNGSource",
]
