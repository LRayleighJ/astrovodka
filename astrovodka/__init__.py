"""
astrovodka - Synthetic HI observations from cosmological SPH simulations.

A modular package for creating mock resolved HI line observations (data cubes)
of smoothed-particle hydrodynamics simulations of galaxies, with:

- Proper cosmological distance handling via astropy.cosmology (default Planck15)
- K13 (Krumholz 2013) volumetric HI/H2 model for accurate gas masses
- Support for IllustrisTNG simulations with local data access
"""

from .martini import Martini, GlobalProfile
from .datacube import DataCube
from .data_products import SubhaloDataProducts
from .L_coords import L_coords
from .__version__ import __version__

__all__ = [
    "Martini",
    "GlobalProfile",
    "DataCube",
    "SubhaloDataProducts",
    "L_coords",
    "__version__",
]
