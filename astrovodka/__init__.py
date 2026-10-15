"""
astrovodka - Synthetic HI observations from cosmological SPH simulations.

A modular package for creating mock resolved HI line observations (data cubes)
of smoothed-particle hydrodynamics simulations of galaxies, with:

- Proper cosmological distance handling via astropy.cosmology (default Planck15)
- Four selectable HI/H2 models from Diemer+2018 (Gnedin & Draine 2014 and
  Krumholz+ 2013, each in a volume and a projected variant), replacing
  MARTINI's pressure-based partition
- Support for IllustrisTNG simulations with local data access

The mock-observation machinery (data cubes, beams, spectral models, SPH
kernels, noise) is MARTINI's, vendored into this package rather than imported
from it, so astrovodka has no runtime dependency on `astromartini`.
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
