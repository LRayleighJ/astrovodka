"""Re-export MARTINI's noise models.

The implementation lives in :mod:`martini.noise`; astrovodka imports it rather
than carrying a copy. This module exists only to keep
``from astrovodka.noise import GaussianNoise`` working for existing code.
"""

from martini.noise import GaussianNoise, _BaseNoise

__all__ = ["_BaseNoise", "GaussianNoise"]
