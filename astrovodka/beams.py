"""Re-export MARTINI's beam models.

The implementation lives in :mod:`martini.beams`; astrovodka imports it so that
the beam model used here is always the one maintained upstream. This module
exists only to keep ``from astrovodka.beams import GaussianBeam`` working for
existing code.
"""

from martini.beams import GaussianBeam, _BaseBeam

__all__ = ["_BaseBeam", "GaussianBeam"]
