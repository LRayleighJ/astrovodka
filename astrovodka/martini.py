"""The mock-observation driver, delegating to :mod:`martini.martini`.

The cube-making pipeline (particle pruning, pixel spectra, beam convolution,
noise, FITS/HDF5 output) is MARTINI's and is inherited unchanged. astrovodka
overrides only the places where a cosmological distance or an explicitly chosen
spectral unit needs to be carried through:

- :meth:`~astrovodka.martini.Martini.reset` rebuilds the cube with its Stokes
  axis, coordinate frame, spectral reference system and spectral unit intact
  (MARTINI's reset drops them);
- :meth:`~astrovodka.martini.Martini.preview` and
  :meth:`~astrovodka.martini.GlobalProfile.insert_source_in_spectrum` present a
  frequency axis to MARTINI in velocity units while those diagnostics run, since
  MARTINI converts the channel width to km/s unconditionally.
"""

from __future__ import annotations

from contextlib import contextmanager

import astropy.units as U
from martini.martini import GlobalProfile as _MartiniGlobalProfile
from martini.martini import Martini as _MartiniMartini
from martini.martini import _BaseMartini

from .datacube import DataCube, _frequency_width_to_velocity

__all__ = ["Martini", "GlobalProfile", "_BaseMartini"]


@contextmanager
def _velocity_axis(datacube):
    """Present ``datacube``'s spectral axis in velocity units for MARTINI.

    MARTINI's ``preview`` and global-profile diagnostics convert
    ``datacube.channel_width`` and ``datacube.spectral_centre`` to km/s without
    a Doppler equivalency, which raises for a frequency-axis cube. Inside this
    context those two attributes are the equivalent velocity width and centre,
    and they are restored on exit (including on error).

    Parameters
    ----------
    datacube : ~astrovodka.datacube.DataCube
        The cube whose spectral axis should be presented in velocity units.
    """
    channel_width = datacube.channel_width
    if U.get_physical_type(channel_width) != "frequency":
        yield
        return
    original = channel_width, datacube.spectral_centre
    datacube.channel_width = _frequency_width_to_velocity(channel_width)
    datacube.spectral_centre = datacube.spectral_centre.to(U.km / U.s)
    try:
        yield
    finally:
        datacube.channel_width, datacube.spectral_centre = original


class _AstrovodkaMartini:
    """Overrides shared by :class:`Martini` and :class:`GlobalProfile`."""

    def reset(self) -> None:
        """Re-initialize the data cube with zero values, keeping its setup."""
        datacube = self._datacube
        self._datacube = DataCube(
            n_px_x=datacube.n_px_x,
            n_px_y=datacube.n_px_y,
            n_channels=datacube.n_channels,
            px_size=datacube.px_size,
            channel_width=datacube.channel_width,
            spectral_centre=datacube.spectral_centre,
            ra=datacube.ra,
            dec=datacube.dec,
            stokes_axis=datacube.stokes_axis,
            coordinate_frame=datacube.coordinate_frame,
            specsys=datacube.specsys,
            spectral_unit=getattr(datacube, "_spectral_unit", None),
        )
        if getattr(self, "beam", None) is not None:
            self._datacube.add_pad(self.beam.needs_pad())
        if hasattr(self, "_spectrum"):
            del self._spectrum
        return

    def preview(self, *args, **kwargs):
        """Produce preview figures; delegates to MARTINI."""
        with _velocity_axis(self._datacube):
            return super().preview(*args, **kwargs)


class Martini(_AstrovodkaMartini, _MartiniMartini):
    """MARTINI's :class:`~martini.martini.Martini` driver for astrovodka cubes.

    Parameters are the same as MARTINI's: ``source``, ``datacube``, ``beam``,
    ``noise``, ``sph_kernel``, ``spectral_model`` and ``quiet``.
    """


class GlobalProfile(_AstrovodkaMartini, _MartiniGlobalProfile):
    """MARTINI's spatially-integrated spectrum driver.

    Parameters are the same as MARTINI's: ``source``, ``spectral_model``,
    ``n_channels``, ``channel_width``, ``spectral_centre`` and ``quiet``.
    """

    def insert_source_in_spectrum(self) -> None:
        """Insert the source; delegates to MARTINI."""
        with _velocity_axis(self._datacube):
            return super().insert_source_in_spectrum()
