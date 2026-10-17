"""MARTINI's data cube with an explicit spectral-axis unit.

:class:`DataCube` is a thin subclass of :class:`martini.datacube.DataCube`.
Everything about the grid, WCS, padding and channel helpers comes from upstream;
astrovodka only adds the ``spectral_unit`` argument, which lets the spectral
axis be built in MHz while ``channel_width`` is still given in km/s (or the
other way round). The conversion between the two uses the *differential* of the
radio-convention Doppler relation, since a channel width is a difference rather
than an absolute velocity.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import astropy.units as U
from astropy import constants as const
from astropy.coordinates import ICRS, SpectralCoord
from martini.datacube import DataCube as _MartiniDataCube
from martini.datacube import HIfreq, _GlobalProfileDataCube

if TYPE_CHECKING:
    from astropy.coordinates import BaseRADecFrame

__all__ = ["DataCube", "HIfreq", "_GlobalProfileDataCube"]


def _velocity_width_to_frequency(
    width: U.Quantity[U.km / U.s], rest_freq: U.Quantity[U.Hz] = HIfreq
) -> U.Quantity[U.Hz]:
    """
    Convert a velocity channel *width* to a frequency channel width.

    :mod:`astropy`'s ``doppler_radio`` equivalency maps an *absolute* velocity to the
    *absolute* frequency at which that velocity is observed, so it cannot be applied to
    a width: doing so returns approximately ``rest_freq`` rather than a small offset.
    The correct width transform is the differential of that relation, which for the
    radio convention ``f = rest_freq * (1 - v / c)`` gives ``|df| = rest_freq * |dv| / c``.

    Parameters
    ----------
    width : ~astropy.units.Quantity[~astropy.units.km / ~astropy.units.s]
        Channel width in velocity units.

    rest_freq : ~astropy.units.Quantity[~astropy.units.Hz]
        Rest frequency of the line. (default: HI 21cm)

    Returns
    -------
    ~astropy.units.Quantity[~astropy.units.Hz]
        Channel width in frequency units.
    """
    return (rest_freq * np.abs(width) / const.c).to(U.Hz)


def _frequency_width_to_velocity(
    width: U.Quantity[U.Hz], rest_freq: U.Quantity[U.Hz] = HIfreq
) -> U.Quantity[U.km / U.s]:
    """
    Convert a frequency channel *width* to a velocity channel width.

    The inverse of :func:`_velocity_width_to_frequency`; see that function for why the
    ``doppler_radio`` equivalency must not be used directly on a width.

    Parameters
    ----------
    width : ~astropy.units.Quantity[~astropy.units.Hz]
        Channel width in frequency units.

    rest_freq : ~astropy.units.Quantity[~astropy.units.Hz]
        Rest frequency of the line. (default: HI 21cm)

    Returns
    -------
    ~astropy.units.Quantity[~astropy.units.km / ~astropy.units.s]
        Channel width in velocity units.
    """
    return (const.c * np.abs(width) / rest_freq).to(U.km / U.s)


class DataCube(_MartiniDataCube):
    """
    A MARTINI data cube with a selectable spectral-axis unit.

    Parameters
    ----------
    n_px_x : int
        Number of pixels along the x (RA) axis.

    n_px_y : int
        Number of pixels along the y (Dec) axis.

    n_channels : int
        Number of channels along the spectral axis.

    px_size : ~astropy.units.Quantity
        Angular scale of one pixel.

    channel_width : ~astropy.units.Quantity
        Step size along the spectral axis, in velocity or frequency units.

    spectral_centre : ~astropy.units.Quantity, optional
        Velocity (or frequency) of the centre along the spectral axis.

    ra : ~astropy.units.Quantity, optional
        Right ascension of the cube centre.

    dec : ~astropy.units.Quantity, optional
        Declination of the cube centre.

    stokes_axis : bool, optional
        Whether to add a Stokes axis to the cube.

    coordinate_frame : ~astropy.coordinates.builtin_frames.baseradec.BaseRADecFrame, optional
        The coordinate frame of the spatial axes.

    specsys : str, optional
        The spectral reference frame (standard of rest), selected from
        ``"gcrs"``, ``"icrs"``, ``"hcrs"``, ``"lsrk"``, ``"lsrd"``, ``"lsr"``.

    spectral_unit : str or None, optional
        Unit of the spectral axis: ``"MHz"`` for a frequency axis (WCS ``CTYPE3`` is
        ``FREQ``), ``"km/s"`` for a velocity axis (``CTYPE3`` is ``VRAD``), or ``None``
        (default) to adopt the units of ``channel_width``. This matches MARTINI's
        behaviour, so a velocity ``channel_width`` gives a velocity axis and a frequency
        ``channel_width`` gives a frequency axis. If the requested unit differs from the
        units of ``channel_width``, the width is converted using the differential of the
        radio-convention Doppler relation (see :func:`_velocity_width_to_frequency`).

    velocity_centre : ~astropy.units.Quantity, deprecated
        Deprecated, use spectral centre instead.
    """

    def __init__(
        self,
        *,
        n_px_x: int,
        n_px_y: int,
        n_channels: int,
        px_size: U.Quantity[U.arcsec],
        channel_width: U.Quantity[U.km * U.s**-1],
        spectral_centre: U.Quantity[U.km * U.s**-1] = 0.0 * U.km * U.s**-1,
        ra: U.Quantity[U.deg] = 0.0 * U.deg,
        dec: U.Quantity[U.deg] = 0.0 * U.deg,
        stokes_axis: bool = False,
        coordinate_frame: "BaseRADecFrame" = ICRS(),
        specsys: str = "icrs",
        spectral_unit: str | None = None,
        velocity_centre: None = None,  # deprecated
    ) -> None:
        if velocity_centre is not None:  # pragma: no cover
            warnings.warn(
                DeprecationWarning(
                    "velocity_centre is deprecated, use spectral_centre instead."
                )
            )
            spectral_centre = velocity_centre

        if spectral_unit not in ("MHz", "km/s", None):
            raise ValueError(
                "spectral_unit must be 'MHz', 'km/s', or None (auto-detect)."
            )

        channel_width_type = U.get_physical_type(channel_width)
        if channel_width_type not in ("frequency", "velocity"):
            raise ValueError("Channel width must have frequency or velocity units.")

        # Auto-detect reproduces MARTINI's behaviour: the spectral axis adopts
        # the units the caller supplied for channel_width.
        if spectral_unit is None:
            spectral_unit = "MHz" if channel_width_type == "frequency" else "km/s"

        if spectral_unit == "MHz":
            if channel_width_type == "velocity":
                channel_width = _velocity_width_to_frequency(channel_width)
            # spectral_centre is an absolute velocity/frequency, so the doppler
            # equivalency applies directly here.
            spectral_centre = SpectralCoord(
                spectral_centre,
                doppler_convention="radio",
                doppler_rest=HIfreq,
            ).to(U.Hz)
        else:  # "km/s"
            if channel_width_type == "frequency":
                channel_width = _frequency_width_to_velocity(channel_width)
            spectral_centre = SpectralCoord(
                spectral_centre,
                doppler_convention="radio",
                doppler_rest=HIfreq,
            ).to(U.km / U.s)

        super().__init__(
            n_px_x=n_px_x,
            n_px_y=n_px_y,
            n_channels=n_channels,
            px_size=px_size,
            channel_width=channel_width,
            spectral_centre=spectral_centre,
            ra=ra,
            dec=dec,
            stokes_axis=stokes_axis,
            coordinate_frame=coordinate_frame,
            specsys=specsys,
        )

        # Remembered so that the cube can be rebuilt consistently (see `copy`
        # and `_BaseMartini.reset`). The unit is also recoverable from the
        # stored channel_width, which is what `from_wcs` and `load_state` rely on.
        self._spectral_unit = spectral_unit
        return

    def copy(self):
        """
        Produce a copy of the :class:`~astrovodka.datacube.DataCube`.

        MARTINI's implementation does not carry over the Stokes axis, the
        coordinate frame, the spectral reference system or the spectral unit,
        so the copy is rebuilt here with all of them.

        Returns
        -------
        ~astrovodka.datacube.DataCube
            Copy of the :class:`~astrovodka.datacube.DataCube` object.
        """
        copy = type(self)(
            n_px_x=self.n_px_x,
            n_px_y=self.n_px_y,
            n_channels=self.n_channels,
            px_size=self.px_size,
            channel_width=self.channel_width,
            spectral_centre=self.spectral_centre,
            ra=self.ra,
            dec=self.dec,
            stokes_axis=self.stokes_axis,
            coordinate_frame=self.coordinate_frame,
            specsys=self.specsys,
            spectral_unit=self._spectral_unit,
        )
        copy.padx, copy.pady = self.padx, self.pady
        copy._wcs = self.wcs.copy()
        copy._freq_channel_mode = self._freq_channel_mode
        copy._channel_edges = self._channel_edges
        copy._channel_mids = self._channel_mids
        copy._array = self._array.copy()
        return copy
