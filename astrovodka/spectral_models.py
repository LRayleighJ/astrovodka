"""Spectral models, with the cosmological flux scaling astrovodka needs.

MARTINI's spectral models normalise the flux of each particle by the square of
the distance at which the particle was placed on the sky. For a source in
astrovodka's redshift mode that is the *angular diameter* distance ``D_A``,
whereas the observed flux of an HI mass ``M`` goes as ``M / D_L**2`` with the
*luminosity* distance ``D_L``. The two differ by ``(1+z)**2``, so without a
correction the recovered HI mass would be suppressed by ``(1+z)**4``.

The classes below are thin subclasses of :mod:`martini.spectral_models`: they
delegate the entire calculation upstream and rescale the result afterwards.
For sources that do not carry separate distances (plain MARTINI sources, or
astrovodka sources in ``distance=`` mode) the correction is exactly ``1``.
"""

from __future__ import annotations

import astropy.units as U
from martini.spectral_models import DiracDeltaSpectrum as _MARTINIDiracDeltaSpectrum
from martini.spectral_models import GaussianSpectrum as _MARTINIGaussianSpectrum
from martini.spectral_models import _BaseSpectrum

__all__ = [
    "_BaseSpectrum",
    "DiracDeltaSpectrum",
    "GaussianSpectrum",
]


def _flux_scale_correction(source):
    """Return ``(D_A / D_L)**2`` for ``source``, or ``1`` if it has one distance.

    Parameters
    ----------
    source : ~martini.sources.sph_source.SPHSource
        The source the spectra are being computed for. astrovodka sources in
        redshift mode expose both distance attributes; everything else is
        treated as having a single distance, for which no correction is needed.

    Returns
    -------
    float
        Multiplicative correction to apply to the spectra.
    """
    distance_luminosity = getattr(source, "luminosity_distance", None)
    if distance_luminosity is None:
        return 1.0
    distance_angular = getattr(source, "angular_diameter_distance", None)
    if distance_angular is None:
        if source.skycoords is None:
            return 1.0
        distance_angular = source.skycoords.distance
    ratio = (distance_angular / distance_luminosity).to_value(U.dimensionless_unscaled)
    return float(ratio) ** 2


class _CosmologicalFluxScale:
    """Rescale spectra from the placement distance to the luminosity distance."""

    def init_spectra(self, source, datacube):
        """
        Pre-compute the spectrum of each particle.

        The calculation is MARTINI's; the result is then corrected for the
        ``D_A``/``D_L`` split described in the module docstring.

        Parameters
        ----------
        source : ~martini.sources.sph_source.SPHSource
            Source object containing arrays of particle properties.

        datacube : ~martini.datacube.DataCube
            Data cube defining the observational parameters.
        """
        super().init_spectra(source, datacube)
        correction = _flux_scale_correction(source)
        if correction != 1.0:
            self.spectra = self.spectra * correction
        return


class GaussianSpectrum(_CosmologicalFluxScale, _MARTINIGaussianSpectrum):
    """MARTINI's Gaussian line model with the cosmological flux scaling."""


class DiracDeltaSpectrum(_CosmologicalFluxScale, _MARTINIDiracDeltaSpectrum):
    """MARTINI's Dirac-delta model with the cosmological flux scaling."""
