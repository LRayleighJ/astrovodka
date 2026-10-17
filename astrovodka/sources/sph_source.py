"""Cosmology-aware SPH source built on :class:`martini.sources.SPHSource`.

astrovodka does not reimplement MARTINI's source machinery. Particle storage,
rotations, sky coordinates, masking and previews all come from
:class:`martini.sources.SPHSource`; this class only adds the redshift-driven
distance mode that is astrovodka's reason to exist:

- in ``redshift=`` mode the flux scale uses the **luminosity** distance ``D_L``
  while particles are placed on the sky at the **angular diameter** distance
  ``D_A`` (they differ by ``(1+z)**2``);
- the systemic velocity is the radio-convention recession velocity
  ``c z / (1+z)``, so the source lands on the cube's spectral axis;
- the velocity gradient across the source is ``H(z) dr``, a differential
  Hubble flow, rather than the linear ``100 h r`` law.

The distance-mode path (``distance=...`` with ``h=...``) is passed straight to
MARTINI and is therefore bit-for-bit upstream behaviour.
"""

from __future__ import annotations

import inspect
import warnings

import astropy.units as U
from astropy.coordinates import (
    CartesianDifferential,
    CartesianRepresentation,
    ICRS,
    SkyCoord,
    SpectralCoord,
    SphericalRepresentation,
)
from scipy.spatial.transform import Rotation

from martini.datacube import HIfreq
from martini.sources import SPHSource as _MartiniSPHSource

from ..cosmology_utils import (
    angular_diameter_distance,
    get_cosmology,
    hubble_parameter,
    luminosity_distance,
    radio_velocity_from_redshift,
)

# MARTINI changed its rotation API in 2.1.16: `rotate(rotation=Rotation)` replaced
# `rotate(axis_angle=..., rotmat=..., L_coords=(...))`. astrovodka supports every
# released version, so the adapter below dispatches on the installed signature
# rather than on a version string.
_UPSTREAM_ROTATE_TAKES_ROTATION = (
    "rotation" in inspect.signature(_MartiniSPHSource.rotate).parameters
)


def _legacy_l_coords(L_coords):
    """Convert astrovodka's ``L_coords`` to the pre-2.1.16 tuple form."""
    if L_coords is None:
        return None
    if hasattr(L_coords, "incl"):
        return (L_coords.incl, L_coords.az_rot, L_coords.pa)
    return tuple(L_coords)


class SPHSource(_MartiniSPHSource):
    """
    Class abstracting HI emission sources consisting of SPH simulation particles.

    Extended from MARTINI's :class:`~martini.sources.SPHSource` with cosmological
    distance support. Can accept either a ``distance`` (Mpc) with Hubble
    parameter ``h`` (backward compatible with MARTINI), or a ``redshift`` with a
    ``cosmology`` instance for proper cosmological distance calculations at high
    redshift.

    Parameters
    ----------
    distance : ~astropy.units.Quantity, optional
        Source distance, also used to set the velocity offset via Hubble's law.
        If provided, ``h`` must also be given. Mutually exclusive with ``redshift``.

    redshift : float, optional
        Cosmological redshift of the source. If provided, distances and velocities
        are computed using ``cosmology``. Mutually exclusive with ``distance``.

    cosmology : astropy.cosmology.Cosmology, optional
        Cosmology instance to use. Defaults to Planck15. Only used when
        ``redshift`` is provided.

    vpeculiar : ~astropy.units.Quantity, optional
        Source peculiar velocity along the direction to the source centre.

    rotation : ~scipy.spatial.transform.Rotation, optional
        A rotation to apply to the source particles.

    L_coords : ~astrovodka.L_coords.L_coords, optional
        A named tuple specifying 3 angles for angular momentum alignment.

    ra : ~astropy.units.Quantity, optional
        Right ascension for the source centroid.

    dec : ~astropy.units.Quantity, optional
        Declination for the source centroid.

    h : float, optional
        Dimensionless Hubble constant. Only used in ``distance`` mode.

    T_g : ~astropy.units.Quantity, optional
        Particle temperature.

    mHI_g : ~astropy.units.Quantity
        Particle HI mass.

    xyz_g : ~astropy.units.Quantity
        Particle position offset from source centroid. 'x' is line-of-sight,
        'y-z' is the sky plane.

    vxyz_g : ~astropy.units.Quantity
        Particle velocity offset from source centroid.

    hsm_g : ~astropy.units.Quantity, optional
        Particle SPH smoothing lengths, defined as FWHM of the smoothing kernel.

    coordinate_axis : int, optional
        Rank of axis corresponding to position/velocity of a single particle.

    coordinate_frame : ~astropy.coordinates.builtin_frames.baseradec.BaseRADecFrame
        The coordinate frame for converting particle coordinates to RA and Dec.
    """

    def __init__(
        self,
        *,
        distance=None,
        redshift=None,
        cosmology=None,
        vpeculiar=0.0 * U.km / U.s,
        rotation=None,
        L_coords=None,
        ra=0.0 * U.deg,
        dec=0.0 * U.deg,
        h=0.7,
        T_g=None,
        mHI_g=None,
        xyz_g=None,
        vxyz_g=None,
        hsm_g=None,
        coordinate_axis=None,
        coordinate_frame=ICRS(),
    ):
        if isinstance(rotation, dict):
            raise ValueError(
                "The method to specify rotations in astrovodka has been updated. "
                "Use scipy.spatial.transform.Rotation or L_coords instead."
            )

        if distance is not None and redshift is not None:
            raise ValueError(
                "Specify either 'distance' (with 'h') or 'redshift' "
                "(with optional 'cosmology'), not both."
            )
        if distance is None and redshift is None:
            raise ValueError(
                "Either 'distance' (with 'h') or 'redshift' must be provided."
            )

        self._use_cosmology = redshift is not None
        self._redshift = redshift
        self._cosmology = get_cosmology(cosmology) if self._use_cosmology else None

        if self._use_cosmology:
            # Two different distances are needed: D_L sets the flux scale
            # (S ~ M_HI / D_L^2), D_A sets every angular scale, because particles
            # are placed on the sky at radius D_A. They differ by (1+z)^2.
            self.luminosity_distance = luminosity_distance(redshift, self._cosmology)
            self.angular_diameter_distance = angular_diameter_distance(
                redshift, self._cosmology
            )
            self._Hz = hubble_parameter(redshift, self._cosmology)
            # `vsys` is compared against the cube's spectral axis, which is built
            # with the radio Doppler convention: v = c z / (1+z).
            self._cosmology_vhubble = radio_velocity_from_redshift(redshift).to(
                U.km / U.s
            )
            rotation_distance = self.luminosity_distance
        else:
            self.luminosity_distance = distance
            self.angular_diameter_distance = distance
            self._Hz = h * 100.0 * U.km / U.s / U.Mpc
            self._cosmology_vhubble = None
            rotation_distance = distance

        # MARTINI places the source and computes its diagnostics with a single
        # `distance`. The luminosity distance is the right one for both; angular
        # scales are corrected in `_init_skycoords`, which would otherwise use it.
        super().__init__(
            distance=rotation_distance,
            vpeculiar=vpeculiar,
            ra=ra,
            dec=dec,
            h=h,
            T_g=T_g,
            mHI_g=mHI_g,
            xyz_g=xyz_g,
            vxyz_g=vxyz_g,
            hsm_g=hsm_g,
            coordinate_axis=coordinate_axis,
            coordinate_frame=coordinate_frame,
        )

        # MARTINI's own __init__ applies its default (identity) rotation; apply
        # the caller's afterwards so that both generations of its API work.
        if rotation is not None or L_coords is not None:
            self.rotate(rotation=rotation, L_coords=L_coords)

        if self._use_cosmology:
            self.h = self._cosmology.H0.to_value(U.km / U.s / U.Mpc) / 100.0
            self.vhubble = self._cosmology_vhubble
            self.vsys = self.vhubble + vpeculiar
            # The distances are floored for z -> 0 (see `cosmology_utils`). Below
            # that floor the angular size is not physically meaningful, so warn
            # rather than silently produce a cube containing almost nothing.
            if self.angular_diameter_distance < 0.1 * U.Mpc:
                warnings.warn(
                    f"Cosmology mode with redshift={redshift:g} gives an angular "
                    f"diameter distance of only {self.angular_diameter_distance:.3g}; "
                    "angular scales are not meaningful this close. To mock-observe a "
                    "z~0 galaxy at a chosen distance, pass `observed_redshift` "
                    "(TNGSource) or use the `distance` argument instead of `redshift`.",
                    UserWarning,
                )
        return

    def rotate(self, rotation=None, *, L_coords=None, axis_angle=None, rotmat=None):
        """
        Rotate the source.

        Accepts both MARTINI rotation interfaces: astrovodka's
        ``rotation=<scipy Rotation>`` (MARTINI 2.1.16 onward) and the older
        ``axis_angle=``/``rotmat=`` form, whichever version is installed.

        Parameters
        ----------
        rotation : ~scipy.spatial.transform.Rotation, optional
            The rotation to apply.

        L_coords : ~astrovodka.L_coords.L_coords, optional
            Angular momentum alignment angles.

        axis_angle : tuple, optional
            Legacy form: ``(axis, angle)``.

        rotmat : array-like, optional
            Legacy form: a 3x3 rotation matrix.
        """
        if rotation is None:
            if axis_angle is not None:
                rotation = Rotation.from_euler(
                    axis_angle[0], axis_angle[1].to_value(U.rad)
                )
            elif rotmat is not None:
                rotation = Rotation.from_matrix(rotmat)

        if _UPSTREAM_ROTATE_TAKES_ROTATION:
            return super().rotate(rotation=rotation, L_coords=L_coords)
        return super().rotate(
            rotmat=None if rotation is None else rotation.as_matrix(),
            L_coords=_legacy_l_coords(L_coords),
        )

    def _init_skycoords(self, _reset=True):
        """
        Initialize the sky coordinates of the particles.

        In ``distance`` mode this is MARTINI's implementation. In ``redshift``
        mode the particles are placed at the *angular diameter distance*, so
        that a transverse separation ``dr`` subtends the angle ``dr / D_A``, and
        the velocities are built from a uniform systemic term plus the
        differential Hubble flow ``H(z) dr`` within the source.

        Parameters
        ----------
        _reset : bool
            If ``True``, return particles to their original positions. Setting to
            ``False`` is only intended for testing.
        """
        if not self._use_cosmology:
            return super()._init_skycoords(_reset=_reset)

        distance_unit_vector = (
            SphericalRepresentation(self.ra, self.dec, 1)
            .represent_as(CartesianRepresentation)
            .xyz
        )
        distance_vector = distance_unit_vector * self.angular_diameter_distance
        # The systemic velocity already includes the Hubble term, so every
        # particle receives the same uniform boost.
        uniform_velocity_vector = distance_unit_vector * self.vsys

        self.rotate(Rotation.from_euler("y", -self.dec.to_value(U.rad)))
        self.rotate(Rotation.from_euler("z", self.ra.to_value(U.rad)))
        self.translate(distance_vector)

        # Differential Hubble flow across the source: the radio-convention
        # velocity difference over a physical separation dr is H(z) * dr.
        # `coordinates_g.xyz` has shape (3, n_particles), the centroid (3,).
        xyz = self.coordinates_g.xyz
        centroid = distance_vector.reshape((-1,) + (1,) * (xyz.ndim - 1))
        vhubble_vectors = (self._Hz * (xyz - centroid).to(U.Mpc)).to(U.km / U.s)

        self.boost(uniform_velocity_vector)
        self.coordinates_g.differentials["s"] = CartesianDifferential(
            self.coordinates_g.differentials["s"].d_xyz + vhubble_vectors
        )
        self.skycoords = SkyCoord(
            self.coordinates_g, frame=self.coordinate_frame, copy=True
        )
        origin_skycoord = SkyCoord(
            x=0 * U.kpc,
            y=0 * U.kpc,
            z=0 * U.kpc,
            v_x=0 * U.km / U.s,
            v_y=0 * U.km / U.s,
            v_z=0 * U.km / U.s,
            representation_type="cartesian",
            differential_type="cartesian",
            frame=self.coordinate_frame,
        )
        self.spectralcoords = SpectralCoord(
            self.skycoords.radial_velocity,
            doppler_convention="radio",
            doppler_rest=HIfreq,
            target=self.skycoords,
            observer=origin_skycoord,
        )
        if _reset:
            self.coordinates_g.differentials["s"] = CartesianDifferential(
                self.coordinates_g.differentials["s"].d_xyz - vhubble_vectors
            )
            self.boost(-uniform_velocity_vector)
            self.translate(-distance_vector)
            self.rotate(Rotation.from_euler("z", -self.ra.to_value(U.rad)))
            self.rotate(Rotation.from_euler("y", self.dec.to_value(U.rad)))
        return
