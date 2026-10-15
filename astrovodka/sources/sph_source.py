"""
Provide the SPHSource class with cosmological distance support.

Modified from MARTINI's sph_source.py to support:
- Redshift-driven distances using astropy.cosmology (Planck15 default)
- Proper cosmological velocity computation instead of simple Hubble law
- Backward-compatible with the original distance + h parameter mode
"""

import numpy as np
import warnings
from scipy.spatial.transform import Rotation
from astropy.coordinates import (
    CartesianRepresentation,
    CartesianDifferential,
    SphericalRepresentation,
    SkyCoord,
    SpectralCoord,
    ICRS,
)
import astropy.units as U
from ._L_align import L_align
from ._extra_cartesian_transforms import (
    translate,
    translate_d,
    affine_transform,
    affine_transform_d,
)
from ..datacube import HIfreq
from ..cosmology_utils import (
    luminosity_distance,
    angular_diameter_distance,
    radio_velocity_from_redshift,
    hubble_parameter,
    get_cosmology,
)

# Extend CartesianRepresentation to allow more transformations
setattr(CartesianRepresentation, "translate", translate)
setattr(CartesianRepresentation, "affine_transform", affine_transform)
setattr(CartesianDifferential, "translate", translate_d)
setattr(CartesianDifferential, "affine_transform", affine_transform_d)

_origin = CartesianRepresentation(
    np.zeros((3, 1)) * U.kpc,
    differentials={"s": CartesianDifferential(np.zeros((3, 1)) * U.km / U.s)},
)

_COORDINATE_TRANSFORM_UNITS = U.Mpc
_VELOCITY_TRANSFORM_UNITS = U.km / U.s


class SPHSource(object):
    """
    Class abstracting HI emission sources consisting of SPH simulation particles.

    Extended from MARTINI's SPHSource with cosmological distance support. Can
    accept either a ``distance`` (Mpc) with Hubble parameter ``h`` (backward
    compatible with MARTINI), or a ``redshift`` with a ``cosmology`` instance
    for proper cosmological distance calculations at high redshift.

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

        # --- Validate distance/redshift mode ---
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
        self._cosmology = get_cosmology(cosmology) if self._use_cosmology else None
        self._redshift = redshift

        if self._use_cosmology:
            # Cosmological distance mode. Two different distances are needed:
            #   - D_L sets the flux scale (S ~ M_HI / D_L^2);
            #   - D_A sets every angular scale, because particles are placed on the sky
            #     at radius D_A, so a transverse separation dr subtends dr / D_A.
            # They differ by (1+z)^2, which is the whole point of this mode.
            self.distance = luminosity_distance(redshift, self._cosmology)
            self.luminosity_distance = self.distance
            self.angular_diameter_distance = angular_diameter_distance(
                redshift, self._cosmology
            )
            self._Hz = hubble_parameter(redshift, self._cosmology)
            # `vsys` is compared against the cube's spectral axis, which is built with
            # the *radio* Doppler convention (see `datacube`), so it must be expressed
            # in that same convention: v = c * z / (1+z) = c * (1 - f_rest/f_obs).
            # Anything else (e.g. H(z)*D_C, which exceeds c beyond z ~ 0.79) silently
            # misplaces the source in the band and eventually yields negative
            # frequencies.
            self.vhubble = radio_velocity_from_redshift(redshift).to(U.km / U.s)
            self.vsys = self.vhubble + vpeculiar
            # h is not used in cosmology mode, but set for interface compatibility
            self.h = self._cosmology.H0.to_value(U.km / U.s / U.Mpc) / 100.0
            # The distances are floored for z -> 0 (see `cosmology_utils`). Below that
            # floor the angular size of the source is not physically meaningful, so
            # warn rather than silently produce a cube containing almost nothing.
            if self.angular_diameter_distance < 0.1 * U.Mpc:
                warnings.warn(
                    f"Cosmology mode with redshift={redshift:g} gives an angular "
                    f"diameter distance of only {self.angular_diameter_distance:.3g}; "
                    "angular scales are not meaningful this close. To mock-observe a "
                    "z~0 galaxy at a chosen distance, pass `observed_redshift` "
                    "(TNGSource) or use the `distance` argument instead of `redshift`.",
                    UserWarning,
                )
        else:
            # Original MARTINI distance mode (backward compatible)
            self.distance = distance
            self.luminosity_distance = distance
            self.angular_diameter_distance = distance  # same in flat space
            self.h = h
            self.vhubble = (self.h * 100.0 * U.km * U.s**-1 * U.Mpc**-1) * self.distance
            self.vsys = self.vhubble + vpeculiar
            self._Hz = self.h * 100.0 * U.km / U.s / U.Mpc

        self.vpeculiar = vpeculiar

        # --- Coordinate handling (same as MARTINI) ---
        if coordinate_axis is None:
            if (xyz_g.shape[0] == 3) and (xyz_g.shape[1] != 3):
                coordinate_axis = 0
            elif (xyz_g.shape[0] != 3) and (xyz_g.shape[1] == 3):
                coordinate_axis = 1
            elif xyz_g.shape == (3, 3):
                raise RuntimeError(
                    "Cannot guess coordinate_axis with shape (3, 3), "
                    "provide explicitly."
                )
            else:
                raise RuntimeError(
                    "Incorrect coordinate shape (not (3, N) or (N, 3))."
                )

        if xyz_g.shape != vxyz_g.shape:
            raise ValueError(
                "xyz_g and vxyz_g must have matching shapes."
            )

        if coordinate_axis == 0:
            xyz_g = xyz_g.T
            vxyz_g = vxyz_g.T

        self.coordinate_frame = coordinate_frame
        self.T_g = T_g
        self.mHI_g = mHI_g
        self.input_mass = mHI_g.sum()
        self.coordinates_g = CartesianRepresentation(
            xyz_g,
            xyz_axis=1,
            differentials={"s": CartesianDifferential(vxyz_g, xyz_axis=1)},
        )
        self.hsm_g = hsm_g

        self.npart = self.coordinates_g.size

        self.ra = ra
        self.dec = dec
        self._coordinate_affine_transform = np.eye(4)
        self._velocity_affine_transform = np.eye(4)
        self.rotate(rotation=rotation, L_coords=L_coords)
        self.skycoords = None
        self.spectralcoords = None
        self.pixcoords = None

    def _init_skycoords(self, _reset=True):
        """
        Initialize the sky coordinates of the particles.

        The particles are placed at the *angular diameter distance*, which is the
        distance at which a transverse separation ``dr`` subtends the angle ``dr/D_A``,
        so that the angular size of the source is correct at any redshift. In the
        original MARTINI distance mode ``D_A`` and ``D_L`` are the same number, so this
        reduces to MARTINI's behaviour.

        Velocities are built from a uniform systemic term plus a per-particle Hubble
        gradient:

        - Distance mode: MARTINI's original scheme - a uniform ``vpeculiar`` boost plus
          the full linear Hubble flow ``100 h r``.
        - Cosmology mode: a uniform ``vsys`` boost (already in the radio convention the
          data cube uses) plus the *differential* Hubble flow within the source,
          ``H(z) dr``, which is the radio-convention velocity difference between two
          points separated by a physical distance ``dr``. This puts the source centroid
          exactly at ``vsys`` and removes any dependence of the gradient on the choice
          of ``D_A`` rather than ``D_L``.
        """
        distance_unit_vector = (
            SphericalRepresentation(self.ra, self.dec, 1)
            .represent_as(CartesianRepresentation)
            .xyz
        )
        distance_vector = distance_unit_vector * self.angular_diameter_distance
        # Uniform velocity added to every particle.
        if self._use_cosmology:
            uniform_velocity_vector = distance_unit_vector * self.vsys
        else:
            uniform_velocity_vector = distance_unit_vector * self.vpeculiar

        self.rotate(Rotation.from_euler("y", -self.dec.to_value(U.rad)))
        self.rotate(Rotation.from_euler("z", self.ra.to_value(U.rad)))
        self.translate(distance_vector)

        if self._use_cosmology:
            # Differential Hubble flow across the source: the radio-convention
            # velocity difference over a physical separation dr is H(z) * dr.
            # `coordinates_g.xyz` has shape (3, n_particles), the centroid (3,).
            xyz = self.coordinates_g.xyz
            centroid = distance_vector.reshape((-1,) + (1,) * (xyz.ndim - 1))
            vhubble_vectors = (self._Hz * (xyz - centroid).to(U.Mpc)).to(U.km / U.s)
        else:
            # Original linear Hubble flow
            vhubble_vectors = (
                self.h * 100.0 * U.km * U.s**-1 * U.Mpc**-1
            ) * self.coordinates_g.xyz

        self.boost(uniform_velocity_vector)
        self.coordinates_g.differentials["s"] = CartesianDifferential(
            self.coordinates_g.differentials["s"].d_xyz + vhubble_vectors
        )
        self.skycoords = SkyCoord(
            self.coordinates_g, frame=self.coordinate_frame, copy=True
        )
        origin_skycoord = SkyCoord(
            x=0 * U.kpc, y=0 * U.kpc, z=0 * U.kpc,
            v_x=0 * U.km / U.s, v_y=0 * U.km / U.s, v_z=0 * U.km / U.s,
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

    def _init_pixcoords(self, datacube, origin=0):
        """Initialize pixel coordinates of the particles."""
        assert self.skycoords is not None, (
            "Initialize source.skycoords before calling _init_pixcoords."
        )
        assert self.spectralcoords is not None, (
            "Initialize source.spectralcoords before calling _init_pixcoords."
        )
        skycoords_df_frame = self.skycoords.transform_to(datacube.coordinate_frame)
        spectralcoords_df_specsys = (
            self.spectralcoords.with_observer_stationary_relative_to(datacube.specsys)
        )
        self.pixcoords = (
            np.vstack(
                datacube.wcs.sub(3).wcs_world2pix(
                    skycoords_df_frame.ra.to(datacube.units[0]),
                    skycoords_df_frame.dec.to(datacube.units[1]),
                    spectralcoords_df_specsys.to(datacube.units[2]),
                    origin,
                )
            )
            * U.pix
        )

    def apply_mask(self, mask):
        """Remove particles from source arrays according to a mask."""
        if mask.size != self.npart:
            raise ValueError("Mask must have same length as particle arrays.")
        mask_sum = np.sum(mask)
        if mask_sum == 0:
            raise RuntimeError("No non-zero mHI source particles in target region.")
        self.npart = mask_sum
        if self.T_g is not None and not self.T_g.isscalar:
            self.T_g = self.T_g[mask]
        if not self.mHI_g.isscalar:
            self.mHI_g = self.mHI_g[mask]
        self.coordinates_g = self.coordinates_g[mask]
        if self.skycoords is not None:
            self.skycoords = self.skycoords[mask]
        if self.spectralcoords is not None:
            self.spectralcoords = self.spectralcoords[mask]
        if self.pixcoords is not None:
            self.pixcoords = self.pixcoords[:, mask]
        if self.hsm_g is not None and not self.hsm_g.isscalar:
            self.hsm_g = self.hsm_g[mask]

    def _append_to_coordinate_affine_transform(self, affine_transform):
        self._coordinate_affine_transform = np.dot(
            affine_transform, self._coordinate_affine_transform
        )

    def _append_to_velocity_affine_transform(self, affine_transform):
        self._velocity_affine_transform = np.dot(
            affine_transform, self._velocity_affine_transform
        )

    def rotate(self, rotation=None, *, L_coords=None):
        """Rotate the source."""
        args_given = (rotation is not None, L_coords is not None)
        if np.sum(args_given) == 0:
            return Rotation.identity()
        elif np.sum(args_given) > 1:
            raise ValueError("Multiple rotations in a single call not allowed.")

        if rotation is not None:
            do_rot = rotation.as_matrix()

        if L_coords is not None:
            do_rot = np.eye(3)
            do_rot = L_align(
                self.coordinates_g.get_xyz(),
                self.coordinates_g.differentials["s"].get_d_xyz(),
                self.mHI_g,
                frac=0.3,
                Laxis="x",
            ).dot(do_rot)
            do_rot = (
                Rotation.from_euler("x", L_coords.az_rot.to_value(U.rad))
                .as_matrix()
                .dot(do_rot)
            )
            do_rot = (
                Rotation.from_euler("y", L_coords.incl.to_value(U.rad))
                .as_matrix()
                .dot(do_rot)
            )
            do_rot = (
                Rotation.from_euler(
                    "x",
                    (
                        L_coords.pa - 90 * U.deg
                        if L_coords.incl >= 0
                        else L_coords.pa - 270 * U.deg
                    ).to_value(U.rad),
                )
                .as_matrix()
                .dot(do_rot)
            )

        affine_transform = np.eye(4)
        affine_transform[:3, :3] = do_rot
        self._append_to_coordinate_affine_transform(affine_transform)
        self._append_to_velocity_affine_transform(affine_transform)
        self.coordinates_g = self.coordinates_g.transform(do_rot)
        return Rotation.from_matrix(do_rot)

    def translate(self, translation_vector):
        """Translate the source."""
        self.coordinates_g = self.coordinates_g.translate(translation_vector)
        affine_transform = np.eye(4)
        affine_transform[:3, 3] = translation_vector.squeeze().to_value(
            _COORDINATE_TRANSFORM_UNITS
        )
        self._append_to_coordinate_affine_transform(affine_transform)

    def boost(self, boost_vector):
        """Apply an offset to the source velocity."""
        self.coordinates_g.differentials["s"] = self.coordinates_g.differentials[
            "s"
        ].translate(boost_vector)
        affine_transform = np.eye(4)
        affine_transform[:3, 3] = boost_vector.squeeze().to_value(
            _VELOCITY_TRANSFORM_UNITS
        )
        self._append_to_velocity_affine_transform(affine_transform)

    @property
    def current_rotation(self):
        """Current rotation matrix of the source."""
        return self._coordinate_affine_transform[:3, :3]

    def save_current_rotation(self, fname):
        """Output current rotation matrix to file."""
        np.savetxt(fname, self._coordinate_affine_transform[:3, :3])

    def save_current_affine_transformations(self, fname):
        """Output current affine transformation matrices to file."""
        np.save(
            fname,
            np.array(
                [self._coordinate_affine_transform, self._velocity_affine_transform]
            ),
        )

    def load_affine_transformations(self, fname):
        """Load affine transformation matrices from file."""
        coordinate_affine_transform, velocity_affine_transform = np.load(fname)
        if not np.allclose(
            coordinate_affine_transform[:3, :3], velocity_affine_transform[:3, :3]
        ):
            raise ValueError(
                "Rotation portion of position and velocity affine transforms "
                "do not match."
            )
        self.coordinates_g = self.coordinates_g.affine_transform(
            coordinate_affine_transform, _COORDINATE_TRANSFORM_UNITS
        )
        self.coordinates_g.differentials["s"] = self.coordinates_g.differentials[
            "s"
        ].affine_transform(velocity_affine_transform, _VELOCITY_TRANSFORM_UNITS)
        self._append_to_coordinate_affine_transform(coordinate_affine_transform)
        self._append_to_velocity_affine_transform(velocity_affine_transform)

    def preview(self, max_points=5000, fig=1, lim=None, vlim=None,
                point_scaling="auto", title="", save=None):
        """Produce a figure showing the source particle coordinates and velocities."""
        import matplotlib.pyplot as plt

        self._init_skycoords()
        assert self.skycoords is not None

        lim = (
            np.max(
                SkyCoord(self.ra, self.dec, frame=self.coordinate_frame)
                .separation(self.skycoords)
                .to_value(U.deg)
            )
            if lim is None
            else lim.to_value(U.deg)
        )
        vlim = (
            np.max(
                np.abs(
                    (self.skycoords.radial_velocity - self.vsys).to_value(U.km / U.s)
                )
            )
            if vlim is None
            else vlim.to_value(U.km / U.s)
        )
        wrapped_ra = self.skycoords.ra.wrap_at(self.ra + 180 * U.deg)
        cmask = np.logical_and.reduce(
            (
                np.abs((wrapped_ra - self.ra)).to_value(U.deg) < lim,
                np.abs((self.skycoords.dec - self.dec)).to_value(U.deg) < lim,
                np.abs((self.skycoords.radial_velocity - self.vsys)).to_value(
                    U.km / U.s
                )
                < vlim,
            )
        )
        nparts = cmask.sum()
        mask = np.arange(self.mHI_g.size)[cmask][:: max(nparts // max_points, 1)]
        hsm_factor = (
            1
            if (self.hsm_g is None or self.hsm_g.isscalar or mask.size <= 1)
            else (1 - (self.hsm_g[mask] / self.hsm_g[mask].max()) ** 0.1).to_value(
                U.dimensionless_unscaled
            )
            + (self.hsm_g[mask].min() == self.hsm_g[mask].max()) * 0.5
        )
        alpha = hsm_factor if point_scaling == "auto" else 1.0
        if self.hsm_g is None:
            size_scale = 1
        else:
            # `lim` is an angle on the sky, so converting it to a physical size on the
            # plot requires the angular diameter distance (not D_L).
            size_scale = (
                self.hsm_g.to_value(U.kpc)
                / (lim * U.deg * self.angular_diameter_distance).to_value(
                    U.kpc, U.dimensionless_angles()
                )
                if (self.hsm_g.isscalar or mask.size <= 1)
                else (
                    self.hsm_g[mask].to_value(U.kpc)
                    / (lim * U.deg * self.angular_diameter_distance).to_value(
                        U.kpc, U.dimensionless_angles()
                    )
                )
            )
        size = 300 * size_scale if point_scaling == "auto" else 10

        if newfig := not isinstance(fig, plt.Figure):
            figure = plt.figure(fig, figsize=(12, 4))
        else:
            figure = fig
        figure.suptitle(title)

        if newfig:
            sp1 = figure.add_subplot(1, 3, 1, aspect="equal")
            sp1.set_facecolor("#222222")
            sp2 = figure.add_subplot(1, 3, 2)
            sp2.set_facecolor("#222222")
            sp3 = figure.add_subplot(1, 3, 3)
            sp3.set_facecolor("#222222")
            sp1.set_xlabel(r"$\mathrm{RA}\,[\mathrm{deg}]$")
            sp1.set_ylabel(r"$\mathrm{Dec}\,[\mathrm{deg}]$")
            sp2.set_xlabel(r"$\mathrm{RA}\,[\mathrm{deg}]$")
            sp2.set_ylabel(r"$v_r\,[\mathrm{km}\,\mathrm{s}^{-1}]$")
            sp3.set_xlabel(r"$\mathrm{Dec}\,[\mathrm{deg}]$")
            sp3.set_ylabel(r"$v_r\,[\mathrm{km}\,\mathrm{s}^{-1}]$")
            vmin, vmax = -vlim, vlim
        else:
            sp1, sp2, sp3, _ = figure.axes
            cb = sp1.collections[0].colorbar
            if cb is not None:
                vmin, vmax = (
                    min(-vlim, cb.mappable.get_clim()[0]),
                    max(vlim, cb.mappable.get_clim()[1]),
                )

        scatter = sp1.scatter(
            wrapped_ra[mask].to_value(U.deg),
            self.skycoords.dec[mask].to_value(U.deg),
            c=(self.skycoords.radial_velocity[mask] - self.vsys).to_value(U.km / U.s),
            marker="o", cmap="coolwarm", edgecolor="None",
            s=size, vmin=vmin, vmax=vmax, alpha=alpha, zorder=0,
        )
        sp1.plot(
            [self.ra.to_value(U.deg)], [self.dec.to_value(U.deg)],
            marker="+", ls="None", mec="white", ms=6, zorder=1,
        )
        if newfig:
            sp1.set_xlim((self.ra.to_value(U.deg) + lim, self.ra.to_value(U.deg) - lim))
            sp1.set_ylim((self.dec.to_value(U.deg) - lim, self.dec.to_value(U.deg) + lim))
            cb = figure.colorbar(mappable=scatter, ax=sp1, orientation="horizontal")
            cb.set_label(r"$v_\mathrm{r} - v_\mathrm{sys}\,[\mathrm{km}\,\mathrm{s}^{-1}]$")

        sp2.scatter(
            wrapped_ra[mask].to_value(U.deg),
            self.skycoords.radial_velocity[mask].to_value(U.km / U.s),
            c="white", edgecolors="None", marker="o",
            s=size, alpha=alpha, zorder=0,
        )
        sp2.plot(
            [self.ra.to_value(U.deg)], [self.vsys.to_value(U.km / U.s)],
            marker="+", ls="None", mec="red", ms=6, zorder=1,
        )
        if newfig:
            sp2.set_xlim((self.ra.to_value(U.deg) + lim, self.ra.to_value(U.deg) - lim))
            sp2.set_ylim((
                self.vsys.to_value(U.km / U.s) - vlim,
                self.vsys.to_value(U.km / U.s) + vlim,
            ))

        sp3.scatter(
            self.skycoords.dec[mask].to_value(U.deg),
            self.skycoords.radial_velocity[mask].to_value(U.km / U.s),
            c="white", edgecolors="None", marker="o",
            s=size, alpha=alpha, zorder=0,
        )
        sp3.plot(
            [self.dec.to_value(U.deg)], [self.vsys.to_value(U.km / U.s)],
            marker="+", ls="None", mec="red", ms=6, zorder=1,
        )
        if newfig:
            sp3.set_xlim((self.dec.to_value(U.deg) - lim, self.dec.to_value(U.deg) + lim))
            sp3.set_ylim((
                self.vsys.to_value(U.km / U.s) - vlim,
                self.vsys.to_value(U.km / U.s) + vlim,
            ))

        figure.subplots_adjust(wspace=0.3)
        if save is not None:
            plt.savefig(save)
        return figure
