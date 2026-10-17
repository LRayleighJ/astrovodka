"""
3D gridded data products for TNG subhalos.

Provides SubhaloDataProducts — a container for 3D arrays of UV field,
DM/HI/H2 mass density (comoving real-space and redshift-space), with
HDF5 save/load support.
"""

import numpy as np
import astropy.units as U
from astropy.cosmology import Planck15

from .himodel import (
    DEFAULT_METHOD,
    HIModelConfig,
    calc_UMW_global_Faucher11,
    compute_hi_mass,
    is_known,
    resolve,
)
from .himodel.uvfield import compute_galaxy_fft_uv
from martini.sources._illustris_tools import (
    getSnapOffsets,
    loadHeader,
    loadSingle,
    loadSubset,
)


class SubhaloDataProducts:
    """Container for 3D gridded data products of a TNG subhalo.

    Contains 7 grids of shape (ngrid, ngrid, ngrid):
    - uv_field: real-space UV radiation field [Draine units]
    - dm_real, dm_rsd: comoving DM mass density [Msun / (ckpc/h)^3]
      (the real-space grids sit on `grid_edges`, the redshift-space ones on the
      longer `grid_edges_rsd`, so the two are not cell-by-cell comparable)
    - hi_real, hi_rsd: comoving HI mass density [Msun / (ckpc/h)^3]
    - h2_real, h2_rsd: comoving H2 mass density [Msun / (ckpc/h)^3]

    Real-space uses comoving coordinates (ckpc/h).
    Redshift-space (RSD) shifts the LOS coordinate by v_los/(a*H(a)), with the
    velocities measured in the object's own frame (the DM centre of mass) and on
    its own grid box, `grid_edges_rsd` -- see `_rsd_los_halfwidth` for why that
    box is longer than the real-space one.
    """

    def __init__(self, ngrid=128):
        self.ngrid = ngrid
        self.uv_field = None
        self.dm_real = None
        self.dm_rsd = None
        self.hi_real = None
        self.hi_rsd = None
        self.h2_real = None
        self.h2_rsd = None
        self.grid_edges = None  # 3-tuple of 1D edge arrays [ckpc/h], real space
        self.grid_edges_rsd = None  # the same for the redshift-space box
        self.metadata = {}

    @classmethod
    def from_subhalo(
        cls,
        basePath,
        snapNum,
        subID,
        ngrid=128,
        hi_model=DEFAULT_METHOD,
        hi_config=None,
        uv_mode="adaptive_fft",
        load_subhalo_only=True,
        cosmology=None,
    ):
        """Load data and compute all 3D gridded products for a subhalo.

        Parameters
        ----------
        basePath : str — TNG data root directory
        snapNum : int — Snapshot number
        subID : int — Subhalo ID
        ngrid : int — Grid resolution per dimension (default 128)
        hi_model : str — one of 'GD14_vol', 'GD14_map', 'K13_vol', 'K13_map'
            (case-insensitive; 'GD14'/'K13' mean the volume variant).
        hi_config : HIModelConfig or dict, optional — HI model settings
        uv_mode : str — 'adaptive_fft' (default) or 'global'; volume models only
        load_subhalo_only : bool — Load subhalo-only (True) or full FOF group
        cosmology : astropy.cosmology.Cosmology or None

        Returns
        -------
        SubhaloDataProducts instance
        """
        if not is_known(hi_model):
            raise ValueError(
                f"unknown hi_model {hi_model!r}; choose one of "
                "['GD14_vol', 'GD14_map', 'K13_vol', 'K13_map']"
            )
        spec = resolve(hi_model)
        is_projected = spec["kind"] == "map"
        if is_projected and uv_mode == "global":
            raise ValueError(
                f"uv_mode='global' does not apply to {spec['label']}: the "
                "projected models take their UV field from the SFR map"
            )

        dp = cls(ngrid=ngrid)

        # --- Load header ---
        header = loadHeader(basePath, snapNum)
        h = header["HubbleParam"]
        a = header["Time"]
        z = header["Redshift"]
        dp.metadata["redshift"] = float(z)
        dp.metadata["hubble_h"] = float(h)
        dp.metadata["scale_factor"] = float(a)
        dp.metadata["snapNum"] = snapNum
        dp.metadata["subID"] = subID
        dp.metadata["hi_model"] = spec["label"]
        dp.metadata["hi_reference"] = spec["reference"]

        # h and z come from this snapshot, not from the caller's config
        overrides = {"h": float(h), "z": float(z)}
        if header.get("BoxSize"):
            overrides["boxsize_ckpc_per_h"] = float(header["BoxSize"])
        cfg = HIModelConfig.build(
            hi_config,
            h=overrides["h"],
            z=overrides["z"],
            boxsize_ckpc_per_h=overrides.get("boxsize_ckpc_per_h"),
        )
        # The product grid resolution also sets the UV grid resolution, so that
        # the stored uv_field and the field that produced the masses agree.
        cfg = cfg.with_(ngrid=int(ngrid))

        # --- Cosmology for RSD ---
        if cosmology is None:
            cosmology = Planck15
        H_a = cosmology.H(z).to_value(U.km / U.s / U.Mpc)  # H(z) in km/s/Mpc

        # --- Load subhalo metadata ---
        data_sub = loadSingle(basePath, snapNum, subhaloID=subID)
        if load_subhalo_only:
            subset_g = getSnapOffsets(basePath, snapNum, subID, "Subhalo")
            subset_dm = getSnapOffsets(basePath, snapNum, subID, "Subhalo")
        else:
            haloID = int(data_sub["SubhaloGrNr"])
            subset_g = getSnapOffsets(basePath, snapNum, haloID, "Group")
            subset_dm = getSnapOffsets(basePath, snapNum, haloID, "Group")

        # --- Load gas particles ---
        gas_fields = [
            "Masses", "Density", "InternalEnergy", "NeutralHydrogenAbundance",
            "StarFormationRate", "Coordinates", "Velocities",
            "GFM_Metals", "GFM_Metallicity",
        ]
        gas_mdi = [None, None, None, None, None, None, None, 0, None]
        data_g = loadSubset(
            basePath, snapNum, "gas", fields=gas_fields,
            subset=subset_g, mdi=gas_mdi,
        )

        # --- Load DM particles ---
        # DM particles have constant mass from Header.MassTable[1]
        dm_fields = ["Coordinates", "Velocities"]
        data_dm = loadSubset(
            basePath, snapNum, "dm", fields=dm_fields, subset=subset_dm,
        )
        # Get DM particle mass from snapshot header (MassTable in snap file)
        import h5py as _h5py
        from martini.sources._illustris_tools import snapPath
        with _h5py.File(snapPath(basePath, snapNum), "r") as _snap:
            _mt = _snap["Header"].attrs.get("MassTable", [0.0, 0.0])
            dm_mass_per_particle = float(_mt[1]) * 1e10 / h  # Msun

        # --- Compute HI/H2 masses with the selected model ---
        hi_kwargs = {}
        if is_projected:
            # The projected models project the galaxy about the subhalo's own
            # position, so they need the subhalo record and the stars. The stars
            # are loaded here as well as below for the source-property grids, but
            # only the fields the model reads are required.
            star_fields_model = [
                "Masses", "Coordinates", "Velocities", "GFM_StellarFormationTime",
            ]
            data_star_model = loadSubset(
                basePath, snapNum, "star", fields=star_fields_model,
                subset=subset_g,
            )
            halfmass = np.asarray(
                data_sub["SubhaloHalfmassRadType"], dtype="float64"
            )
            hi_kwargs["stars"] = data_star_model
            hi_kwargs["sub"] = {
                "pos": np.asarray(data_sub["SubhaloPos"], dtype="float64"),
                "vel": np.asarray(data_sub["SubhaloVel"], dtype="float64"),
                "rhalf_gas": float(halfmass[0]),
                "rhalf_star": float(halfmass[4]),
            }

        hi_result = compute_hi_mass(
            data_g, spec["label"], cfg, uv_mode=uv_mode, **hi_kwargs
        )
        m_HI, m_H2 = hi_result.m_hi, hi_result.m_h2
        # TNG code units -> Msun: * 1e10 / h
        m_HI_Msun = m_HI * 1e10 / h
        m_H2_Msun = m_H2 * 1e10 / h
        dp.metadata["m_HI_Msun"] = float(np.sum(m_HI_Msun))
        dp.metadata["m_H2_Msun"] = float(np.sum(m_H2_Msun))

        # UV field grid, for the volume models only (the projected models work
        # on an SFR surface-density map instead, which has its own geometry).
        uv_edges_phys = None
        if not is_projected:
            if uv_mode == "adaptive_fft":
                # The grid is the raw convolved FFT field -- it does NOT include
                # the global-background floor that the per-particle values get,
                # and its axes are in physical kpc, unlike every other grid here.
                _umw, uv_grid, uv_edges = compute_galaxy_fft_uv(
                    data_g, cfg, return_grid=True
                )
                dp.uv_field = uv_grid
                uv_edges_phys = uv_edges
            else:
                dp.uv_field = np.full(
                    (ngrid, ngrid, ngrid),
                    calc_UMW_global_Faucher11(z),
                    dtype="float32",
                )

        # --- Get comoving coordinates ---
        # Gas coords: ckpc/h (comoving)
        gas_coords = data_g["Coordinates"][:].astype("float64")  # (N, 3) ckpc/h
        gas_vel = data_g["Velocities"][:].astype("float64")  # (N, 3) km/s * sqrt(a)
        gas_vlos = gas_vel[:, 0] * np.sqrt(a)  # LOS peculiar velocity [km/s]

        # DM coords
        dm_coords = data_dm["Coordinates"][:].astype("float64")  # ckpc/h
        dm_mass_arr = np.full(data_dm["count"], dm_mass_per_particle, dtype="float64")
        dm_vel = data_dm["Velocities"][:].astype("float64")
        dm_vlos = dm_vel[:, 0] * np.sqrt(a)  # LOS peculiar velocity [km/s]

        # --- Define comoving grid centered on DM center of mass ---
        # Compute DM center of mass (equal-mass particles → mean)
        dm_com = dm_coords.mean(axis=0)  # ckpc/h

        # Extents: max distance from DM COM to any particle (gas + DM), per axis
        all_coords = np.vstack([gas_coords, dm_coords])
        max_dist = np.max(np.abs(all_coords - dm_com), axis=0)  # per-axis max offset
        half_box = np.maximum(1.1 * max_dist, 5.0)  # at least ±5 ckpc/h per axis
        center = dm_com

        dp.metadata["dm_centre_ckpc_h"] = [float(c) for c in dm_com]

        edges = [
            np.linspace(
                center[i] - half_box[i],
                center[i] + half_box[i],
                ngrid + 1,
            )
            for i in range(3)
        ]
        dp.grid_edges = tuple(edges)
        dp.metadata["grid_bounds_ckpc_h"] = [
            float(edges[0][0]), float(edges[0][-1]),
            float(edges[1][0]), float(edges[1][-1]),
            float(edges[2][0]), float(edges[2][-1]),
        ]
        dp.metadata["grid_unit"] = "ckpc/h"
        dp.metadata["mass_unit"] = "Msun"

        # --- Redshift-space positions ---
        # s_x = x + v_los / (a H(a))  [all in comoving ckpc/h]
        # v_los [km/s] / H_a [km/s/Mpc] = physical Mpc
        # * 1000 -> physical kpc; * h / a -> comoving ckpc/h
        #
        # The velocities are taken **relative to the object's own frame** (the DM
        # centre of mass -- the same reference the box is centred on). Applying
        # the shift to absolute velocities instead leaves it dominated by the
        # bulk motion: at z = 0 one km/s is 10 ckpc/h, so a 600 km/s bulk
        # velocity displaces every particle by 6 Mpc, which is far outside a
        # galaxy-sized box. Measured on subhalo 113, the redshift-space grids
        # retained 0.5% of the mass before this was fixed.
        shift_unit = 1.0e3 * h / (a * H_a)  # ckpc/h per km/s of v_los
        v_bulk_los = float(dm_vlos.mean())
        coords_rsd = {}
        for name, coords, vlos in (
            ("dm", dm_coords, dm_vlos),
            ("gas", gas_coords, gas_vlos),
        ):
            shifted = coords.copy()
            shifted[:, 0] += (vlos - v_bulk_los) * shift_unit
            coords_rsd[name] = shifted

        # --- Redshift-space box ---
        # Same sky-plane extent, but the line-of-sight axis has to hold the
        # shifted particles as well, and that extent comes from the velocity
        # range rather than from the geometry. The real-space box above is left
        # exactly as it was: keeping the two separate is what stops RSD from
        # stretching -- and so degrading -- grids that are not about redshift
        # space at all. Note `half_box` must not be reused as a starting point
        # without this max(): the helper pads by 1.05 per axis, while the
        # real-space box above uses 1.1.
        los_shifted = np.concatenate([coords_rsd["gas"][:, 0], coords_rsd["dm"][:, 0]])
        _, los_half_rsd = _rsd_los_halfwidth(
            all_coords[:, 0], los_shifted, center[0], pad=1.1
        )
        half_box_rsd = np.array(half_box, dtype="float64")
        half_box_rsd[0] = max(float(los_half_rsd), float(half_box[0]))
        edges_rsd = [
            np.linspace(
                center[i] - half_box_rsd[i],
                center[i] + half_box_rsd[i],
                ngrid + 1,
            )
            for i in range(3)
        ]
        dp.grid_edges_rsd = tuple(edges_rsd)
        dp.metadata["grid_bounds_ckpc_h_rsd"] = [
            float(edges_rsd[0][0]), float(edges_rsd[0][-1]),
            float(edges_rsd[1][0]), float(edges_rsd[1][-1]),
            float(edges_rsd[2][0]), float(edges_rsd[2][-1]),
        ]
        dp.metadata["rsd_velocity_frame"] = "dm_com_bulk_subtracted"

        # --- Grid: DM real-space ---
        dm_real, _ = np.histogramdd(dm_coords, bins=edges, weights=dm_mass_arr)
        dp.dm_real = dm_real / _cell_volume(edges)  # Msun / (ckpc/h)^3

        # --- Grid: DM redshift-space ---
        dm_rsd, _ = np.histogramdd(coords_rsd["dm"], bins=edges_rsd, weights=dm_mass_arr)
        dp.dm_rsd = dm_rsd / _cell_volume(edges_rsd)

        # --- Grid: HI real-space ---
        hi_real, _ = np.histogramdd(gas_coords, bins=edges, weights=m_HI_Msun)
        dp.hi_real = hi_real / _cell_volume(edges)

        # --- Grid: HI redshift-space ---
        hi_rsd, _ = np.histogramdd(coords_rsd["gas"], bins=edges_rsd, weights=m_HI_Msun)
        dp.hi_rsd = hi_rsd / _cell_volume(edges_rsd)

        # --- Grid: H2 real-space ---
        h2_real, _ = np.histogramdd(gas_coords, bins=edges, weights=m_H2_Msun)
        dp.h2_real = h2_real / _cell_volume(edges)

        # --- Grid: H2 redshift-space ---
        h2_rsd, _ = np.histogramdd(coords_rsd["gas"], bins=edges_rsd, weights=m_H2_Msun)
        dp.h2_rsd = h2_rsd / _cell_volume(edges_rsd)

        # What the redshift-space box actually holds. It is sized to fit, so
        # these should be 1.0; a value below it means the box no longer follows
        # the particles and mass is being lost silently.
        vol, vol_rsd = _cell_volume(edges), _cell_volume(edges_rsd)
        for name, real_grid, rsd_grid in (
            ("dm", dp.dm_real, dp.dm_rsd),
            ("hi", dp.hi_real, dp.hi_rsd),
            ("h2", dp.h2_real, dp.h2_rsd),
        ):
            dp.metadata[f"f_mass_in_rsd_box_{name}"] = float(
                rsd_grid.sum() * vol_rsd
            ) / max(float(real_grid.sum() * vol), 1e-30)

        # --- UV field ---
        # dp.uv_field was filled in above (or left None for the projected models,
        # whose UV field lives on the projection grid, not this one).
        if dp.uv_field is not None:
            # unlike every other grid here, uv_field lives on a physical-kpc grid
            dp.metadata["uv_grid_physical_kpc"] = True
            if uv_edges_phys is not None:
                # The UV grid's own bounds, in the same
                # [x_lo, x_hi, y_lo, y_hi, z_lo, z_hi] order as
                # grid_bounds_ckpc_h -- but NOT the same grid: uv_field sits on
                # a physical-kpc box that follows the gas extent, while every
                # mass grid uses the comoving ckpc/h DM-centred box.
                bounds = []
                for axis in uv_edges_phys:
                    bounds += [float(axis[0]), float(axis[-1])]
                dp.metadata["uv_grid_bounds_physical_kpc"] = bounds

        # Convert to float32 for storage
        for attr in ["uv_field", "dm_real", "dm_rsd", "hi_real", "hi_rsd",
                     "h2_real", "h2_rsd"]:
            val = getattr(dp, attr)
            if val is not None and val.dtype == np.float64:
                setattr(dp, attr, val.astype("float32"))

        return dp

    def save_hdf5(self, filename):
        """Save all data products to an HDF5 file.

        Parameters
        ----------
        filename : str — Output HDF5 file path
        """
        import h5py

        with h5py.File(filename, "w") as f:
            # Grid edges
            for i, axis in enumerate(["x", "y", "z"]):
                f.create_dataset(
                    f"grid_edges_{axis}", data=self.grid_edges[i].astype("float32")
                )
                if self.grid_edges_rsd is not None:
                    f.create_dataset(
                        f"grid_edges_{axis}_rsd",
                        data=self.grid_edges_rsd[i].astype("float32"),
                    )

            # Data grids
            for name in ["uv_field", "dm_real", "dm_rsd", "hi_real",
                         "hi_rsd", "h2_real", "h2_rsd"]:
                val = getattr(self, name)
                if val is not None:
                    f.create_dataset(
                        name, data=val, compression="lzf",
                        chunks=(min(32, self.ngrid),) * 3,
                    )

            # Metadata
            for key, val in self.metadata.items():
                if isinstance(val, list):
                    f.attrs[key] = np.array(val)
                else:
                    f.attrs[key] = val
            f.attrs["ngrid"] = self.ngrid

    @classmethod
    def load_hdf5(cls, filename):
        """Load data products from an HDF5 file.

        Parameters
        ----------
        filename : str — Path to HDF5 file

        Returns
        -------
        SubhaloDataProducts instance
        """
        import h5py

        with h5py.File(filename, "r") as f:
            ngrid = f.attrs["ngrid"]
            dp = cls(ngrid=ngrid)

            # Grid edges
            dp.grid_edges = tuple(
                f[f"grid_edges_{axis}"][:] for axis in ["x", "y", "z"]
            )
            if "grid_edges_x_rsd" in f:
                dp.grid_edges_rsd = tuple(
                    f[f"grid_edges_{axis}_rsd"][:] for axis in ["x", "y", "z"]
                )

            # Data grids
            for name in ["uv_field", "dm_real", "dm_rsd", "hi_real",
                         "hi_rsd", "h2_real", "h2_rsd"]:
                if name in f:
                    setattr(dp, name, f[name][:])

            # Metadata
            dp.metadata = dict(f.attrs.items())

        return dp


def _rsd_los_halfwidth(los_real, los_shifted, centre, pad=1.05):
    """Line-of-sight half-widths for the real- and redshift-space boxes.

    The real-space extent follows the particles. The redshift-space one follows
    the *velocity* range, which is a different thing entirely: at z = 0 one km/s
    is 10 ckpc/h, so even a modest velocity dispersion displaces particles by
    Mpc and the redshift-space box has to be correspondingly longer. Returning
    both keeps the two boxes distinct, so that a redshift-space choice cannot
    stretch the real-space grids.

    Parameters
    ----------
    los_real, los_shifted : array-like
        Line-of-sight coordinates of the particles in real and redshift space
        [ckpc/h].
    centre : float
        Line-of-sight coordinate of the box centre [ckpc/h].
    pad : float, optional
        Margin factor applied to the extents (default 1.05).

    Returns
    -------
    real_half, rsd_half : float
        Half-widths [ckpc/h], with ``rsd_half >= real_half``.
    """
    real = pad * float(
        np.max(np.abs(np.asarray(los_real, dtype="float64") - float(centre)))
    )
    shifted = pad * float(
        np.max(np.abs(np.asarray(los_shifted, dtype="float64") - float(centre)))
    )
    return real, max(real, shifted)


def _cell_volume(edges):
    """Compute volume of histogram cells from edge arrays.

    Parameters
    ----------
    edges : tuple of 3 arrays — Bin edges

    Returns
    -------
    float — Cell volume (product of spacings for uniform grid)
    """
    dx = edges[0][1] - edges[0][0]
    dy = edges[1][1] - edges[1][0]
    dz = edges[2][1] - edges[2][0]
    return float(dx * dy * dz)
