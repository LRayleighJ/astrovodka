"""
3D gridded data products for TNG subhalos.

Provides SubhaloDataProducts — a container for 3D arrays of UV field,
DM/HI/H2 mass density (comoving real-space and redshift-space), with
HDF5 save/load support.
"""

import numpy as np
import astropy.units as U
import astropy.constants as const
from astropy.cosmology import Planck15

from .sources._k13_model import compute_HI_mass_particle, compute_galaxy_fft_uv
from .sources._illustris_tools import loadHeader, loadSingle, loadSubset, getSnapOffsets


class SubhaloDataProducts:
    """Container for 3D gridded data products of a TNG subhalo.

    Contains 7 grids of shape (ngrid, ngrid, ngrid):
    - uv_field: real-space UV radiation field [Draine units]
    - dm_real, dm_rsd: comoving DM mass density [Msun / (ckpc/h)^3]
    - hi_real, hi_rsd: comoving HI mass density [Msun / (ckpc/h)^3]
    - h2_real, h2_rsd: comoving H2 mass density [Msun / (ckpc/h)^3]

    Real-space uses comoving coordinates (ckpc/h).
    Redshift-space (RSD) shifts the LOS coordinate by v_los/(a*H(a)).
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
        self.grid_edges = None  # 3-tuple of 1D edge arrays [ckpc/h]
        self.metadata = {}

    @classmethod
    def from_subhalo(
        cls,
        basePath,
        snapNum,
        subID,
        ngrid=128,
        uv_mode="global",
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
        uv_mode : str — "global" (UVB) or "adaptive_fft"
        load_subhalo_only : bool — Load subhalo-only (True) or full FOF group
        cosmology : astropy.cosmology.Cosmology or None

        Returns
        -------
        SubhaloDataProducts instance
        """
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
        from .sources._illustris_tools import snapPath
        with _h5py.File(snapPath(basePath, snapNum), "r") as _snap:
            _mt = _snap["Header"].attrs.get("MassTable", [0.0, 0.0])
            dm_mass_per_particle = float(_mt[1]) * 1e10 / h  # Msun

        # --- Compute K13 HI/H2 masses ---
        from .sources._k13_model import calc_UMW_global_Faucher11
        U_MW_global = calc_UMW_global_Faucher11(z)

        if uv_mode == "adaptive_fft":
            U_MW, uv_grid, uv_edges = compute_galaxy_fft_uv(
                data_g, redshift=z, hubble_h=h, ngrid=ngrid, return_grid=True,
            )
            dp.uv_field = uv_grid  # Physical kpc grid
            uv_edges_phys = uv_edges  # edges in physical kpc
        else:
            U_MW = U_MW_global
            uv_edges_phys = None

        m_HI, m_H2 = compute_HI_mass_particle(
            data_g, U_MW, redshift=z, hubble_h=h,
        )
        # TNG code units -> Msun: * 1e10 / h
        m_HI_Msun = m_HI * 1e10 / h
        m_H2_Msun = m_H2 * 1e10 / h

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

        # --- Define comoving grid ---
        # Use gas + DM coords for bounds
        all_coords = np.vstack([gas_coords, dm_coords])
        coord_min = all_coords.min(axis=0)
        coord_max = all_coords.max(axis=0)
        coord_range = coord_max - coord_min
        # Per-axis box: 1.1× each axis's range, with 10 ckpc/h floor
        box_size = np.maximum(1.1 * coord_range, 10.0)
        center = (coord_min + coord_max) / 2.0

        edges = [
            np.linspace(
                center[i] - box_size[i] / 2.0,
                center[i] + box_size[i] / 2.0,
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

        # --- Grid: DM real-space ---
        dm_real, _ = np.histogramdd(dm_coords, bins=edges, weights=dm_mass_arr)
        dp.dm_real = dm_real / _cell_volume(edges)  # Msun / (ckpc/h)^3

        # --- Grid: DM redshift-space ---
        # RSD: s_x = x + v_los / (a * H(a))  [all in comoving ckpc/h]
        # v_los [km/s] / H_a [km/s/Mpc] = physical Mpc
        # * 1000 -> physical kpc; * h / a -> comoving ckpc/h
        rsd_shift_dm = dm_vlos / (a * H_a) * 1000.0 * h  # ckpc/h
        dm_coords_rsd = dm_coords.copy()
        dm_coords_rsd[:, 0] += rsd_shift_dm
        dm_rsd, _ = np.histogramdd(dm_coords_rsd, bins=edges, weights=dm_mass_arr)
        dp.dm_rsd = dm_rsd / _cell_volume(edges)

        # --- Grid: HI real-space ---
        hi_real, _ = np.histogramdd(gas_coords, bins=edges, weights=m_HI_Msun)
        dp.hi_real = hi_real / _cell_volume(edges)

        # --- Grid: HI redshift-space ---
        rsd_shift_gas = gas_vlos / (a * H_a) * 1000.0 * h
        gas_coords_rsd = gas_coords.copy()
        gas_coords_rsd[:, 0] += rsd_shift_gas
        hi_rsd, _ = np.histogramdd(gas_coords_rsd, bins=edges, weights=m_HI_Msun)
        dp.hi_rsd = hi_rsd / _cell_volume(edges)

        # --- Grid: H2 real-space ---
        h2_real, _ = np.histogramdd(gas_coords, bins=edges, weights=m_H2_Msun)
        dp.h2_real = h2_real / _cell_volume(edges)

        # --- Grid: H2 redshift-space ---
        h2_rsd, _ = np.histogramdd(gas_coords_rsd, bins=edges, weights=m_H2_Msun)
        dp.h2_rsd = h2_rsd / _cell_volume(edges)

        # --- UV field ---
        if uv_mode == "adaptive_fft" and dp.uv_field is not None:
            # uv_field was computed in physical kpc grid; store center positions
            dp.metadata["uv_grid_physical_kpc"] = True
        elif uv_mode == "global":
            dp.uv_field = np.full((ngrid, ngrid, ngrid), U_MW_global, dtype="float32")

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

            # Data grids
            for name in ["uv_field", "dm_real", "dm_rsd", "hi_real",
                         "hi_rsd", "h2_real", "h2_rsd"]:
                if name in f:
                    setattr(dp, name, f[name][:])

            # Metadata
            dp.metadata = dict(f.attrs.items())

        return dp


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
