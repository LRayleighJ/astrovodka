"""
Provides the TNGSource class for astrovodka.

Loads IllustrisTNG simulation data from local disk and computes HI/H2
masses using the K13 (Krumholz 2013) volumetric model. Uses cosmological
distances from astropy.cosmology (default Planck15).

Unlike MARTINI's TNGSource, this class:
1. Reads TNG data from a local directory (no web API dependency)
2. Computes HI/H2 masses using the K13 model instead of a simple pressure model
3. Derives distances and velocities from the simulation redshift + cosmology
"""

import os
import numpy as np
from scipy.spatial.transform import Rotation
from astropy import units as U, constants as C
from astropy.coordinates import ICRS
from .sph_source import SPHSource
from ..sph_kernels import _CubicSplineKernel, find_fwhm
from ..L_coords import L_coords
from ..cosmology_utils import (
    luminosity_distance,
    angular_diameter_distance,
    get_cosmology,
)
from ._k13_model import (
    compute_HI_mass_particle,
    compute_galaxy_fft_uv,
    calc_UMW_global_Faucher11,
)
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from astropy.coordinates import BaseRADecFrame
    from astropy.cosmology import Cosmology


class TNGSource(SPHSource):
    """
    HI source for IllustrisTNG simulations with K13 HI modeling.

    Loads gas particle data from local TNG snapshot files, computes
    per-particle HI and H2 masses using the K13 volumetric model, and
    initializes the SPHSource with cosmological distances.

    Parameters
    ----------
    basePath : str
        Path to the TNG simulation output directory (e.g.,
        '/data/hw2/projects/FAST/HD2_high/TNG/TNG-100/').

    snapNum : int
        Snapshot number (e.g., 99 for z=0, 67 for z=0.5, 50 for z=2).

    subID : int
        Subhalo ID of the target object. All particles in the parent
        FOF group are used.

    cosmology : astropy.cosmology.Cosmology, optional
        Cosmology instance. Defaults to Planck15.

    uv_mode : str, optional
        UV field mode: 'global' (default) uses a global UVB from
        Faucher-Giguere+2011, 'adaptive_fft' computes a per-galaxy
        adaptive FFT UV field from the SFR distribution.

    observed_redshift : float, optional
        Observed redshift for cosmological distance calculations.
        If None (default), uses the snapshot redshift. Set this to
        place the source at a different cosmological distance while
        keeping the simulation physics (K13 unit conversions) at
        the snapshot redshift.

    load_subhalo_only : bool, optional
        If True (default), load only the subhalo particles. If False,
        load the entire parent FOF group (MARTINI behavior). Subhalo-only
        loading is much faster but may miss foreground/background gas
        from other subhalos in the same group.

    vpeculiar : ~astropy.units.Quantity, optional
        Source peculiar velocity along the direction to the source centre.

    rotation : ~scipy.spatial.transform.Rotation, optional
        A rotation to apply to the source particles.

    L_coords : L_coords, optional
        Named tuple specifying inclination, position angle, and azimuthal
        rotation for angular momentum alignment.

    ra : ~astropy.units.Quantity, optional
        Right ascension for the source centroid.

    dec : ~astropy.units.Quantity, optional
        Declination for the source centroid.

    coordinate_frame : BaseRADecFrame, optional
        The coordinate frame for converting particle coordinates.
    """

    def __init__(
        self,
        basePath: str,
        snapNum: int,
        subID: int,
        *,
        cosmology: "Cosmology | None" = None,
        uv_mode: str = "global",
        observed_redshift: float | None = None,
        load_subhalo_only: bool = True,
        vpeculiar: U.Quantity[U.km / U.s] = 0 * U.km / U.s,
        rotation: Rotation | None = None,
        L_coords: L_coords | None = None,
        ra: U.Quantity[U.deg] = 0.0 * U.deg,
        dec: U.Quantity[U.deg] = 0.0 * U.deg,
        coordinate_frame: "BaseRADecFrame" = ICRS(),
    ):
        from ._illustris_tools import (
            loadHeader,
            loadSingle,
            loadSubset,
            getSnapOffsets,
        )

        # --- Load snapshot header ---
        print(f"Loading TNG snapshot {snapNum}, subhalo {subID}...")
        data_header = loadHeader(basePath, snapNum)
        h = data_header["HubbleParam"]
        a = data_header["Time"]
        z = data_header["Redshift"]

        print(f"  Redshift: {z:.4f}, h: {h:.4f}, scale factor: {a:.4f}")

        # Use snapshot redshift for physics (K13 model), observed_redshift for distance
        z_obs = observed_redshift if observed_redshift is not None else z
        if observed_redshift is not None:
            print(f"  Observed redshift (for distance): {z_obs:.4f}")

        # --- Set up cosmology ---
        cosmo = get_cosmology(cosmology)
        D_L = luminosity_distance(z_obs, cosmo)
        D_A = angular_diameter_distance(z_obs, cosmo)
        print(f"  Luminosity distance: {D_L:.1f}")
        print(f"  Angular diameter distance: {D_A:.1f}")

        # --- Load subhalo + group data ---
        data_sub = loadSingle(basePath, snapNum, subhaloID=subID)
        haloID = int(data_sub["SubhaloGrNr"])
        if load_subhalo_only:
            subset_g = getSnapOffsets(basePath, snapNum, subID, "Subhalo")
        else:
            subset_g = getSnapOffsets(basePath, snapNum, haloID, "Group")

        # Full gas fields needed for K13 + coordinate computation
        full_fields_g = [
            "Masses",
            "Velocities",
            "InternalEnergy",
            "ElectronAbundance",
            "NeutralHydrogenAbundance",
            "StarFormationRate",
            "Density",
            "CenterOfMass",
            "GFM_Metals",
            "GFM_Metallicity",
        ]
        mdi_full = [None, None, None, None, None, None, None, None, 0, None]

        print(f"  Loading gas particles ({'subhalo' if load_subhalo_only else 'FOF group'})...")
        data_g = loadSubset(
            basePath, snapNum, "gas",
            fields=full_fields_g,
            subset=subset_g,
            mdi=mdi_full,
        )
        n_gas = data_g["count"]
        print(f"  Loaded {n_gas} gas particles")

        # --- Compute K13 HI/H2 masses ---
        print(f"  Computing K13 HI/H2 masses (uv_mode={uv_mode})...")

        if uv_mode == "adaptive_fft":
            U_MW = compute_galaxy_fft_uv(data_g, redshift=z, hubble_h=h)
            print(f"    Adaptive FFT UV field: U_MW range "
                  f"[{U_MW.min():.3f}, {U_MW.max():.3f}]")
        else:
            U_MW = calc_UMW_global_Faucher11(z)
            print(f"    Global UVB: U_MW = {U_MW:.4f}")

        m_HI, m_H2 = compute_HI_mass_particle(data_g, U_MW, redshift=z, hubble_h=h)
        m_HI_sum = np.sum(m_HI) / h  # code units -> 1e10 Msun
        m_H2_sum = np.sum(m_H2) / h
        print(f"    Total HI mass: {m_HI_sum:.4e} x 1e10 Msun")
        print(f"    Total H2 mass: {m_H2_sum:.4e} x 1e10 Msun")

        # --- Compute physical properties ---
        X_H_g = data_g["GFM_Metals"]  # hydrogen mass fraction (1D, mdi=0)
        xe_g = data_g["ElectronAbundance"]
        rho_g = data_g["Density"] << 1e10 / h * U.Msun * (a / h * U.kpc) ** -3
        u_g = data_g["InternalEnergy"] << (U.km / U.s) ** 2

        # Temperature from internal energy (same formula as MARTINI TNGSource)
        mu_g = 4 / (1 + 3 * X_H_g + 4 * X_H_g * xe_g)
        gamma = 5.0 / 3.0
        T_g = (gamma - 1) / C.k_B * (u_g * mu_g * C.m_p) << U.K

        # Mass in physical units (Msun)
        m_g = data_g["Masses"] << 1e10 / h * U.Msun

        # Coordinates: comoving -> physical
        xyz_g = data_g["CenterOfMass"] * (a / h) << U.kpc
        vxyz_g = data_g["Velocities"] * np.sqrt(a) << U.km / U.s

        # Center on subhalo
        xyz_centre = data_sub["SubhaloPos"] * (a / h) << U.kpc
        xyz_g = xyz_g - xyz_centre
        vxyz_centre = data_sub["SubhaloVel"] * np.sqrt(a) << U.km / U.s
        vxyz_g = vxyz_g - vxyz_centre

        # Smoothing lengths from Voronoi cell volumes
        V_cell = (
            data_g["Masses"] / data_g["Density"] * (a / h * U.kpc) ** 3
        )
        r_cell = (3.0 * V_cell / (4.0 * np.pi)) ** (1.0 / 3.0)
        r_cell = r_cell.to(U.kpc)
        hsm_g = 2.5 * r_cell * find_fwhm(_CubicSplineKernel().kernel)

        # Convert m_HI from code units to Msun
        mHI_g = m_HI * (1e10 / h) << U.Msun

        # Store H2 mass for reference
        self.m_H2_g = m_H2 * (1e10 / h) << U.Msun
        self.source_redshift = z
        self.snapshot_hubble = h
        self.D_L = D_L
        self.D_A = D_A

        # --- Call SPHSource with redshift mode ---
        super().__init__(
            redshift=z_obs,
            cosmology=cosmo,
            vpeculiar=vpeculiar,
            rotation=rotation,
            L_coords=L_coords,
            ra=ra,
            dec=dec,
            h=h,
            T_g=T_g,
            mHI_g=mHI_g,
            xyz_g=xyz_g,
            vxyz_g=vxyz_g,
            hsm_g=hsm_g,
        )

        print(f"  TNGSource initialized: {self.npart} particles")
        print(f"  vsys = {self.vsys.to(U.km/U.s):.1f}")
