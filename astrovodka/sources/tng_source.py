"""
Provides the TNGSource class for astrovodka.

Loads IllustrisTNG simulation data from local disk and computes HI/H2
masses using one of four Diemer+2018 models. Uses cosmological distances
from astropy.cosmology (default Planck15).

Unlike MARTINI's TNGSource, this class:
1. Reads TNG data from a local directory (no web API dependency)
2. Computes HI/H2 masses with a physically motivated model instead of a
   nearby-galaxy pressure calibration
3. Derives distances and velocities from the simulation redshift + cosmology
"""

from __future__ import annotations

import numpy as np
from astropy import units as U
from astropy.coordinates import ICRS
from scipy.spatial.transform import Rotation

from ..L_coords import L_coords
from ..cosmology_utils import (
    angular_diameter_distance,
    get_cosmology,
    luminosity_distance,
)
from ..himodel import (
    DEFAULT_METHOD,
    HIModelConfig,
    compute_hi_mass,
    is_known,
    resolve,
)
from ..himodel.physics import calc_temperature
from ..sph_kernels import _CubicSplineKernel, find_fwhm
from .sph_source import SPHSource
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from astropy.coordinates import BaseRADecFrame
    from astropy.cosmology import Cosmology


def _m200c_msun(basePath, snapNum, haloID, h):
    """M_200c of the parent FOF group in Msun, or None if unavailable.

    Only the projected K13 model uses this, and only when
    ``HIModelConfig(rho_dm_mode='nfw')`` asks for it (the default is ``'none'``,
    because the paper notes rho_dm is subdominant wherever molecular gas
    matters). The group value is a fair proxy for a central subhalo and an
    overestimate for a satellite, which is worth knowing before switching the
    option on.
    """
    try:
        from ._illustris_tools import loadSingle

        group = loadSingle(basePath, snapNum, haloID=int(haloID))
        m200 = group.get("Group_M_Crit200")
        if m200 is None or not np.isfinite(m200) or m200 <= 0:
            return None
        return float(m200) * 1e10 / h
    except Exception:
        return None


class TNGSource(SPHSource):
    """
    HI source for IllustrisTNG simulations with Diemer+2018 HI modeling.

    Loads gas particle data from local TNG snapshot files, computes
    per-particle HI and H2 masses with one of four selectable models, and
    initializes the SPHSource with cosmological distances.

    Parameters
    ----------
    basePath : str
        Path to the TNG simulation output directory, i.e. the directory that
        contains ``groups_<snapNum>/`` and ``snapdir_<snapNum>/``
        (e.g. ``'/path/to/TNG/TNG-100/'``).

    snapNum : int
        Snapshot number (e.g., 99 for z=0, 67 for z=0.5, 50 for z=2).

    subID : int
        Subhalo ID of the target object. All particles in the parent
        FOF group are used.

    cosmology : astropy.cosmology.Cosmology, optional
        Cosmology instance. Defaults to Planck15.

    hi_model : str, optional
        Which HI/H2 model to use. One of:

        ``'K13_vol'`` (default)
            Krumholz+2013, per gas cell (Diemer+2018 Appendix C.4).
        ``'K13_map'``
            K13 on a face-on projection of the galaxy (§2.3).
        ``'GD14_vol'``
            Gnedin & Draine 2014, per gas cell (Appendix C.3).
        ``'GD14_map'``
            GD14 on the face-on projection.

        The names are case-insensitive, and the bare model names ``'K13'``
        and ``'GD14'`` mean the volume variant. The projected models need the
        stellar particles as well as the gas, so they do more I/O.

    hi_config : HIModelConfig or dict, optional
        Numerical settings and cosmology for the HI model. Anything not
        supplied takes its documented default; ``h`` and ``z`` are always
        overridden with the snapshot header's values, because mixing a
        snapshot's data with another snapshot's units is never what the caller
        wants. Pass a dict to override individual fields, or a full
        ``HIModelConfig`` to control everything else. The ``HIModelConfig``
        actually used is kept as ``self.hi_config``.

    uv_mode : str, optional
        How the volume models get their UV radiation field:

        ``'adaptive_fft'`` (default)
            A per-galaxy FFT convolution of the SFR distribution, floored at
            the Faucher-Giguere+2011 background (Diemer+2018 §2.4).
        ``'global'``
            Just the uniform Faucher-Giguere+2011 background.

        Ignored by the projected models, whose UV field is part of the model
        definition. Note the default changed in v0.3.0: the previous default
        was ``'global'``, and the previous ``'adaptive_fft'`` was silently
        identical to it (its interpolator was built with one axis' grid reused
        for all three, so the FFT field came back as zero everywhere and only
        the background survived).

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
        hi_model: str = DEFAULT_METHOD,
        hi_config: "HIModelConfig | dict | None" = None,
        uv_mode: str = "adaptive_fft",
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

        # --- Resolve the HI model before touching the disk: a typo should not
        # cost a multi-GB load.
        if not is_known(hi_model):
            raise ValueError(
                f"unknown hi_model {hi_model!r}; choose one of "
                f"['GD14_vol', 'GD14_map', 'K13_vol', 'K13_map'] "
                "(or the bare model name 'GD14'/'K13' for the volume variant)"
            )
        self.hi_model_spec = resolve(hi_model)
        self.hi_model = self.hi_model_spec["label"]
        self.is_projected = self.hi_model_spec["kind"] == "map"
        if uv_mode not in ("adaptive_fft", "global"):
            raise ValueError(
                f"uv_mode must be 'adaptive_fft' or 'global', got {uv_mode!r}"
            )
        if self.is_projected and uv_mode == "global":
            raise ValueError(
                f"uv_mode='global' does not apply to {self.hi_model}: the "
                "projected models derive their UV field from the SFR map, so "
                "there is no particle-level field to replace. Drop the "
                "argument, or use a *_vol model to control the UV field."
            )
        self.uv_mode = uv_mode

        # --- Load snapshot header ---
        print(f"Loading TNG snapshot {snapNum}, subhalo {subID}...")
        data_header = loadHeader(basePath, snapNum)
        h = data_header["HubbleParam"]
        a = data_header["Time"]
        z = data_header["Redshift"]

        print(f"  Redshift: {z:.4f}, h: {h:.4f}, scale factor: {a:.4f}")

        # h and z always come from this snapshot, whatever the caller passed:
        # the model converts TNG code units, so using another snapshot's values
        # would silently scale every mass and density.
        self.hi_config = HIModelConfig.build(
            hi_config, h=h, z=z, boxsize_ckpc_per_h=data_header.get("BoxSize")
        )
        print(f"  HI model: {self.hi_model}  ({self.hi_config.summary()})")

        # Use snapshot redshift for physics (HI model), observed_redshift for distance
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
            "Coordinates",
            "CenterOfMass",
            "GFM_Metals",
            "GFM_Metallicity",
        ]
        mdi_full = [None, None, None, None, None, None, None, None, None, 0, None]

        print(f"  Loading gas particles ({'subhalo' if load_subhalo_only else 'FOF group'})...")
        data_g = loadSubset(
            basePath, snapNum, "gas",
            fields=full_fields_g,
            subset=subset_g,
            mdi=mdi_full,
        )
        n_gas = data_g["count"]
        print(f"  Loaded {n_gas} gas particles")

        # --- Load DM particles and compute DM center of mass ---
        import h5py as _h5py
        from ._illustris_tools import snapPath

        print("  Loading DM particles...")
        if load_subhalo_only:
            subset_dm = getSnapOffsets(basePath, snapNum, subID, "Subhalo")
        else:
            subset_dm = getSnapOffsets(basePath, snapNum, haloID, "Group")
        dm_fields = ["Coordinates", "Velocities"]
        data_dm = loadSubset(
            basePath, snapNum, "dm",
            fields=dm_fields,
            subset=subset_dm,
        )
        n_dm = data_dm["count"]
        print(f"  Loaded {n_dm} DM particles")

        dm_coords = data_dm["Coordinates"][:].astype("float64")  # ckpc/h
        dm_vel = data_dm["Velocities"][:].astype("float64")  # km/s * sqrt(a)

        # DM particle mass from snapshot MassTable (constant for all DM particles)
        with _h5py.File(snapPath(basePath, snapNum), "r") as _snap:
            _mt = _snap["Header"].attrs.get("MassTable", [0.0] * 6)
            dm_mass_per_particle = float(_mt[1]) * 1e10 / h  # Msun

        # DM COM: all DM particles have equal mass → plain mean
        dm_com_ckpc = dm_coords.mean(axis=0)  # ckpc/h comoving
        dm_com_vel = dm_vel.mean(axis=0)  # km/s * sqrt(a)

        # Convert to physical kpc
        dm_com_phys = dm_com_ckpc * (a / h)  # physical kpc
        dm_com_vel_phys = dm_com_vel * np.sqrt(a)  # physical km/s

        print(f"  DM center of mass [ckpc/h]: "
              f"[{dm_com_ckpc[0]:.1f}, {dm_com_ckpc[1]:.1f}, {dm_com_ckpc[2]:.1f}]")

        # --- Load star particles ---
        print("  Loading star particles...")
        if load_subhalo_only:
            subset_star = getSnapOffsets(basePath, snapNum, subID, "Subhalo")
        else:
            subset_star = getSnapOffsets(basePath, snapNum, haloID, "Group")
        star_fields = ["Masses", "Coordinates", "Velocities",
                       "GFM_StellarFormationTime", "GFM_InitialMass",
                       "GFM_Metallicity", "GFM_Metals"]
        # GFM_Metals: 0 -> H fraction, GFM_Metallicity: raw
        star_mdi = [None, None, None, None, None, None, 0]
        data_star = loadSubset(
            basePath, snapNum, "star",
            fields=star_fields,
            subset=subset_star,
            mdi=star_mdi,
        )
        n_star = data_star["count"]
        print(f"  Loaded {n_star} star particles")

        # --- Compute HI/H2 masses with the selected model ---
        print(f"  Computing HI/H2 masses (model={self.hi_model})...")

        hi_kwargs = {}
        if self.is_projected:
            # The projected models register the subhalo's own position, velocity
            # and half-mass radii: they build a face-on map centred on the
            # subhalo, which is a different frame from the cube's DM-centred one.
            halfmass = np.asarray(data_sub["SubhaloHalfmassRadType"], dtype="float64")
            hi_kwargs["sub"] = {
                "pos": np.asarray(data_sub["SubhaloPos"], dtype="float64"),
                "vel": np.asarray(data_sub["SubhaloVel"], dtype="float64"),
                "rhalf_gas": float(halfmass[0]),
                "rhalf_star": float(halfmass[4]),
                "m200c_msun": _m200c_msun(basePath, snapNum, haloID, h),
            }
            hi_kwargs["stars"] = data_star
            print(f"    Projecting {subID}: R_half(gas)={halfmass[0]:.1f} ckpc/h, "
                  f"R_half(star)={halfmass[4]:.1f} ckpc/h")

        hi_result = compute_hi_mass(
            data_g,
            self.hi_model,
            self.hi_config,
            uv_mode=uv_mode,
            **hi_kwargs,
        )
        m_HI, m_H2 = hi_result.m_hi, hi_result.m_h2

        m_HI_sum = np.sum(m_HI) / h  # code units -> 1e10 Msun
        m_H2_sum = np.sum(m_H2) / h
        print(f"    Total HI mass: {m_HI_sum:.4e} x 1e10 Msun")
        print(f"    Total H2 mass: {m_H2_sum:.4e} x 1e10 Msun")
        if self.is_projected:
            diag = hi_result.diagnostics
            print(
                f"    Kept {diag['n_cells_kept']}/{diag['n_cells_total']} cells "
                f"inside the projection box (L={diag['L_map_kpc']:.1f} kpc, "
                f"pixel={diag['pix_kpc'] * 1e3:.0f} pc); "
                f"f_gas_in_map={diag['f_gas_in_map']:.3f}; "
                f"align_success={diag['align_success']} via {diag['align_branch']}"
            )
            if diag["f_gas_in_map"] < 0.9:
                print(f"    WARNING: only {100 * diag['f_gas_in_map']:.1f}% of the "
                      "gas mass is inside the projection box; the projected model "
                      "assigns zero HI to the rest. Use a volume model, or set "
                      "HIModelConfig(map_clip_outside=True), if the cube needs "
                      "the whole source.")
        else:
            u = hi_result.u_mw
            print(f"    UV field ({uv_mode}): U_MW range [{u.min():.3f}, {u.max():.3f}]")

        # --- Compute physical properties ---
        X_H_g = data_g["GFM_Metals"]  # hydrogen mass fraction (1D, mdi=0)
        xe_g = data_g["ElectronAbundance"]

        # Temperature from the internal energy (the same expression MARTINI's TNG
        # source uses; shared with the group products through himodel.physics).
        T_g = calc_temperature(data_g["InternalEnergy"], X_H_g, xe_g) << U.K

        # Mass in physical units (Msun)
        m_g = data_g["Masses"] << 1e10 / h * U.Msun

        # Coordinates: comoving -> physical
        # Use "Coordinates" field if available, fall back to "CenterOfMass"
        if "Coordinates" in data_g:
            coords_raw = data_g["Coordinates"]
        else:
            coords_raw = data_g["CenterOfMass"]
        xyz_g = coords_raw * (a / h) << U.kpc
        vxyz_g = data_g["Velocities"] * np.sqrt(a) << U.km / U.s

        # Center on DM center of mass (not subhalo position)
        xyz_centre = dm_com_phys << U.kpc
        vxyz_centre = dm_com_vel_phys << U.km / U.s
        xyz_g = xyz_g - xyz_centre
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

        # ==================================================================
        # Store all physical attributes
        # ==================================================================

        # --- HI / H2 ---
        self.m_H2_g = m_H2 * (1e10 / h) << U.Msun
        self.m_HI_total = np.sum(mHI_g)  # total HI mass [Msun]
        # Per-gas-cell molecular and neutral fractions, and the diagnostics of
        # the model evaluation (UV field, alignment, projection truncation).
        self.hi_result = hi_result
        self.f_mol_g = hi_result.f_mol
        self.f_neutral_g = hi_result.cells.f_neutral
        self.u_mw_g = hi_result.u_mw
        self.hi_diagnostics = hi_result.diagnostics

        # --- Gas ---
        self.n_gas = n_gas
        self.m_gas_total = np.sum(m_g)  # total gas mass [Msun]
        self.m_gas = m_g  # per-particle gas mass [Msun]
        gas_sfr = data_g["StarFormationRate"][:].astype("float64")  # code units
        self.sfr_total = float(np.sum(gas_sfr))  # total SFR [code units Msun/yr]
        self.sfr_g = gas_sfr  # per-particle SFR

        # --- DM ---
        self.n_dm = n_dm
        self.m_dm_total = n_dm * dm_mass_per_particle << U.Msun
        self.m_dm_per_particle = dm_mass_per_particle << U.Msun
        self.dm_coords = dm_coords  # comoving ckpc/h
        self.dm_vel = dm_vel  # km/s * sqrt(a)

        # --- Stars ---
        self.n_star = n_star
        star_masses = data_star["Masses"][:].astype("float64")  # code units
        self.m_star_total = np.sum(star_masses) * 1e10 / h << U.Msun
        self.m_star_g = star_masses * 1e10 / h << U.Msun  # per-particle
        self.star_coords = data_star["Coordinates"][:].astype("float64")  # ckpc/h
        self.star_vel = data_star["Velocities"][:].astype("float64")  # km/s * sqrt(a)
        if n_star > 0:
            self.star_age = data_star["GFM_StellarFormationTime"][:].astype("float64")
            self.star_initial_mass = data_star["GFM_InitialMass"][:].astype("float64")
            self.star_metallicity = data_star["GFM_Metallicity"][:].astype("float64")
            self.star_X_H = data_star["GFM_Metals"][:].astype("float64")  # H fraction
        else:
            self.star_age = np.array([])
            self.star_initial_mass = np.array([])
            self.star_metallicity = np.array([])
            self.star_X_H = np.array([])

        # --- Group catalog metadata ---
        self.subhalo_mass = float(data_sub.get("SubhaloMass", np.nan))
        self.subhalo_sfr = float(data_sub.get("SubhaloSFR", np.nan))
        self.subhalo_mass_type = data_sub.get("SubhaloMassType", np.zeros(6))

        # --- Cosmology / distances ---
        self.source_redshift = z
        self.snapshot_hubble = h
        self.snapshot_scale_factor = a
        self.D_L = D_L
        self.D_A = D_A
        self.cosmology = cosmo

        # --- DM center of mass ---
        self.dm_centre = xyz_centre  # physical kpc
        self.dm_velocity = vxyz_centre  # physical km/s
        self.dm_centre_ckpc = dm_com_ckpc  # comoving ckpc/h

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
