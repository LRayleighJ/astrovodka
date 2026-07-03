"""
K13 (Krumholz 2013) volumetric HI/H2 model for IllustrisTNG simulations.

Implements the Diemer+2018 full volumetric version of the K13 model
with adaptive per-galaxy FFT UV field computation. Extracted and adapted
from the standalone K13 catalog builder script.
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator


# --- Global UV background ---

def calc_UMW_global_Faucher11(redshift):
    """
    Compute the global UV background radiation field following
    Faucher-Giguere et al. (2011).

    Returns U_MW in dimensionless units relative to the Milky Way value
    (Draine field: 3.43e-8 erg/s/cm^2).

    Parameters
    ----------
    redshift : float
        Redshift.

    Returns
    -------
    float
        U_MW, the dimensionless UV radiation field.
    """
    log10_UMW = (
        -0.004292 * redshift ** 4
        + 0.061446 * redshift ** 3
        - 0.358413 * redshift ** 2
        + 1.204787 * redshift
        - 2.846354
    )
    return 10.0 ** log10_UMW


# --- f_neutral estimation ---

def get_f_neutral_H(nH_code, sfr_code):
    """
    Estimate the neutral hydrogen fraction f_neutral = f_HI + f_H2
    from TNG snapshot fields.

    TNG NeutralHydrogenAbundance = f_HI (after Rahmati+2013 self-shielding),
    NOT f_neutral = f_HI + f_H2.

    - SF particles (SFR > 0): f_neutral = 1 (cold phase ~ fully neutral).
    - Non-SF particles: H2 negligible, f_neutral ~ f_HI (bias ~10%).

    Parameters
    ----------
    nH_code : ndarray
        TNG NeutralHydrogenAbundance field (f_HI).
    sfr_code : ndarray
        TNG StarFormationRate field.

    Returns
    -------
    ndarray
        Estimated f_neutral per particle.
    """
    f_neutral = nH_code.copy()
    f_neutral[sfr_code > 0.0] = 1.0
    return f_neutral


# --- K13 volumetric model ---

def f_mol_K13_vol(rho_code, u_code, Z, nH_code, sfr_code, X_H, U_MW,
                  redshift=0.0, hubble_h=0.6774, mgas=None):
    """
    Compute molecular hydrogen fraction using the K13 model,
    Diemer+2018 full volumetric version.

    Parameters
    ----------
    rho_code : ndarray
        TNG comoving density in code units.
    u_code : ndarray
        TNG internal energy in code units [(km/s)^2].
    Z : ndarray
        Gas metallicity (absolute mass fraction).
    nH_code : ndarray
        TNG NeutralHydrogenAbundance (f_HI).
    sfr_code : ndarray
        TNG StarFormationRate.
    X_H : ndarray
        Hydrogen mass fraction per particle.
    U_MW : ndarray or float
        UV radiation field in Draine/MW units. Can be scalar (global UVB)
        or per-particle array (adaptive FFT).
    redshift : float
        Redshift for unit conversions.
    hubble_h : float
        Dimensionless Hubble parameter.
    mgas : ndarray or None
        Gas particle masses (unused, kept for interface compatibility).

    Returns
    -------
    f_mol : ndarray
        H2 fraction of total hydrogen [0, 1].
    f_HI : ndarray
        HI fraction of total hydrogen.
    f_H2 : ndarray
        H2 fraction of total hydrogen.
    """
    from ..cosmology_utils import (
        rho_code_to_nH,
        calc_Sigma_neutral_physical,
        Z_SUN,
    )

    # Estimate f_neutral
    f_neutral = get_f_neutral_H(nH_code, sfr_code)

    # Hydrogen number density
    n_H = np.maximum(rho_code_to_nH(rho_code, redshift, hubble_h, X_H), 1e-10)
    n_CNM_gas = np.maximum(n_H * f_neutral, 1e-10)

    # Metallicity normalized to solar
    D_MW = np.maximum(Z / Z_SUN, 1e-10)

    # Two-phase CNM density (K13 Eq. C14)
    n_CNM_2p = 23.0 * U_MW * 4.1 / (1.0 + 3.1 * D_MW ** 0.365)
    n_CNM = np.maximum(n_CNM_gas, n_CNM_2p)

    # Radiation field parameter
    chi = 7.2 * U_MW * (n_CNM / 10.0) ** (-1.0)
    chi = np.clip(chi, 1e-10, 1e10)

    # Surface density via Jeans length, clumping factor f_c=5
    Sigma_neutral = calc_Sigma_neutral_physical(
        u_code, rho_code, f_neutral * X_H, redshift, hubble_h
    )
    Sigma_Msun_pc2 = Sigma_neutral * 1e-6
    f_c = 5.0
    tau_c = 0.066 * f_c * D_MW * np.maximum(Sigma_Msun_pc2, 1e-20)

    # s parameter (K13 Eq. 6)
    s = np.log(1.0 + 0.6 * chi + 0.01 * chi ** 2) / (0.6 * np.maximum(tau_c, 1e-20))

    # f_mol (K13 Eq. 2)
    f_mol = np.where(s < 2.0, 1.0 - 3.0 * s / (4.0 + s), 0.0)
    f_mol = np.clip(f_mol, 0.0, 1.0)

    f_HI = f_neutral * (1.0 - f_mol)
    f_H2 = f_neutral * f_mol

    return f_mol, f_HI, f_H2


# --- Per-particle HI/H2 mass computation ---

def compute_HI_mass_particle(gas_data, U_MW, redshift=0.0, hubble_h=0.6774):
    """
    Compute per-particle HI and H2 masses using the K13 volumetric model.

    Parameters
    ----------
    gas_data : dict
        Dictionary with TNG gas particle fields:
        'Masses', 'Density', 'InternalEnergy', 'StarFormationRate',
        'NeutralHydrogenAbundance', 'GFM_Metals', 'GFM_Metallicity'.
    U_MW : ndarray or float
        UV radiation field per particle, or scalar for global UVB.
    redshift : float
        Redshift.
    hubble_h : float
        Dimensionless Hubble parameter.

    Returns
    -------
    m_HI : ndarray
        HI mass per particle in code mass units (1e10 M_sun / h).
    m_H2 : ndarray
        H2 mass per particle in code mass units.
    """
    rho_code = gas_data['Density'][:].astype('float64')
    u_code = gas_data['InternalEnergy'][:].astype('float64')
    nH_code = gas_data['NeutralHydrogenAbundance'][:].astype('float64')
    sfr_code = gas_data['StarFormationRate'][:].astype('float64')
    # GFM_Metals may be 1D (if loaded with mdi=0) or 2D
    gfm = gas_data['GFM_Metals'][:].astype('float64')
    if gfm.ndim == 2:
        X_H = gfm[:, 0]
    else:
        X_H = gfm
    Z = gas_data['GFM_Metallicity'][:].astype('float64')
    m_gas = gas_data['Masses'][:].astype('float64')

    f_mol, f_HI, f_H2 = f_mol_K13_vol(
        rho_code, u_code, Z, nH_code, sfr_code, X_H, U_MW,
        redshift=redshift, hubble_h=hubble_h, mgas=m_gas,
    )

    m_HI = m_gas * X_H * f_HI
    m_H2 = m_gas * X_H * f_H2

    return m_HI, m_H2


# --- Adaptive per-galaxy FFT UV field ---

def compute_galaxy_fft_uv(gas_data, redshift=0.0, hubble_h=0.6774,
                          ngrid=128, f_esc=0.1, box_expand=1.1,
                          return_grid=False):
    """
    Compute the adaptive per-galaxy FFT UV radiation field.

    Builds a 128^3 SFR grid, convolves with a 1/r^2 kernel via FFT,
    then interpolates back to per-particle U_MW values in the
    Diemer+2018 / Starburst99 / Draine calibration frame.
    Always applies a floor of the global UVB.

    Parameters
    ----------
    gas_data : dict
        Dictionary with keys 'Coordinates' and 'StarFormationRate'.
    redshift : float
        Redshift.
    hubble_h : float
        Dimensionless Hubble parameter.
    ngrid : int
        Grid resolution per dimension (default 128).
    f_esc : float
        Escape fraction for UV photons (default 0.1).
    box_expand : float
        Multiplier for adaptive box size (default 1.1).

    Returns
    -------
    ndarray
        Per-particle U_MW (dimensionless, in Draine units),
        with a floor of the global UVB.
    """
    coords = gas_data['Coordinates'][:].astype('float64')  # ckpc/h
    sfr = gas_data['StarFormationRate'][:].astype('float64')

    # TNG coords: ckpc/h -> physical kpc
    coords_kpc = coords / hubble_h / (1.0 + redshift)

    # Adaptive box
    coord_min = coords_kpc.min(axis=0)
    coord_max = coords_kpc.max(axis=0)
    coord_range = coord_max - coord_min
    L_box = box_expand * np.max(coord_range)
    L_box = max(L_box, 10.0)
    center = (coord_min + coord_max) / 2.0
    dx = L_box / ngrid
    half_box = L_box / 2.0

    # Build SFR grid
    sf_mask = sfr > 0.0
    edges = [
        np.linspace(center[i] - half_box, center[i] + half_box, ngrid + 1)
        for i in range(3)
    ]

    if np.sum(sf_mask) > 0:
        sfr_grid, _ = np.histogramdd(
            coords_kpc[sf_mask], bins=edges, weights=sfr[sf_mask]
        )
    else:
        sfr_grid = np.zeros((ngrid, ngrid, ngrid))

    # Build 1/r^2 kernel in real space
    x_1d = (np.arange(ngrid) - ngrid // 2) * dx
    XK, YK, ZK = np.meshgrid(x_1d, x_1d, x_1d, indexing='ij')
    RK = np.sqrt(XK ** 2 + YK ** 2 + ZK ** 2)
    RK_soft = np.maximum(RK, dx / 4.0)

    # SFR-to-FUV conversion: 3.3e-6 photons/s/cm^2/Hz per M_sun/yr
    # (calibrated to Starburst99 models)
    K_grid_real = 3.3e-6 * (1.0 / RK_soft) ** 2

    # FFT convolution
    sfr_k = np.fft.fftn(sfr_grid)
    K_k = np.fft.fftn(np.fft.ifftshift(K_grid_real))
    F_grid = np.fft.ifftn(sfr_k * K_k).real

    # U_MW (Diemer+2018 calibration frame)
    # 3.43e-8 is the Draine field normalization in erg/s/cm^2
    UMW_raw = f_esc * F_grid / 3.43e-8

    # Interpolate to particle positions
    x_c = edges[0][:-1] + dx / 2.0
    interp = RegularGridInterpolator(
        (x_c, x_c, x_c), UMW_raw,
        bounds_error=False, fill_value=0.0,
    )

    coords_clipped = np.clip(
        coords_kpc,
        center - half_box + dx,
        center + half_box - dx,
    )
    UMW_per_particle = interp(coords_clipped)

    # Apply a floor of the global UVB so U_MW is never zero
    UMW_UVB = calc_UMW_global_Faucher11(redshift)
    UMW_per_particle = UMW_per_particle + UMW_UVB

    if return_grid:
        return UMW_per_particle, UMW_raw, edges
    return UMW_per_particle
