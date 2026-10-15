"""
Cosmology utilities for astrovodka.

Provides helper functions for computing cosmological distances using
astropy.cosmology, converting TNG code units to physical units, and
computing systemic velocities from redshift for the HI 21-cm line.
"""

import numpy as np
import astropy.units as U
import astropy.constants as const
from astropy.cosmology import Planck15, Cosmology

# HI rest frequency
HIfreq: U.Quantity[U.Hz] = 1.420405751e9 * U.Hz

# Physical constants in CGS
PROTONMASS = const.m_p.cgs.value
MSUN = const.M_sun.cgs.value
MPC = (1.0 * U.Mpc).to(U.cm).value
BOLTZMANN = const.k_B.cgs.value
GAMMA = 5.0 / 3.0
GAMMA_MINUS1 = GAMMA - 1.0
HYDROGEN_MASSFRAC = 0.76
Z_SUN = 0.0127
kpc_to_cm = 3.08568e21


def get_cosmology(cosmo=None):
    """
    Get a cosmology instance, defaulting to Planck15.

    Parameters
    ----------
    cosmo : astropy.cosmology.Cosmology or None
        A cosmology instance. If None, Planck15 is used.

    Returns
    -------
    astropy.cosmology.Cosmology
        The cosmology instance.
    """
    if cosmo is None:
        return Planck15
    if isinstance(cosmo, Cosmology):
        return cosmo
    raise TypeError(f"cosmo must be an astropy Cosmology instance, got {type(cosmo)}")


def luminosity_distance(redshift, cosmo=None):
    """
    Compute luminosity distance at given redshift.

    For z < 1e-6, uses the Hubble law as a minimum distance floor
    to avoid division-by-zero issues in flux computations.

    Parameters
    ----------
    redshift : float
        Redshift.
    cosmo : astropy.cosmology.Cosmology or None
        Cosmology instance. Defaults to Planck15.

    Returns
    -------
    astropy.units.Quantity
        Luminosity distance in Mpc.
    """
    c = get_cosmology(cosmo)
    dl = c.luminosity_distance(redshift).to(U.Mpc)
    # Floor at Hubble law distance for z=0 (D = cz/H0)
    if dl.value < 0.01:
        dl_hubble = (const.c * max(redshift, 1e-6) / c.H0).to(U.Mpc)
        dl = max(dl, dl_hubble)
    return dl


def angular_diameter_distance(redshift, cosmo=None):
    """
    Compute angular diameter distance at given redshift.

    For z < 1e-6, uses the Hubble law as a minimum distance floor.

    Parameters
    ----------
    redshift : float
        Redshift.
    cosmo : astropy.cosmology.Cosmology or None
        Cosmology instance. Defaults to Planck15.

    Returns
    -------
    astropy.units.Quantity
        Angular diameter distance in Mpc.
    """
    c = get_cosmology(cosmo)
    da = c.angular_diameter_distance(redshift).to(U.Mpc)
    # Floor at Hubble law distance / (1+z)^2 for z=0
    if da.value < 0.01:
        z_eff = max(redshift, 1e-6)
        dl_hubble = (const.c * z_eff / c.H0).to(U.Mpc)
        da = dl_hubble / (1 + z_eff) ** 2
    return da


def comoving_distance(redshift, cosmo=None):
    """
    Compute comoving distance at given redshift.

    Parameters
    ----------
    redshift : float
        Redshift.
    cosmo : astropy.cosmology.Cosmology or None
        Cosmology instance. Defaults to Planck15.

    Returns
    -------
    astropy.units.Quantity
        Comoving distance in Mpc.
    """
    c = get_cosmology(cosmo)
    return c.comoving_distance(redshift).to(U.Mpc)


def hubble_parameter(redshift, cosmo=None):
    """
    Compute Hubble parameter H(z) at given redshift.

    Parameters
    ----------
    redshift : float
        Redshift.
    cosmo : astropy.cosmology.Cosmology or None
        Cosmology instance. Defaults to Planck15.

    Returns
    -------
    astropy.units.Quantity
        H(z) in km/s/Mpc.
    """
    c = get_cosmology(cosmo)
    return c.H(redshift)


def systemic_velocity_from_redshift(redshift):
    """
    Compute systemic recession velocity from cosmological redshift.

    Uses the special-relativistic (kinematic) Doppler formula:
    v/c = ((1+z)^2 - 1) / ((1+z)^2 + 1)

    .. note::
        This is the *optical/kinematic* convention. It is NOT the convention used by
        the data cube's spectral axis, which follows the radio convention - see
        :func:`radio_velocity_from_redshift`.

    Parameters
    ----------
    redshift : float
        Cosmological redshift.

    Returns
    -------
    astropy.units.Quantity
        Recession velocity in km/s.
    """
    c_light = const.c.to(U.km / U.s)
    z = redshift
    v_over_c = ((1 + z) ** 2 - 1) / ((1 + z) ** 2 + 1)
    return (v_over_c * c_light).to(U.km / U.s)


def radio_velocity_from_redshift(redshift, rest_freq=None):
    """
    Compute the recession velocity of a line at ``redshift`` in the radio convention.

    The radio (or "radio definition") Doppler convention relates the observed frequency
    of a line to the velocity of its source as ``f_obs = f_rest * (1 - v/c)``, so

    .. math::
        v = c \\left(1 - \\frac{f_{\\rm obs}}{f_{\\rm rest}}\\right)
          = c \\frac{z}{1 + z}.

    This is the convention used throughout the WCS-based spectral machinery
    (:class:`~astropy.coordinates.SpectralCoord` with ``doppler_convention="radio"``,
    FITS ``VRAD`` axes and :func:`astropy.units.doppler_radio`), so a source systemic
    velocity intended to line up with a data cube's spectral axis must be computed
    with this function rather than :func:`systemic_velocity_from_redshift`.

    Unlike ``H(z) * D_C`` this never exceeds the speed of light: it tends to ``c``
    as ``z -> inf``.

    Parameters
    ----------
    redshift : float
        Cosmological redshift.

    rest_freq : astropy.units.Quantity or None
        Rest frequency of the line. Defaults to the HI 21cm line.

    Returns
    -------
    astropy.units.Quantity
        Recession velocity in km/s.
    """
    if rest_freq is None:
        rest_freq = HIfreq
    f_obs = observed_frequency(redshift, rest_freq)
    return (const.c * (1.0 - f_obs / rest_freq)).to(U.km / U.s)


def observed_frequency(redshift, rest_freq=None):
    """
    Compute observed frequency of HI line at given redshift.

    f_obs = f_rest / (1 + z)

    Parameters
    ----------
    redshift : float
        Redshift.
    rest_freq : astropy.units.Quantity or None
        Rest frequency. Defaults to HI 21cm line.

    Returns
    -------
    astropy.units.Quantity
        Observed frequency.
    """
    if rest_freq is None:
        rest_freq = HIfreq
    return rest_freq / (1 + redshift)


def rho_code_to_physical(rho_code, redshift, hubble_h):
    """
    Convert TNG comoving density code units to physical g/cm^3.

    TNG code unit for density: 1e10 M_sun / kpc^3 * h^2 (comoving).
    Physical = comoving * (1+z)^3.

    Parameters
    ----------
    rho_code : ndarray
        Density in TNG code units.
    redshift : float
        Redshift.
    hubble_h : float
        Dimensionless Hubble parameter (h = H0 / 100).

    Returns
    -------
    ndarray
        Physical density in g/cm^3.
    """
    rho_com = rho_code * 1e10 * const.M_sun * hubble_h**2 / U.kpc**3
    rho_phys = rho_com.to(U.g / U.cm**3) * (1.0 + redshift) ** 3
    return rho_phys.value


def rho_code_to_nH(rho_code, redshift, hubble_h, X_H=HYDROGEN_MASSFRAC):
    """
    Convert TNG comoving density code units to hydrogen number density nH.

    Parameters
    ----------
    rho_code : ndarray
        Density in TNG code units.
    redshift : float
        Redshift.
    hubble_h : float
        Dimensionless Hubble parameter.
    X_H : float
        Hydrogen mass fraction.

    Returns
    -------
    ndarray
        Hydrogen number density in cm^-3.
    """
    rho_phys = rho_code_to_physical(rho_code, redshift, hubble_h) * U.g / U.cm**3
    n_H = (X_H * rho_phys / const.m_p).to(1.0 / U.cm**3)
    return n_H.value


def calc_Jeans_length_physical(u_code, rho_code, redshift, hubble_h):
    """
    Compute Jeans length in physical kpc.

    Parameters
    ----------
    u_code : ndarray
        Internal energy in TNG code units [(km/s)^2].
    rho_code : ndarray
        Density in TNG code units.
    redshift : float
        Redshift.
    hubble_h : float
        Dimensionless Hubble parameter.

    Returns
    -------
    ndarray
        Jeans length in physical kpc.
    """
    u_cgs = u_code * 1e10  # (km/s)^2 -> cm^2/s^2
    rho_cgs = rho_code_to_physical(rho_code, redshift, hubble_h)
    L_J_cm = np.sqrt(
        GAMMA * GAMMA_MINUS1 * u_cgs / (const.G.cgs.value * np.maximum(rho_cgs, 1e-40))
    )
    return L_J_cm / kpc_to_cm


def calc_Sigma_neutral_physical(u_code, rho_code, f_neutral_XH, redshift, hubble_h):
    """
    Compute neutral gas surface density Sigma = rho_phys * L_J.

    Parameters
    ----------
    u_code : ndarray
        Internal energy in TNG code units.
    rho_code : ndarray
        Density in TNG code units.
    f_neutral_XH : ndarray
        f_neutral * X_H (neutral hydrogen mass fraction correction).
    redshift : float
        Redshift.
    hubble_h : float
        Dimensionless Hubble parameter.

    Returns
    -------
    ndarray
        Surface density in M_sun / kpc^2.
    """
    L_J_kpc = calc_Jeans_length_physical(u_code, rho_code, redshift, hubble_h)
    rho_msun_kpc3 = rho_code * 1e10 * hubble_h**2 * (1.0 + redshift) ** 3
    return np.maximum(L_J_kpc * rho_msun_kpc3 * f_neutral_XH, 1e-20)
