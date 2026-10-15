"""
Unit conversions and the volume (per-cell) molecular-fraction models.

Ported line by line from the verified `TNGHIcalc` package, which is itself a
verbatim restatement of Diemer+2018's reference implementation
``TNGHI_GD14model/gd14_model.py``. **The physical expressions must stay
bit-for-bit identical to that reference**, operand order included; the only
change made here is that the runtime parameters arrive as an explicit
`~astrovodka.himodel.config.HIModelConfig` instead of through a global
accessor. See the module docstring of `config.py` for why.
"""

from __future__ import annotations

import astropy.constants as const
import astropy.units as u
import numpy as np

from .config import HIModelConfig
from .constants import (
    G_CGS,
    GAMMA,
    GAMMA_MINUS1,
    HYDROGEN_MASSFRAC,
    Z_SUN,
    kpc_to_cm,
)


def rho_code_to_physical(rho_code, cfg: HIModelConfig):
    """TNG comoving density code unit -> physical g/cm^3 at the snapshot z.

    The conversion runs through astropy units rather than a hand-rolled
    constant on purpose: the reference implementation does, and astropy's kpc
    (3.0856775814913673e21 cm) differs from this module's ``kpc_to_cm``
    (3.08568e21 cm) by ~8e-7 relative. Substituting one for the other would be
    a silent 1e-6-level change to every density, so keep the quantity chain.
    """
    rho_com = rho_code * 1e10 * const.M_sun * cfg.h**2 / u.kpc**3
    return rho_com.to(u.g / u.cm**3).value * (1.0 + cfg.z) ** 3


def rho_code_to_msun_kpc3(rho_code, cfg: HIModelConfig):
    """TNG comoving density code unit -> physical Msun/kpc^3 at the snapshot z."""
    return rho_code * 1e10 * cfg.h**2 * (1.0 + cfg.z) ** 3


def rho_code_to_nH(rho_code, X_H=HYDROGEN_MASSFRAC, *, cfg: HIModelConfig):
    """TNG comoving density code unit -> hydrogen number density [cm^-3].

    ``X_H`` is the hydrogen mass fraction; for TNG gas cells this is
    ``GFM_Metals[:, 0]``.
    """
    rho_phys = rho_code_to_physical(rho_code, cfg) * u.g / u.cm**3
    return (X_H * rho_phys / const.m_p).to(1.0 / u.cm**3).value


def calc_Jeans_length_physical(u_code, rho_code, cfg: HIModelConfig):
    """Jeans length [physical kpc]: ``lambda_J = sqrt(gamma(gamma-1)u / (G rho))``.

    Diemer+2018 Eq. (3). ``u`` is the specific internal energy (code unit =
    1e10 (km/s)^2); ``rho`` is the physical density.
    """
    u_cgs = u_code * 1e10
    rho_cgs = rho_code_to_physical(rho_code, cfg)
    L_J_cm = np.sqrt(GAMMA * GAMMA_MINUS1 * u_cgs / (G_CGS * rho_cgs))
    return L_J_cm / kpc_to_cm


def calc_Sigma_neutral_physical(u_code, rho_code, f_neutral_XH, cfg: HIModelConfig):
    """Neutral surface density ``Sigma_HI+H2 = lambda_J * rho_HI+H2`` [Msun/kpc^2].

    Diemer+2018 §2.3: the volume version converts a volume density into a
    surface density through the Jeans length, ``Sigma_X = lambda_J rho_X``.
    Pass ``f_neutral * X_H``, i.e. the hydrogen (HI+H2) mass fraction.

    Note that the paper's Appendix C.3/C.4 ``Sigma_HI+H2`` is the **hydrogen
    nucleus** surface density (no He), with no extra factor of 1/X_H.
    """
    L_J_kpc = calc_Jeans_length_physical(u_code, rho_code, cfg)
    rho_msun_kpc3 = rho_code_to_msun_kpc3(rho_code, cfg)
    return np.maximum(L_J_kpc * rho_msun_kpc3 * f_neutral_XH, 1e-20)


def calc_Lcell_pc(mgas_code, rho_code, cfg: HIModelConfig):
    """GD14 resolution scale ``L_cell = (m_cell / rho_cell)^(1/3)`` [pc].

    Diemer+2018 Appendix C.3: moving-mesh cells have no well-defined edge
    length, so take ``L_cell = V_cell^(1/3) = (m_cell/rho_cell)^(1/3)``; the
    projected version uses the pixel size instead. This is **not** the Jeans
    length.
    """
    rho_msun_kpc3 = rho_code_to_msun_kpc3(rho_code, cfg)
    m_gas_msun = mgas_code * 1e10 / cfg.h
    L_cell_kpc = np.maximum(
        (m_gas_msun / np.maximum(rho_msun_kpc3, 1e-20)) ** (1.0 / 3.0), 1e-10
    )
    return L_cell_kpc * 1e3


def calc_UMW_global_Faucher11(redshift):
    """Faucher-Giguere et al. (2011) global UV background, in Draine units."""
    log10_UMW = (
        -0.004292 * redshift**4
        + 0.061446 * redshift**3
        - 0.358413 * redshift**2
        + 1.204787 * redshift
        - 2.846354
    )
    return float(10.0**log10_UMW)


# ============================================================
# Neutral fraction
# ============================================================
def get_f_neutral_H(nH_code, sfr_code):
    """Estimate ``f_neutral = f_HI + f_H2``, the neutral hydrogen mass fraction.

    TNG's ``NeutralHydrogenAbundance`` is **f_HI** (already including the
    Rahmati+2013 self-shielding correction), not ``f_neutral``. Star-forming
    cells (SFR > 0) are cold phase and essentially fully neutral, so take 1;
    for non-star-forming cells H2 is negligible and ``f_neutral ~ f_HI``
    (~10% bias).
    """
    f_neutral = np.array(nH_code, dtype="float64", copy=True)
    f_neutral[sfr_code > 0.0] = 1.0
    return f_neutral


def f_mol_GD14_vol(
    rho_code, u_code, Z, nH_code, sfr_code, X_H, U_MW, cfg: HIModelConfig, mgas=None
):
    """GD14 (Gnedin & Draine 2014) molecular fraction, Diemer+2018 Appendix C.3.

        R_mol = (Sigma_HI+H2 / Sigma_c)^alpha ,   f_mol = R_mol / (1 + R_mol)
        alpha   = 0.5 + 1 / (1 + sqrt(U_MW * D_MW^2 / 600))
        Sigma_c = 5e7 * sqrt(0.001 + 0.1 U_MW)
                  / ( g * (1 + 1.69 sqrt(0.001 + 0.1 U_MW)) )        [Msun/kpc^2]
        g       = sqrt(D_MW^2 + D_*^2)
        D_*     = 0.17 * (2 + S^5) / (1 + S^5),   S = L_cell / 100 pc
        L_cell  = (m_cell / rho_phys)^(1/3)

    Returns ``(f_mol, f_HI, f_H2)``, all fractions of the **hydrogen mass**:
    ``f_HI = f_neutral * (1 - f_mol)``, ``f_H2 = f_neutral * f_mol``. The caller
    then forms ``m_HI = m_gas * X_H * f_HI``.
    """
    f_neutral = get_f_neutral_H(nH_code, sfr_code)  # dimensionless, [0, 1]
    D_MW = np.maximum(Z / Z_SUN, 1e-10)  # dust-to-gas ratio, solar units

    # Neutral hydrogen surface density [Msun/kpc^2] (same Jeans approximation as K13 vol)
    Sigma_neutral = calc_Sigma_neutral_physical(
        u_code, rho_code, f_neutral * X_H, cfg
    )
    Sigma_neutral = np.maximum(Sigma_neutral, 1e-20)

    # Resolution parameter S = L_cell/100 pc (required by GD14; Jeans length is the fallback)
    if mgas is not None:
        L_cell_pc = calc_Lcell_pc(mgas, rho_code, cfg)
    else:
        L_cell_pc = np.maximum(
            calc_Jeans_length_physical(u_code, rho_code, cfg) * 1e3, 1.0
        )
    S = L_cell_pc / 100.0
    D_star = 0.17 * (2.0 + S**5) / (1.0 + S**5)  # (0.17, 0.34]

    U_pos = np.maximum(U_MW, 0.0)
    g = np.maximum(np.sqrt(D_MW**2 + D_star**2), 1e-10)
    alpha = 0.5 + 1.0 / (1.0 + np.sqrt(U_pos * D_MW**2 / 600.0))  # (0.5, 1.5]

    sqrt_u = np.sqrt(0.001 + 0.1 * U_pos)
    Sigma_c = np.maximum(5.0e7 * sqrt_u / (g * (1.0 + 1.69 * sqrt_u)), 1e-20)

    R_mol = np.clip((Sigma_neutral / Sigma_c) ** alpha, 0.0, 1e30)
    f_mol = np.clip(R_mol / (1.0 + R_mol), 0.0, 1.0)

    return f_mol, f_neutral * (1.0 - f_mol), f_neutral * f_mol


def f_mol_K13_vol(
    rho_code, u_code, Z, nH_code, sfr_code, X_H, U_MW, cfg: HIModelConfig, mgas=None
):
    """K13 (Krumholz, McKee & Tumlinson 2013) molecular fraction, Diemer+2018 App. C.4.

        n_CNM    = max(n_H f_neutral, n_CNM,2p)
        n_CNM,2p = 23 U_MW * 4.1 / (1 + 3.1 D_MW^0.365)              [cm^-3]
        chi      = 7.2 U_MW (n_CNM/10 cm^-3)^-1
        tau_c    = 0.066 f_c D_MW Sigma_HI+H2 [Msun/pc^2],  f_c = 5
        s        = ln(1 + 0.6 chi + 0.01 chi^2) / (0.6 tau_c)
        f_mol    = 1 - 3s/(4+s)  (s < 2),  0 (s >= 2)

    Returns ``(f_mol, f_HI, f_H2)`` as fractions of the hydrogen mass.
    """
    f_neutral = get_f_neutral_H(nH_code, sfr_code)
    D_MW = np.maximum(Z / Z_SUN, 1e-10)

    n_H = np.maximum(rho_code_to_nH(rho_code, X_H, cfg=cfg), 1e-10)
    n_CNM_gas = np.maximum(n_H * f_neutral, 1e-10)
    n_CNM_2p = 23.0 * U_MW * 4.1 / (1.0 + 3.1 * D_MW**0.365)  # K13 Eq. C14
    n_CNM = np.maximum(n_CNM_gas, n_CNM_2p)

    chi = np.clip(7.2 * U_MW * (n_CNM / 10.0) ** (-1.0), 1e-10, 1e10)

    Sigma_neutral = calc_Sigma_neutral_physical(
        u_code, rho_code, f_neutral * X_H, cfg
    )
    Sigma_Msun_pc2 = Sigma_neutral * 1e-6  # kpc^-2 -> pc^-2
    f_c = 5.0
    tau_c = 0.066 * f_c * D_MW * np.maximum(Sigma_Msun_pc2, 1e-20)

    s = np.log(1.0 + 0.6 * chi + 0.01 * chi**2) / (0.6 * np.maximum(tau_c, 1e-20))
    f_mol = np.where(s < 2.0, 1.0 - 3.0 * s / (4.0 + s), 0.0)
    f_mol = np.clip(f_mol, 0.0, 1.0)

    return f_mol, f_neutral * (1.0 - f_mol), f_neutral * f_mol


def hydrogen_fraction(gfm_metals):
    """Hydrogen mass fraction from a TNG ``GFM_Metals`` array.

    ``loadSubset(..., mdi=0)`` — the way every astrovodka TNG path loads this
    field — returns a 1-D array holding only the hydrogen column. A full-width
    2-D array (element index 0 = H) is equally valid, so accept both rather
    than making the caller remember which one it has.
    """
    gfm = np.asarray(gfm_metals, dtype="float64")
    if gfm.ndim == 2:
        return gfm[:, 0]
    if gfm.ndim == 1:
        return gfm
    raise ValueError(f"GFM_Metals must be 1-D or 2-D, got {gfm.ndim}-D")


def _gas_columns(gas):
    """Pull the fields every volume model needs out of a TNG gas dict."""
    rho_code = gas["Density"][:].astype("float64")
    u_code = gas["InternalEnergy"][:].astype("float64")
    nH_code = gas["NeutralHydrogenAbundance"][:].astype("float64")
    sfr_code = gas["StarFormationRate"][:].astype("float64")
    X_H = hydrogen_fraction(gas["GFM_Metals"])
    Z = gas["GFM_Metallicity"][:].astype("float64")
    m_gas = gas["Masses"][:].astype("float64")
    return rho_code, u_code, nH_code, sfr_code, X_H, Z, m_gas


__all__ = [
    "rho_code_to_physical",
    "rho_code_to_msun_kpc3",
    "rho_code_to_nH",
    "calc_Jeans_length_physical",
    "calc_Sigma_neutral_physical",
    "calc_Lcell_pc",
    "calc_UMW_global_Faucher11",
    "get_f_neutral_H",
    "f_mol_GD14_vol",
    "f_mol_K13_vol",
]
