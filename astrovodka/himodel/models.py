"""
Projected (map) molecular-fraction models and the registry of the four methods.

Ported line by line from the verified `TNGHIcalc` package (see `physics.py` for
the provenance and the bit-exactness contract).

The four methods are the two-by-two of {Gnedin & Draine 2014, Krumholz+ 2013}
x {volume, projected}, all as implemented by Diemer+2018:

=============  ==========================================================
``GD14_vol``   Appendix C.3: per-cell surface density ``Sigma = lambda_J rho``
``GD14_map``   §2.3 projection + Appendix C.3
``K13_vol``    Appendix C.4: per-cell
``K13_map``    §2.3 projection + Appendix C.4
=============  ==========================================================
"""

from __future__ import annotations

import numpy as np

from .config import K13_VARIANTS, HIModelConfig
from .constants import (
    BOLTZMANN,
    DM_C_A0,
    DM_C_A1,
    DM_C_B,
    G_CGS,
    K13_ALPHA_HYDRO,
    K13_C_W_KMS,
    K13_F_C,
    K13_FWT,
    K13_ITER_INIT,
    K13_ITER_MAX,
    K13_ITER_TOL,
    K13_ITER_WEIGHT,
    K13_T_CNM_MAX,
    K13_ZETA_D,
    MPC,
    MSUN,
    MSUN_PC2_TO_G_CM2,
    MSUN_PC3_TO_G_CM3,
    PC_TO_CM,
    SIGMA_STAR_FLOOR_KMS,
    Z_SUN,
)
from .physics import f_mol_GD14_vol, f_mol_K13_vol


def f_mol_GD14_map(Sigma_n, Z_map, U_MW, pix_pc):
    """GD14 projected model (Diemer+2018 Appendix C.3, last paragraph).

    Term for term identical to `f_mol_GD14_vol`, with **one** difference:

        vol:  S = L_cell/100 pc, L_cell = (m_cell/rho_cell)^(1/3), per cell
        map:  S = pixel size / 100 pc, a single scalar for the whole map

    Since the pixel size is usually >= 469 pc, ``S >= 4.7`` and hence
    ``D_* == 0.170``: the map version *is* the resolution-free limit of GD14
    (measured ``D_*`` in the volume version spans [0.170, 0.183]).
    """
    Sigma_n = np.maximum(np.asarray(Sigma_n, dtype="float64"), 1e-20)
    D_MW = np.maximum(np.asarray(Z_map, dtype="float64") / Z_SUN, 1e-10)
    U_pos = np.maximum(np.asarray(U_MW, dtype="float64"), 0.0)

    S = float(pix_pc) / 100.0
    D_star = 0.17 * (2.0 + S**5) / (1.0 + S**5)
    g = np.maximum(np.sqrt(D_MW**2 + D_star**2), 1e-10)

    alpha = 0.5 + 1.0 / (1.0 + np.sqrt(U_pos * D_MW**2 / 600.0))
    sqrt_u = np.sqrt(0.001 + 0.1 * U_pos)
    Sigma_c = np.maximum(5.0e7 * sqrt_u / (g * (1.0 + 1.69 * sqrt_u)), 1e-20)

    R_mol = np.clip((Sigma_n / Sigma_c) ** alpha, 0.0, 1e30)
    return np.clip(R_mol / (1.0 + R_mol), 0.0, 1.0)


def _k13_fmol_core(n_CNM, Sigma_n, D_MW, U_MW):
    """The chi / tau_c / s stretch of K13, line-for-line with `f_mol_K13_vol`."""
    chi = np.clip(7.2 * U_MW * (n_CNM / 10.0) ** (-1.0), 1e-10, 1e10)
    Sigma_pc2 = np.maximum(Sigma_n, 1e-20) * 1e-6  # Msun/kpc^2 -> Msun/pc^2
    tau_c = 0.066 * K13_F_C * D_MW * np.maximum(Sigma_pc2, 1e-20)
    s = np.log(1.0 + 0.6 * chi + 0.01 * chi**2) / (0.6 * np.maximum(tau_c, 1e-20))
    f = np.where(s < 2.0, 1.0 - 3.0 * s / (4.0 + s), 0.0)
    return np.clip(f, 0.0, 1.0)


def _k13_ncnm_2p(U_MW, D_MW):
    """Two-phase equilibrium CNM density (K13 eq. C14)."""
    return 23.0 * U_MW * 4.1 / (1.0 + 3.1 * D_MW**0.365)


def _k13_ncnm_hydro(Sigma_n, R_mol, rho_sd):
    """CNM density floor from Ostriker+10 hydrostatic pressure [cm^-3].

    ``P_th = pi G Sigma^2/(4 alpha)
             (1 + 2R + sqrt((1+2R)^2 + 32 zeta_d alpha f_w c_w^2 rho_sd/(pi G Sigma^2)))``
    ``n_CNM,hydro = P_th / (1.1 k_B T_CNM,max)``

    This is the patch for ``n_CNM,2p`` tending to 0 as ``U_MW -> 0``:
    hydrostatic equilibrium demands a finite pressure.
    """
    Sig = np.maximum(Sigma_n, 1e-20) * MSUN_PC2_TO_G_CM2  # g/cm^2
    rho = np.maximum(rho_sd, 0.0) * MSUN_PC3_TO_G_CM3  # g/cm^3
    cw = K13_C_W_KMS * 1e5
    base = np.maximum(np.pi * G_CGS * Sig**2, 1e-300)
    root = np.sqrt(
        np.maximum(
            (1.0 + 2.0 * R_mol) ** 2
            + 32.0 * K13_ZETA_D * K13_ALPHA_HYDRO * K13_FWT * cw**2 * rho / base,
            0.0,
        )
    )
    P_th = base / (4.0 * K13_ALPHA_HYDRO) * (1.0 + 2.0 * R_mol + root)
    return P_th / (1.1 * BOLTZMANN * K13_T_CNM_MAX)


def calc_rho_star_map(Sigma_star, Sigma_n, sigma_star):
    """Stellar mid-plane density ``rho_* = Sigma_*/(2 h_*)`` (paper Appendix C.4).

    ``h_* = sigma_*^2 / (pi G (Sigma_* + Sigma_n))``. Surface densities in
    Msun/kpc^2 and ``sigma_star`` in km/s; returns Msun/pc^3.
    """
    Sn = np.maximum(Sigma_n, 0.0) * MSUN_PC2_TO_G_CM2
    Ss = np.maximum(Sigma_star, 0.0) * MSUN_PC2_TO_G_CM2
    sig = np.maximum(np.asarray(sigma_star, dtype="float64"), SIGMA_STAR_FLOOR_KMS) * 1e5
    # The 1e-300 floor guards against dividing by zero, but it is far too small to
    # stop the division overflowing to inf for the (Sigma_* + Sigma_n) -> 0 pixels
    # that a map always has at its edges. The clamp below then makes those pixels'
    # scale height infinite and their mid-plane density zero, which is the right
    # limit -- so the overflow is harmless, and silencing it here keeps the
    # arithmetic bit-identical while sparing callers a RuntimeWarning.
    with np.errstate(over="ignore", invalid="ignore"):
        h_pc = np.maximum(
            sig**2 / np.maximum(np.pi * G_CGS * (Ss + Sn), 1e-300) / PC_TO_CM, 1.0
        )
    return np.maximum(Sigma_star, 0.0) * 1e-6 / (2.0 * h_pc)


def calc_rho_dm_nfw_map(R_kpc, m200c_msun, cfg: HIModelConfig):
    """NFW dark-matter mid-plane density [Msun/pc^3] (paper Appendix C.4).

    The concentration uses the Dutton & Maccio (2014) fit. The paper uses
    Diemer & Kravtsov (2015), which needs a sigma(M,z) solver and a
    collapse-redshift root find; that is not available in this environment, and
    the paper itself notes rho_dm is subdominant wherever molecular gas is
    significant, so the substitution is cheap.
    """
    m200 = max(float(m200c_msun), 1e-30)
    c = 10.0 ** (
        DM_C_A0 + DM_C_A1 * float(cfg.z) + DM_C_B * np.log10(m200 * cfg.h / 1e12)
    )
    # R_200c from M_200c = 200 rho_c (4 pi/3) R^3, rho_c(z) = 3H(z)^2/(8 pi G)
    hz = cfg.h * np.sqrt(0.3089 * (1.0 + cfg.z) ** 3 + 1.0 - 0.3089)
    rho_c = 3.0 * (hz * 100.0 * 1e5 / MPC) ** 2 / (8.0 * np.pi * G_CGS)
    r200_cm = (3.0 * m200 * MSUN / (4.0 * np.pi * 200.0 * rho_c)) ** (1.0 / 3.0)
    rs_cm = r200_cm / c
    rs_pc = rs_cm / PC_TO_CM
    rho_s = m200 * MSUN / (4.0 * np.pi * rs_cm**3 * (np.log(1.0 + c) - c / (1.0 + c)))
    rho_s_pc3 = rho_s * PC_TO_CM**3 / MSUN
    x = np.maximum(np.asarray(R_kpc, dtype="float64") * 1e3 / rs_pc, 1e-6)
    return rho_s_pc3 / (x * (1.0 + x) ** 2)


def f_mol_K13_map(
    Sigma_n,
    Z_map,
    U_MW,
    Sigma_star=None,
    sigma_star=None,
    rho_dm=None,
    variant=None,
    *,
    cfg: HIModelConfig,
):
    """K13 projected model (Diemer+2018 Appendix C.4).

    Note that K13 has **no** resolution term (only GD14 carries ``D_*``), so
    this function takes no pixel-size argument. That is a real difference
    between the two models in projection, not an omission.

    ``variant`` controls how ``n_CNM`` is obtained, which is the one
    unresolved difference between the paper and the existing code:

    ``'hydro'``
        ``n_CNM = max(n_CNM,2p, P_th/(1.1 k_B 243 K))`` with
        ``rho_sd = rho_* + rho_dm`` -- the paper's Appendix C.4 prescription.
        Needs iteration because ``P_th`` depends on ``R_mol``.
    ``'hydro_nostar'``
        As above but with ``rho_sd = 0``, to isolate the effect of the
        mid-plane density model.
    ``'2p'``
        Only ``n_CNM,2p``. The volume version uses
        ``max(n_H f_neutral, n_CNM,2p)``; a projected model has no volume
        density, so it degenerates to the 2p term alone.
    """
    variant = cfg.k13_variant if variant is None else variant
    if variant not in K13_VARIANTS:
        raise ValueError(f"unknown K13 variant: {variant!r}")

    Sigma_n = np.maximum(np.asarray(Sigma_n, dtype="float64"), 0.0)
    D_MW = np.maximum(np.asarray(Z_map, dtype="float64") / Z_SUN, 1e-10)
    U_MW = np.maximum(np.asarray(U_MW, dtype="float64"), 0.0)
    n_2p = _k13_ncnm_2p(U_MW, D_MW)

    if variant == "2p":
        return _k13_fmol_core(n_2p, Sigma_n, D_MW, U_MW)

    rho_sd = np.zeros_like(Sigma_n)
    if variant == "hydro":
        if rho_dm is not None:
            rho_sd = rho_sd + np.maximum(np.asarray(rho_dm, dtype="float64"), 0.0)
        if Sigma_star is not None and sigma_star is not None:
            rho_sd = rho_sd + calc_rho_star_map(Sigma_star, Sigma_n, sigma_star)

    active = Sigma_n > 0
    f = np.where(active, K13_ITER_INIT, 0.0)
    for _ in range(K13_ITER_MAX):
        R_mol = f / np.maximum(1.0 - f, 1e-6)
        n_cnm = np.maximum(n_2p, _k13_ncnm_hydro(Sigma_n, R_mol, rho_sd))
        f_new = _k13_fmol_core(n_cnm, Sigma_n, D_MW, U_MW)
        f_next = np.where(active, K13_ITER_WEIGHT * f_new + (1.0 - K13_ITER_WEIGHT) * f, 0.0)
        # Pixels with Sigma_n = 0 must not take part in the convergence test,
        # or the max is pinned by them forever
        delta = float(np.max(np.abs(f_next - f)[active])) if np.any(active) else 0.0
        f = f_next
        if delta < K13_ITER_TOL:
            break
    return np.clip(f, 0.0, 1.0)


# ============================================================
# Registry of the four selectable methods
# ============================================================
# ``kind`` picks the driver: 'vol' works per gas cell, 'map' projects first and
# then reads each cell's f_mol off the map. The three extra entries are the
# per-method knobs that the projected path needs and the volume path ignores.
MODEL_SPEC = {
    "GD14_vol": {
        "model": "GD14",
        "kind": "vol",
        "fmol": f_mol_GD14_vol,
        "label": "GD14_vol",
        "description": "Gnedin & Draine (2014) volumetric, per-cell Jeans surface density",
        "reference": "Diemer+2018 (arXiv:1806.02341) Appendix C.3",
        "needs_stars": False,
    },
    "GD14_map": {
        "model": "GD14",
        "kind": "map",
        "fmol": f_mol_GD14_map,
        "label": "GD14_map",
        "description": "GD14 projected; D_* == 0.17 resolution-free limit",
        "reference": "Diemer+2018 (arXiv:1806.02341) §2.3 + Appendix C.3",
        "needs_stars": True,
    },
    "K13_vol": {
        "model": "K13",
        "kind": "vol",
        "fmol": f_mol_K13_vol,
        "label": "K13_vol",
        "description": "Krumholz+ (2013) volumetric, per-cell Jeans surface density",
        "reference": "Diemer+2018 (arXiv:1806.02341) Appendix C.4",
        "needs_stars": False,
    },
    "K13_map": {
        "model": "K13",
        "kind": "map",
        "fmol": f_mol_K13_map,
        "label": "K13_map",
        "description": "K13 projected; n_CNM variant set by cfg.k13_variant",
        "reference": "Diemer+2018 (arXiv:1806.02341) §2.3 + Appendix C.4",
        "needs_stars": True,
    },
}

METHODS = tuple(MODEL_SPEC)
VOL_METHODS = tuple(m for m, s in MODEL_SPEC.items() if s["kind"] == "vol")
MAP_METHODS = tuple(m for m, s in MODEL_SPEC.items() if s["kind"] == "map")

#: The method used when the caller does not name one. `K13_vol` is the direct
#: replacement for the single model astrovodka shipped before, so it keeps
#: existing scripts working while giving them the corrected FFT UV field.
DEFAULT_METHOD = "K13_vol"


def resolve(method: str) -> dict:
    """Look up a method by name, with a helpful error for a wrong one.

    Also accepts the bare model names ``'GD14'``/``'K13'`` (meaning the volume
    variant) and the ``('GD14', 'map')`` two-argument spelling used by the
    production scripts.
    """
    if isinstance(method, (tuple, list)):
        if len(method) != 2:
            raise ValueError(f"method tuple must be (model, kind), got {method!r}")
        method = f"{method[0]}_{method[1]}"
    if not isinstance(method, str):
        raise TypeError(f"method must be a string, got {type(method).__name__}")

    key = method.strip()
    if key in MODEL_SPEC:
        return MODEL_SPEC[key]

    upper = key.upper()
    for name, spec in MODEL_SPEC.items():
        if name.upper() == upper:
            return spec
    # Bare model name -> the volume variant
    aliases = {spec["model"]: name for name, spec in MODEL_SPEC.items() if spec["kind"] == "vol"}
    if upper in aliases:
        return MODEL_SPEC[aliases[upper]]

    raise ValueError(
        f"unknown HI model {method!r}; choose one of {list(METHODS)} "
        f"(or the bare model name {sorted(aliases)})"
    )


def is_known(method) -> bool:
    """Whether ``method`` names one of the four models."""
    try:
        resolve(method)
    except (ValueError, TypeError):
        return False
    return True


__all__ = [
    "MODEL_SPEC",
    "METHODS",
    "VOL_METHODS",
    "MAP_METHODS",
    "DEFAULT_METHOD",
    "resolve",
    "is_known",
    "f_mol_GD14_map",
    "f_mol_K13_map",
    "calc_rho_star_map",
    "calc_rho_dm_nfw_map",
]
