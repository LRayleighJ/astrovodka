"""
Physical constants for the astrovodka HI/H2 models.

Pure constants only — **this module contains no h**. The Hubble parameter is a
per-snapshot quantity and lives in `HIModelConfig` (see `config.py`); a module
level `h` is the single easiest way to silently mix two snapshots' units.

The numerical values and the expressions that consume them are ported verbatim
from the verified `TNGHIcalc` package (itself a line-by-line restatement of
Diemer+2018's reference implementation `TNGHI_GD14model/gd14_model.py`). Every
expression here is checked bit-for-bit against that reference — see
`tests/test_hi_models.py` and `docs/en/verification.md`. When
editing, keep the arithmetic *literally* as written, including operand order.
"""

import astropy.constants as const
import astropy.units as u

PROTONMASS = const.m_p.cgs.value  # g
MSUN = const.M_sun.cgs.value  # g
MPC = (1.0 * u.Mpc).to(u.cm).value  # cm
G_CGS = const.G.cgs.value  # cm^3/g/s^2
BOLTZMANN = const.k_B.cgs.value  # erg/K

GAMMA = 5.0 / 3.0
GAMMA_MINUS1 = GAMMA - 1.0

HYDROGEN_MASSFRAC = 0.76  # solar-abundance hydrogen mass fraction (X_H fallback)
Z_SUN = 0.0127  # solar metallicity (mass fraction)
kpc_to_cm = 3.08568e21
PC_TO_CM = kpc_to_cm / 1e3

# --- FUV normalisation (Diemer+2018 §2.4 / Appendix A.2) ---------------------
Q_UV_CAL = 3.3e-6  # photons/s/cm^2/Hz per Msun/yr, Starburst99 spectrum at 1 kpc
U_MW0_CAL = 3.43e-8  # Draine (1978) field flux at 1000 A -> U_MW = 1
UVB_Z0 = 0.0014245  # Faucher-Giguere+2011 U_MW at z = 0 (documentation/selftest)

# --- per-galaxy FFT UV field defaults (matched to the K13 production chain) --
NGRID = 128  # 3D grid edge length
F_ESC = 0.1  # UV escape fraction, Diemer+2018 §4.2 calibration
BOX_EXPAND = 1.1  # adaptive box = BOX_EXPAND x coordinate extent
BOX_MIN_KPC = 10.0  # adaptive box minimum edge length
R_SOFT_FACTOR = 0.25  # 1/r^2 kernel softening = R_SOFT_FACTOR x dx
# (Appendix A.2 says "assume a distance of ~(r_cell/2)*sqrt(3)"; the text
#  extraction is ambiguous, so the K13 chain's dx/4 is kept. Measured: it puts
#  GD14 within 1.01 of the official catalog.)

# --- projection (map) defaults, Diemer+2018 §2.3 + Appendix C ---------------
# The official `hih2_galaxy_0XX.hdf5` config group (attrs, not datasets) records
# the parameters Diemer used, and these defaults match it item for item:
#   map_r_npix=128, map_range_rgas=2.0, map_range_rstr=6.0, map_range_min=30.0,
#   map_smoothing_factor=1.0, uv_ngrid=128, uv_escape_frac=0.1
NGRID_MAP = 128  # projection grid edge length (= map_r_npix)
MAP_SUB = 4  # supersampling: deposit on (sub*n)^2, then block-sum down to n^2
MAP_SIGMA_BINS = 12  # log-spaced bins for the per-cell Gaussian kernel width
R_SOFT_FACTOR_MAP = 0.25  # 2D 1/r^2 softening = factor x dx_pix (same ratio as 3D)
MAP_FACTOR_RGAS = 2.0  # L_map = 2 x max(FACTOR_RGAS*R_half_gas, FACTOR_RSTR*R_half_*, MIN)
MAP_FACTOR_RSTR = 6.0
MAP_MIN_KPC = 30.0  # absolute floor on L_map [kpc] (= map_range_min)
MIN_ALIGN_GAS = 50  # minimum gas cells to measure angular momentum (paper §2.3)
MIN_ALIGN_STAR = 50  # minimum star particles for the fallback branch
ALIGN_RGAS_FACTOR = 0.5  # gas L radius = factor x R_half_gas ("within half ...")
ALIGN_RSTR_FACTOR = 2.0  # stellar fallback radius = factor x R_half_star
MAP_ALIGN_RTOL = 0.5  # warn when |r| exceeds BoxSize/2 (periodic wrap guard)
STAR_SIGMA_FACTOR = 0.5  # star deposition kernel width = factor x dx_pix
SIGMA_STAR_FLOOR_KMS = 1.0  # stellar velocity dispersion floor (paper App. C.1)

# --- K13 hydrostatic-pressure term (paper Appendix C.4 / Ostriker+10) -------
K13_ALPHA_HYDRO = 5.0  # alpha: folds turbulent/magnetic/CR support into P_th
K13_ZETA_D = 0.33  # zeta_d: geometric factor
K13_C_W_KMS = 8.0  # c_w: warm neutral medium sound speed [km/s]
K13_FWT = 0.5  # f_w: thermal velocity dispersion / c_w
K13_T_CNM_MAX = 243.0  # T_CNM,max [K] (Wolfire+03)
K13_F_C = 5.0  # clumping factor at kpc scales
MSUN_PC2_TO_G_CM2 = 2.089e-4  # 1 Msun/pc^2 -> g/cm^2
MSUN_PC3_TO_G_CM3 = 6.77e-23  # 1 Msun/pc^3 -> g/cm^3
K13_ITER_INIT = 0.5  # initial f_mol for the hydrostatic iteration
K13_ITER_WEIGHT = 0.7  # weight on the new solution (paper: "weight of 0.7")
K13_ITER_MAX = 20
K13_ITER_TOL = 1e-4

# Dutton & Maccio (2014) concentration fit. The paper uses Diemer & Kravtsov
# (2015), which needs a sigma(M,z) solver and a collapse-redshift root find
# (colossus is not available here). The paper itself notes rho_dm is
# subdominant wherever molecular gas matters, so this substitution is cheap.
DM_C_A0, DM_C_A1 = -0.101, 0.026
DM_C_B = -0.109

K13_RHO_DM_CONST = 0.01  # [Msun/pc^3], for rho_dm_mode='const'

# TNG snapshot fields each method needs.
#
# vol: volume density + internal energy (the Jeans-length surface-density step).
# map: no InternalEnergy -- the surface density is measured directly -- but it
#      does need Velocities for the angular-momentum alignment.
GAS_FIELDS_MIN = [
    "Masses",
    "Density",
    "InternalEnergy",
    "StarFormationRate",
    "NeutralHydrogenAbundance",
    "GFM_Metals",
    "GFM_Metallicity",
    "Coordinates",
]
GAS_FIELDS_FULL = GAS_FIELDS_MIN + ["Velocities"]
MAP_GAS_FIELDS = [
    "Masses",
    "Density",
    "StarFormationRate",
    "NeutralHydrogenAbundance",
    "GFM_Metals",
    "GFM_Metallicity",
    "Coordinates",
    "Velocities",
]
# GFM_StellarFormationTime is used to drop wind particles (negative formation time)
MAP_STAR_FIELDS = [
    "Coordinates",
    "Velocities",
    "Masses",
    "GFM_StellarFormationTime",
]
