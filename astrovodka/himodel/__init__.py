"""
HI and H2 mass models for IllustrisTNG gas cells.

Four selectable methods, the two-by-two of {Gnedin & Draine 2014,
Krumholz+ 2013} x {volume, projected}, all as implemented by
Diemer+2018 (arXiv:1806.02341):

===============  ============================================================
``GD14_vol``     Appendix C.3: per-cell surface density ``Sigma = lambda_J rho``
``GD14_map``     §2.3 projection + Appendix C.3
``K13_vol``      Appendix C.4: per-cell
``K13_map``      §2.3 projection + Appendix C.4
===============  ============================================================

All four replace MARTINI's pressure-based partition, which is calibrated on
nearby galaxies and has no redshift dependence.

Quick start::

    from astrovodka.himodel import HIModelConfig, compute_hi_mass

    cfg = HIModelConfig(h=0.6774, z=0.0)
    res = compute_hi_mass(gas, model='K13_map', cfg=cfg, stars=stars, sub=sub)
    res.m_hi          # per gas cell, TNG code units (1e10 Msun/h)
    res.m_hi_msun     # total, Msun

Notes
-----
The volume methods derive a surface density from the cell's Jeans length
(Diemer+2018 Eq. 3); the projected methods build a face-on map of the galaxy
and read each cell's molecular fraction off it, so they need the stellar
particles and the subhalo record as well.

For the projected methods, cells that fall outside the projection box get zero
HI and H2: the model's galaxy *is* that box, and the catalog total for the
method is the in-box total. The mass that falls outside is reported as
``f_gas_in_map`` in `~astrovodka.himodel.galaxy.HiMassResult.diagnostics`.
Set ``HIModelConfig(map_clip_outside=True)`` to keep it instead (mass
conserving, at the cost of piling it onto the edge pixel).

The physical expressions are ported verbatim from the `TNGHIcalc` package,
which is verified bit-for-bit against Diemer+2018's reference implementation;
see `docs/en/verification.md` for that chain of evidence. Editing anything in this
subpackage means re-running the bit-exactness check.
"""

from .config import (
    BOXSIZE_CKPC_PER_H_DEFAULT,
    HUBBLE_DEFAULT,
    HIModelConfig,
    K13_VARIANTS,
    RHO_DM_MODES,
    R_SOFT_MAP_MODES,
)
from .constants import GAS_FIELDS_FULL, GAS_FIELDS_MIN, MAP_GAS_FIELDS, MAP_STAR_FIELDS
from .galaxy import (
    Cells,
    HiMassResult,
    compute_galaxy_hi_h2_map,
    compute_hi_mass,
)
from .models import (
    DEFAULT_METHOD,
    MAP_METHODS,
    METHODS,
    MODEL_SPEC,
    VOL_METHODS,
    is_known,
    resolve,
)
from .physics import calc_UMW_global_Faucher11

__all__ = [
    # configuration
    "HIModelConfig",
    "HUBBLE_DEFAULT",
    "BOXSIZE_CKPC_PER_H_DEFAULT",
    "K13_VARIANTS",
    "RHO_DM_MODES",
    "R_SOFT_MAP_MODES",
    # the four methods
    "METHODS",
    "VOL_METHODS",
    "MAP_METHODS",
    "DEFAULT_METHOD",
    "MODEL_SPEC",
    "resolve",
    "is_known",
    # entry points
    "compute_hi_mass",
    "compute_galaxy_hi_h2_map",
    "calc_UMW_global_Faucher11",
    # result containers
    "HiMassResult",
    "Cells",
    # TNG field lists
    "GAS_FIELDS_MIN",
    "GAS_FIELDS_FULL",
    "MAP_GAS_FIELDS",
    "MAP_STAR_FIELDS",
]
