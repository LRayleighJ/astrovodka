"""
Per-galaxy drivers: the single entry point that dispatches to all four methods.

Ported line by line from the verified `TNGHIcalc` package (see `physics.py` for
the provenance and the bit-exactness contract).

The volume methods work per gas cell: convert the cell's density and internal
energy into a Jeans-length surface density, convert that into ``f_mol``, and
split the cell's hydrogen mass. The projected methods build face-on maps first,
read ``f_mol`` off the map for each cell, and split the cell's neutral mass by
that value -- which is what makes the projected models usable for a data cube
rather than only for a catalog total.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import map_coordinates

from .config import HIModelConfig
from .constants import (
    K13_RHO_DM_CONST,
    MAP_FACTOR_RGAS,
    MAP_FACTOR_RSTR,
    MAP_MIN_KPC,
    SIGMA_STAR_FLOOR_KMS,
    Z_SUN,
)
from .models import DEFAULT_METHOD, resolve
from .physics import calc_Lcell_pc, get_f_neutral_H, hydrogen_fraction
from .uvfield import (
    _map_grid,
    _rotation_to_z,
    compute_align_vector,
    compute_galaxy_fft_uv,
    compute_galaxy_fft_uv_map,
    deposit_gaussian_maps,
)


@dataclass
class Cells:
    """Per-gas-cell HI/H2 masses, in TNG code units (1e10 Msun/h).

    Every mass and fraction array has the **same length as the caller's gas
    arrays**, for both families of models. Where a model excludes a cell -- the
    projected models drop those outside the projection box -- the arrays get
    zero rather than the cell being omitted, so a result can be handed straight
    to `~astrovodka.sources.SPHSource` without an index dance.

    ``row_index`` is the list of cells that actually contributed, for callers
    that want to know which ones those were; it is *not* needed to interpret
    the arrays.
    """

    row_index: np.ndarray
    f_mol: np.ndarray
    f_neutral: np.ndarray
    m_hi: np.ndarray
    m_h2: np.ndarray
    u_mw: np.ndarray
    meta: dict = field(default_factory=dict)

    @property
    def m_neutral(self):
        # m_hi + m_h2 equals the neutral hydrogen mass **exactly** under both
        # families: the volume models give f_hi + f_h2 = f_neutral, the
        # projected ones give (1 - f) + f = 1.
        return self.m_hi + self.m_h2


def compute_HI_mass_particle(gas, U_MW, model="GD14", cfg: HIModelConfig | None = None):
    """Per-cell ``(m_HI, m_H2)`` in the same units as ``Masses`` (1e10 Msun/h).

    Thin wrapper around a volume model's ``f_mol``; only meaningful for the
    ``*_vol`` methods, which take the UV field as an argument rather than
    deriving it from a map.
    """
    from .physics import _gas_columns

    rho_code, u_code, nH_code, sfr_code, X_H, Z, m_gas = _gas_columns(gas)
    f_mol_func = resolve(model)["fmol"]
    _, f_HI, f_H2 = f_mol_func(
        rho_code, u_code, Z, nH_code, sfr_code, X_H, U_MW, cfg, mgas=m_gas
    )
    return m_gas * X_H * f_HI, m_gas * X_H * f_H2


def compute_galaxy_hi_h2_map(
    gas,
    stars,
    sub,
    cfg: HIModelConfig,
    model="GD14",
    k13_variant=None,
    rho_dm_mode=None,
    z_weight="gas",
    star_sigma="half-pixel",
    align=True,
    boxsize=None,
    _cells_out=None,
):
    """Projected HI/H2 masses for one subhalo (the only entry the map path needs).

    Parameters
    ----------
    gas, stars : dict
        TNG particle data for the subhalo.
    sub : dict
        ``pos`` (3,) in raw ckpc/h, ``vel`` (3,) in km/s, ``rhalf_gas`` and
        ``rhalf_star`` in raw ckpc/h, and optionally ``m200c_msun``.
    cfg : HIModelConfig
        Runtime configuration.
    model : {'GD14', 'K13'}
        Which projected model to use.
    boxsize : float, optional
        Periodic box length in physical kpc, enabling a wrap guard for
        subhalos near a box edge. Off by default (matching the reference).

    Returns
    -------
    dict
        ``m_hi``, ``m_h2``, ``m_neutral`` in TNG code units (1e10 Msun/h), the
        pre-rotation alignment vector ``align_ang_mom_vector``, the alignment
        diagnostics, ``f_gas_in_map`` (the geometric truncation diagnostic) and
        ``pix_kpc`` / ``L_map_kpc``.
    """
    model = model.upper()
    rho_dm_mode = cfg.rho_dm_mode if rho_dm_mode is None else rho_dm_mode

    ngas = int(gas.get("count", 0)) if gas is not None else 0
    nstar = int(stars.get("count", 0)) if stars is not None else 0

    # ---- 1. face-on alignment ---------------------------------------------
    if align and ngas + nstar > 0:
        vhat, ok, ainfo = compute_align_vector(
            gas, stars, sub["pos"], sub["vel"], sub["rhalf_gas"], sub["rhalf_star"]
        )
    else:
        vhat, ok = np.array([0.0, 0.0, 1.0]), False
        ainfo = {"n_gas_align": ngas, "n_star_align": nstar, "branch": "none"}
    R = _rotation_to_z(vhat)

    # ---- 2. projection geometry (physical kpc here) ------------------------
    conv = cfg.h * (1.0 + cfg.z)  # raw ckpc/h -> physical kpc
    rgas_kpc = sub["rhalf_gas"] / conv
    rstar_kpc = sub["rhalf_star"] / conv
    L_map = 2.0 * max(
        MAP_FACTOR_RGAS * rgas_kpc, MAP_FACTOR_RSTR * rstar_kpc, MAP_MIN_KPC
    )
    dx = L_map / cfg.ngrid_map
    pix_pc = dx * 1e3
    n = cfg.ngrid_map

    out = {
        "align_ang_mom_vector": vhat,
        "align_success": int(ok),
        "n_gas_align": ainfo["n_gas_align"],
        "n_star_align": ainfo["n_star_align"],
        "align_branch": ainfo["branch"],
        "pix_kpc": dx,
        "L_map_kpc": L_map,
    }

    if ngas == 0 and nstar == 0:
        out.update({"m_hi": 0.0, "m_h2": 0.0, "m_neutral": 0.0, "f_gas_in_map": 0.0})
        return out

    m_gas_all = np.asarray(gas["Masses"], dtype="float64") if ngas else np.zeros(0)

    # ---- 3. rotate face-on and take relative coordinates -------------------
    def _xy(part, cnt):
        if cnt == 0:
            return np.zeros((0, 2))
        r = (np.asarray(part["Coordinates"], dtype="float64") - sub["pos"]) / conv
        if boxsize is not None:
            r = (r + boxsize / 2.0) % boxsize - boxsize / 2.0  # wrap guard
        return (r @ R.T)[:, :2]

    xy_g = _xy(gas, ngas)
    xy_s = _xy(stars, nstar)

    # ---- 4. per-point weights ---------------------------------------------
    maps = {}
    sigma_eff = 0.0
    if ngas:
        X_H = hydrogen_fraction(gas["GFM_Metals"])
        Zc = np.asarray(gas["GFM_Metallicity"], dtype="float64")
        nHc = np.asarray(gas["NeutralHydrogenAbundance"], dtype="float64")
        sfrc = np.asarray(gas["StarFormationRate"], dtype="float64")
        f_neu = get_f_neutral_H(nHc, sfrc)
        rho_c = np.asarray(gas["Density"], dtype="float64")
        # The projected models need no InternalEnergy -- there is no Jeans step
        m_n = m_gas_all * X_H * f_neu  # neutral hydrogen mass
        w = {
            "Sigma_n": m_n,
            "Sigma_gas": m_gas_all,
            "Sigma_sfr": sfrc,
            "MZ_gas": m_gas_all * Zc,
            "MZ_n": m_n * Zc,
        }
        sigma_cell = 0.5 * calc_Lcell_pc(m_gas_all, rho_c, cfg) / 1e3  # kpc
        # The neutral-mass-weighted median deposition width = the map's own
        # resolution scale in kpc, used by r_soft_map_mode = 'sigma'/'cell'
        if np.sum(m_n) > 0:
            _o = np.argsort(sigma_cell)
            _cw = np.cumsum(m_n[_o]) / np.sum(m_n)
            sigma_eff = float(sigma_cell[_o][np.searchsorted(_cw, 0.5)])
        maps.update(deposit_gaussian_maps(xy_g, w, sigma_cell, L_map, cfg, n))
        if _cells_out is not None:  # per-cell extension (for cubes / particle files)
            _cells_out["X_H"] = X_H
            _cells_out["f_neu"] = f_neu
            _cells_out["m_n"] = m_n

    if nstar:
        sft = np.asarray(stars["GFM_StellarFormationTime"], dtype="float64")
        keep = sft > 0.0
        m_s = np.asarray(stars["Masses"], dtype="float64")
        vz = np.asarray(stars["Velocities"], dtype="float64")[:, 2]
        if star_sigma == "half-pixel":
            sig_s = np.full(nstar, 0.5 * dx)
        elif star_sigma == "gas":
            sig_s = (
                np.full(
                    nstar,
                    0.5
                    * np.median(
                        calc_Lcell_pc(
                            m_gas_all,
                            np.asarray(gas["Density"], dtype="float64"),
                            cfg,
                        )
                        / 1e3
                    ),
                )
                if ngas
                else np.full(nstar, 0.5 * dx)
            )
        else:
            sig_s = np.zeros(nstar)
        sm = m_s * keep
        maps.update(
            deposit_gaussian_maps(
                xy_s,
                {"Sigma_star": sm, "Mvz": sm * vz, "Mvz2": sm * vz**2},
                sig_s,
                L_map,
                cfg,
                n,
            )
        )

    Sigma_n = maps.get("Sigma_n", np.zeros((n, n))).copy()
    Sigma_gas = maps.get("Sigma_gas", np.zeros((n, n))).copy()
    Sigma_star = maps.get("Sigma_star", np.zeros((n, n))).copy()
    Sigma_sfr = maps.get("Sigma_sfr", np.zeros((n, n))).copy()
    pix_area = dx**2

    # ---- 5. surface densities and the per-pixel metallicity ----------------
    with np.errstate(invalid="ignore", divide="ignore"):
        if z_weight == "neutral":
            Z_map = np.where(
                maps.get("Sigma_n", 0) > 0,
                maps.get("MZ_n", 0) / np.maximum(Sigma_n, 1e-30),
                Z_SUN,
            )
        else:
            Z_map = np.where(
                Sigma_gas > 0,
                maps.get("MZ_gas", 0) / np.maximum(Sigma_gas, 1e-30),
                Z_SUN,
            )
    Z_map = np.clip(np.nan_to_num(Z_map, nan=Z_SUN), 0.0, 0.1)

    # Units: the deposited maps are in TNG code units (1e10 Msun/h) per pixel,
    # while every f_mol formula wants physical Msun/kpc^2 (GD14's
    # Sigma_c = 5e7 Msun/kpc^2, K13's tau_c = 0.066 f_c D_MW Sigma/[Msun pc^-2]).
    # Skipping this makes Sigma small by 1e10/h and f_mol identically zero.
    # StarFormationRate is already physical Msun/yr, so it needs no conversion.
    mass_code_to_msun = 1e10 / cfg.h
    Sig_n = Sigma_n / pix_area * mass_code_to_msun
    Sig_star = Sigma_star / pix_area * mass_code_to_msun
    Sig_sfr = Sigma_sfr / pix_area

    # ---- 6. per-pixel velocity dispersion (two passes) ---------------------
    sigma_star_map = np.full((n, n), SIGMA_STAR_FLOOR_KMS)
    if nstar and "Mvz" in maps:
        mean_vz = maps["Mvz"] / np.maximum(Sigma_star, 1e-30)
        ij = _map_grid(xy_s, L_map, n)
        rel = np.asarray(stars["Velocities"], dtype="float64")[:, 2] - mean_vz[
            ij[:, 0], ij[:, 1]
        ]
        w2 = {"Mvz2": m_s * (sft > 0.0) * rel**2}
        var = deposit_gaussian_maps(xy_s, w2, sig_s, L_map, cfg, n)["Mvz2"]
        sigma_star_map = np.maximum(
            np.sqrt(var / np.maximum(Sigma_star, 1e-30)), SIGMA_STAR_FLOOR_KMS
        )

    # ---- 7. UV field and f_mol --------------------------------------------
    # The softening of the 2D 1/r^2 kernel: the paper gives the self-contribution
    # prescription for the 3D case only. The 0.25*dx default makes the
    # self-contribution diverge as pixels shrink, whereas physically it is capped
    # by the clump (cell) scale.
    r_soft_uv = None
    if cfg.r_soft_map_mode == "sigma":
        r_soft_uv = max(cfg.r_soft_factor_map * dx, sigma_eff)
    elif cfg.r_soft_map_mode == "cell":
        r_soft_uv = max(sigma_eff, 1e-6)
    U_MW = compute_galaxy_fft_uv_map(Sig_sfr, dx, cfg, r_soft_kpc=r_soft_uv)

    rho_dm = None
    if rho_dm_mode == "const":
        rho_dm = np.full((n, n), K13_RHO_DM_CONST)
    elif rho_dm_mode == "nfw" and sub.get("m200c_msun"):
        from .models import calc_rho_dm_nfw_map

        yy, xx = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
        Rk = (
            np.sqrt(
                (xx - (n - 1) / 2.0) ** 2 + (yy - (n - 1) / 2.0) ** 2
            )
            * dx
        )
        rho_dm = calc_rho_dm_nfw_map(Rk, sub["m200c_msun"], cfg)

    if model == "GD14":
        from .models import f_mol_GD14_map

        f_mol = f_mol_GD14_map(Sig_n, Z_map, U_MW, pix_pc)
    elif model == "K13":
        from .models import f_mol_K13_map

        f_mol = f_mol_K13_map(
            Sig_n,
            Z_map,
            U_MW,
            Sigma_star=Sig_star,
            sigma_star=sigma_star_map,
            rho_dm=rho_dm,
            variant=k13_variant,
            cfg=cfg,
        )
    else:
        raise ValueError(f"the projected models are GD14 and K13, got {model!r}")

    # ---- 8. integrate to a galaxy total ------------------------------------
    m_neutral = float(np.sum(Sigma_n))
    m_h2 = float(np.sum(f_mol * Sigma_n))
    m_hi = m_neutral - m_h2
    m_gas_in = float(np.sum(Sigma_gas))
    m_gas_tot = float(np.sum(m_gas_all)) if ngas else 0.0

    if _cells_out is not None:
        _cells_out.update(dict(xy=xy_g, L_map=L_map, n=n, f_mol=f_mol))
    out.update(
        {
            "m_hi": m_hi,
            "m_h2": m_h2,
            "m_neutral": m_neutral,
            "f_gas_in_map": (m_gas_in / m_gas_tot) if m_gas_tot > 0 else 0.0,
            "u_mw_map": U_MW,
            "sigma_eff_kpc": sigma_eff,
        }
    )
    return out


def _cells_vol(gas, model, cfg: HIModelConfig, uv_mode="adaptive_fft", uv=None):
    """Volume method: one f_mol per gas cell (Jeans approximation)."""
    from .physics import _gas_columns

    rho_code, u_code, nH_code, sfr_code, X_H, Z, m_gas = _gas_columns(gas)
    if uv is None:
        uv = compute_galaxy_fft_uv(gas, cfg, uvb_only=(uv_mode == "global"))
    f_mol_func = resolve(model)["fmol"]
    f_mol, f_HI, f_H2 = f_mol_func(
        rho_code, u_code, Z, nH_code, sfr_code, X_H, uv, cfg, mgas=m_gas
    )
    n = m_gas.size
    return Cells(
        row_index=np.arange(n, dtype=np.int64),
        f_mol=f_mol,
        f_neutral=get_f_neutral_H(nH_code, sfr_code),
        m_hi=m_gas * X_H * f_HI,
        m_h2=m_gas * X_H * f_H2,
        u_mw=np.asarray(uv, dtype="float64"),
    )


def _cells_map(
    gas,
    stars,
    sub,
    model,
    cfg: HIModelConfig,
    cell_fmol="bilinear",
    keep_outside=None,
    **kw,
):
    """Projected method: each cell takes the f_mol of **its pixel**."""
    keep_outside = cfg.map_clip_outside if keep_outside is None else bool(keep_outside)
    sink = {}
    out = compute_galaxy_hi_h2_map(
        gas, stars, sub, cfg, model=model, _cells_out=sink, **kw
    )
    n_cells = int(gas.get("count", 0)) if gas is not None else 0
    if "f_mol" not in sink:  # no gas cells at all
        z = np.zeros(n_cells)
        return Cells(
            row_index=np.zeros(0, dtype=np.int64),
            f_mol=z,
            f_neutral=z,
            m_hi=z,
            m_h2=z,
            u_mw=z,
            meta=out,
        )

    xy = sink["xy"]
    L_map = sink["L_map"]
    n = sink["n"]
    fmap = sink["f_mol"]
    m_n = sink["m_n"]
    f_neu = sink["f_neu"]

    half = L_map / 2.0
    inside = (np.abs(xy[:, 0]) < half) & (np.abs(xy[:, 1]) < half)
    keep = np.ones(xy.shape[0], dtype=bool) if keep_outside else inside
    idx = np.where(keep)[0]

    f_cell = _sample_map(fmap, xy[idx], L_map, n, cell_fmol)
    f_cell = np.clip(f_cell, 0.0, 1.0)

    # Build full-length arrays so they line up with the caller's gas arrays.
    # Cells outside the projection box get zero mass: the projected model's
    # galaxy *is* the projection box, and the catalog total is the in-box total.
    f_mol_full = np.zeros(n_cells)
    f_neutral_full = np.zeros(n_cells)
    m_hi_full = np.zeros(n_cells)
    m_h2_full = np.zeros(n_cells)
    f_mol_full[idx] = f_cell
    f_neutral_full[idx] = f_neu[idx]
    m_h2_full[idx] = m_n[idx] * f_cell
    m_hi_full[idx] = m_n[idx] * (1.0 - f_cell)

    meta = dict(out)
    meta.update(
        {
            "cell_fmol": cell_fmol,
            "n_cells_total": int(xy.shape[0]),
            "n_cells_kept": int(idx.size),
            "keep_outside": keep_outside,
        }
    )
    if out.get("m_h2"):
        meta["closure_m_h2"] = float(np.sum(m_h2_full) / out["m_h2"])
        meta["closure_m_hi"] = float(np.sum(m_hi_full) / out["m_hi"])
    return Cells(
        row_index=idx.astype(np.int64),
        f_mol=f_mol_full,
        f_neutral=f_neutral_full,
        m_hi=m_hi_full,
        m_h2=m_h2_full,
        u_mw=np.zeros(n_cells),
        meta=meta,
    )


def _sample_map(fmap, xy, L_map, n, cell_fmol):
    """Read f_mol off the map for the given positions."""
    if cell_fmol == "nearest":
        ij = _map_grid(xy, L_map, n)
        return fmap[ij[:, 0], ij[:, 1]]
    if cell_fmol == "bilinear":
        grid = (xy + L_map / 2.0) / (L_map / n) - 0.5  # pixel coords, centre-aligned
        return map_coordinates(
            fmap, [grid[:, 0], grid[:, 1]], order=1, mode="nearest"
        )
    raise ValueError(f"cell_fmol must be 'bilinear' or 'nearest', got {cell_fmol!r}")


@dataclass
class HiMassResult:
    """Per-gas-cell HI/H2 masses plus the diagnostics of one model evaluation.

    Attributes
    ----------
    model : str
        The method that produced this result.
    m_hi, m_h2 : ndarray
        Per-cell masses in TNG code units (1e10 Msun/h), the same unit as the
        ``Masses`` field. Always full length: cells the method excludes get
        zero, never a shorter array, so this can be handed straight to
        `~astrovodka.sources.SPHSource`.
    cells : Cells
        The same data in per-cell form.
    """

    model: str
    cells: Cells
    cfg: HIModelConfig
    diagnostics: dict = field(default_factory=dict)

    @property
    def m_hi(self):
        return self.cells.m_hi

    @property
    def m_h2(self):
        return self.cells.m_h2

    @property
    def m_neutral(self):
        return self.cells.m_neutral

    @property
    def m_hi_msun(self):
        return float(np.sum(self.cells.m_hi) * 1e10 / self.cfg.h)

    @property
    def m_h2_msun(self):
        return float(np.sum(self.cells.m_h2) * 1e10 / self.cfg.h)

    @property
    def m_neutral_msun(self):
        return float(np.sum(self.cells.m_neutral) * 1e10 / self.cfg.h)

    @property
    def f_mol(self):
        return self.cells.f_mol

    @property
    def u_mw(self):
        return self.cells.u_mw

    def __repr__(self):
        return (
            f"HiMassResult({self.model}, n_cells={self.cells.m_hi.size}, "
            f"m_hi={self.m_hi_msun:.4e} Msun, m_h2={self.m_h2_msun:.4e} Msun)"
        )


def compute_hi_mass(
    gas,
    model=DEFAULT_METHOD,
    cfg: HIModelConfig | None = None,
    *,
    stars=None,
    sub=None,
    uv_mode="adaptive_fft",
    uv=None,
    **kwargs,
):
    """Per-gas-cell HI and H2 masses for one subhalo, under any of the four models.

    This is the single keyword-selectable entry point. ``model`` is one of
    ``GD14_vol``, ``GD14_map``, ``K13_vol``, ``K13_map`` (case-insensitive; the
    bare model names ``GD14``/``K13`` mean the volume variant).

    Parameters
    ----------
    gas : dict
        TNG gas fields, at minimum ``Masses``, ``Density``, ``InternalEnergy``
        (volume methods), ``StarFormationRate``, ``NeutralHydrogenAbundance``,
        ``GFM_Metals``, ``GFM_Metallicity`` and ``Coordinates``.
    model : str, optional
        Which of the four methods to use. Defaults to ``K13_vol``.
    cfg : HIModelConfig, optional
        Runtime configuration. Defaults to ``HIModelConfig()`` with
        ``h = HUBBLE_DEFAULT`` and ``z = 0``.
    stars : dict, optional
        Required by the ``*_map`` methods (they need the stellar surface
        density and velocity dispersion).
    sub : dict, optional
        Required by the ``*_map`` methods: ``pos``, ``vel``, ``rhalf_gas``,
        ``rhalf_star``, optionally ``m200c_msun``.
    uv_mode : {'adaptive_fft', 'global'}
        Volume methods only. ``'global'`` skips the FFT and uses the uniform
        Faucher-Giguere+2011 background. Ignored by the projected methods,
        whose UV field is part of the method definition.
    uv : ndarray, optional
        A precomputed per-cell ``U_MW``, bypassing the UV calculation.

    Returns
    -------
    HiMassResult
    """
    spec = resolve(model)
    cfg = HIModelConfig() if cfg is None else cfg

    if spec["kind"] == "vol":
        if uv_mode not in ("adaptive_fft", "global"):
            raise ValueError(
                f"uv_mode must be 'adaptive_fft' or 'global', got {uv_mode!r}"
            )
        cells = _cells_vol(gas, spec["model"], cfg, uv_mode=uv_mode, uv=uv)
        diag = {"uv_mode": uv_mode}
    else:
        if stars is None or sub is None:
            raise ValueError(
                f"the {spec['label']} model projects the galaxy, so it needs the "
                "stellar particle data (stars=) and the subhalo record (sub=)"
            )
        cells = _cells_map(gas, stars, sub, spec["model"], cfg, **kwargs)
        diag = dict(cells.meta)

    return HiMassResult(model=spec["label"], cells=cells, cfg=cfg, diagnostics=diag)


__all__ = [
    "Cells",
    "HiMassResult",
    "compute_hi_mass",
    "compute_HI_mass_particle",
    "compute_galaxy_hi_h2_map",
]
