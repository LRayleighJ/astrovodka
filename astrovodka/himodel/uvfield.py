"""
UV radiation field (3D and 2D FFT convolutions) and projection helpers.

Ported line by line from the verified `TNGHIcalc` package (see `physics.py` for
the provenance and the bit-exactness contract); runtime parameters arrive as an
explicit `~astrovodka.himodel.config.HIModelConfig`.

The K13 production chain this code replaces had a bug here that is worth
remembering, because it failed *silently*: the 3D interpolation built its
``RegularGridInterpolator`` from ``edges[0]`` reused for all three axes. Each
axis' edge array is centred on that axis' own box centre, so for essentially
every particle two of the three coordinates fell outside the interpolator's
bounds and ``fill_value=0.0`` returned 0. The adaptive FFT field therefore
collapsed to the global UV background everywhere -- ``max(U_MW) == 0`` before
the UVB floor -- so ``uv_mode='adaptive_fft'`` and ``uv_mode='global'`` gave
numerically identical masses. `compute_galaxy_fft_uv` builds its three axis
grids separately, and `tests/test_hi_models.py` pins the off-centre regression
that catches a reintroduction.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator

from .config import HIModelConfig
from .constants import (
    ALIGN_RGAS_FACTOR,
    ALIGN_RSTR_FACTOR,
    BOX_MIN_KPC,
    MIN_ALIGN_GAS,
    MIN_ALIGN_STAR,
    Q_UV_CAL,
    R_SOFT_FACTOR,
    U_MW0_CAL,
)
from .physics import calc_UMW_global_Faucher11


def _fft_uv_box(coords_kpc, sfr, center, half_box, cfg: HIModelConfig, n=None):
    """FFT-convolve the SFR field with a 1/r^2 kernel in one given box.

    Returns ``(UMW grid, dx, edges)``. Factored out so the volume model's second,
    higher-resolution grid (paper Appendix A.2) can reuse the same bookkeeping
    rather than duplicating the FFT/softening/zero-lag details -- duplicating
    them is exactly how the K13 chain acquired its interpolation bug.
    """
    n = cfg.ngrid if n is None else int(n)
    dx = 2.0 * half_box / n
    edges = [
        np.linspace(center[i] - half_box, center[i] + half_box, n + 1) for i in range(3)
    ]

    sf_mask = sfr > 0.0
    if np.any(sf_mask):
        sfr_grid, _ = np.histogramdd(coords_kpc[sf_mask], bins=edges, weights=sfr[sf_mask])
    else:
        sfr_grid = np.zeros((n, n, n))

    x_1d = (np.arange(n) - n // 2) * dx
    XK, YK, ZK = np.meshgrid(x_1d, x_1d, x_1d, indexing="ij")
    RK = np.sqrt(XK**2 + YK**2 + ZK**2)
    RK_soft = np.maximum(RK, R_SOFT_FACTOR * dx)
    K_grid_real = Q_UV_CAL * (1.0 / RK_soft) ** 2  # photons/s/cm^2/Hz per Msun/yr

    F_grid = np.fft.ifftn(
        np.fft.fftn(sfr_grid) * np.fft.fftn(np.fft.ifftshift(K_grid_real))
    ).real
    return cfg.f_esc * F_grid / U_MW0_CAL, dx, edges


def _interp_grid(umw_grid, center, half_box, dx, pts, clip_inside=True):
    """Linearly interpolate a UV grid onto ``pts``.

    ``clip_inside=True`` clamps points into the box first (current behaviour,
    for the coarse grid); ``False`` returns 0 outside the box, which lets the
    fine grid fall back to the coarse one through a ``max``.
    """
    n = umw_grid.shape[0]
    x_c = [center[i] - half_box + (np.arange(n) + 0.5) * dx for i in range(3)]
    interp = RegularGridInterpolator(
        tuple(x_c), umw_grid, bounds_error=False, fill_value=0.0
    )
    if clip_inside:
        pts = np.clip(pts, center - half_box + dx, center + half_box - dx)
    return interp(pts)


def compute_galaxy_fft_uv(
    gas_data,
    cfg: HIModelConfig,
    center_kpc=None,
    rhalf_gas_kpc=None,
    return_grid=False,
    uvb_only=False,
):
    """Per-galaxy adaptive FFT UV field -> per-particle ``U_MW`` (Draine units).

    The recipe (Diemer+2018 §2.4 / Appendix A.2):

    1. TNG coordinates ckpc/h -> physical kpc: ``x_phys = x_code / h / (1+z)``
    2. adaptive cubic box: ``L_box = box_expand`` x the coordinate extent,
       floored at ``BOX_MIN_KPC``
    3. distribute the star-forming cells' SFR onto a 128^3 grid
    4. FFT-convolve with the optically thin 1/r^2 Green function (the kernel is
       ``Q_UV_CAL`` at 1 kpc)
    5. interpolate back to the particle positions to get the local ``U_MW``
    6. take the **maximum** with the global UVB (paper Appendix A, opening:
       "set the UV field to the maximum of two components")

    Step 5 must build its three interpolation axes separately; see the module
    docstring for what happens otherwise.

    Parameters
    ----------
    gas_data : dict
        TNG gas fields, needing ``Coordinates`` (or ``CenterOfMass``) and
        ``StarFormationRate``.
    cfg : HIModelConfig
        Runtime configuration.
    center_kpc : array-like, optional
        Centre of the second, finer UV box [physical kpc]. Only used when
        ``cfg.vol_fine_grid > 0``.
    rhalf_gas_kpc : float, optional
        Gas half-mass radius [physical kpc], setting the fine box's half-width.
    return_grid : bool, optional
        Also return the raw ``(ngrid, ngrid, ngrid)`` field and its edges.
    uvb_only : bool, optional
        Skip the FFT entirely and return the uniform global UV background. This
        is the ``uv_mode="global"`` path.

    Returns
    -------
    U_MW : ndarray
        Per-particle UV field in Draine units.
    grid, edges
        Only when ``return_grid`` is True.
    """
    if uvb_only:
        uvb = calc_UMW_global_Faucher11(cfg.z)
        field = np.full(gas_data["StarFormationRate"][:].shape, uvb, dtype="float64")
        if return_grid:
            n = cfg.ngrid
            lo, hi = 0.0, float(BOX_MIN_KPC)
            edges = [np.linspace(lo, hi, n + 1) for _ in range(3)]
            return field, np.full((n, n, n), uvb), edges
        return field

    coords = _gas_coords(gas_data)
    sfr = gas_data["StarFormationRate"][:].astype("float64")

    # ckpc/h -> physical kpc
    coords_kpc = coords / cfg.h / (1.0 + cfg.z)

    # Adaptive box (coarse grid: spans all the gas bound to this subhalo)
    coord_min = coords_kpc.min(axis=0)
    coord_max = coords_kpc.max(axis=0)
    L_box = max(cfg.box_expand * float(np.max(coord_max - coord_min)), BOX_MIN_KPC)
    center = (coord_min + coord_max) / 2.0
    half_box = L_box / 2.0

    UMW_raw, dx, edges = _fft_uv_box(coords_kpc, sfr, center, half_box, cfg, cfg.ngrid)
    # Interpolate back to the particle positions -- one grid per axis!
    UMW_per_particle = _interp_grid(
        UMW_raw, center, half_box, dx, coords_kpc, clip_inside=True
    )

    # ---- second, higher-resolution grid (paper App. A.2, off by default) ----
    # "We perform this calculation twice on two 128^3 grids, one grid that spans
    #  all gas cells bound to the galaxy and one higher-resolution grid that
    #  spans two gas half-mass radii. Without the high-resolution grid, the UV at
    #  the center of the galaxy can be underestimated significantly."
    # The paper does not say how to combine the two, so take the max: the fine
    # grid returns 0 outside its box and therefore falls back automatically.
    if cfg.vol_fine_grid > 0.0 and center_kpc is not None and rhalf_gas_kpc:
        c_fine = np.asarray(center_kpc, dtype="float64")
        half_fine = cfg.vol_fine_grid * float(rhalf_gas_kpc)
        if half_fine > 0.0:
            UMW_fine, dxf, _ = _fft_uv_box(
                coords_kpc, sfr, c_fine, half_fine, cfg, cfg.ngrid
            )
            UMW_per_particle = np.maximum(
                UMW_per_particle,
                _interp_grid(UMW_fine, c_fine, half_fine, dxf, coords_kpc, clip_inside=False),
            )

    # The global UVB is a floor (the paper takes a max, not a sum)
    UMW_UVB = calc_UMW_global_Faucher11(cfg.z)
    UMW_per_particle = np.maximum(UMW_per_particle, UMW_UVB)

    if return_grid:
        return UMW_per_particle, UMW_raw, edges
    return UMW_per_particle


def _gas_coords(gas_data):
    """Gas coordinates in raw ckpc/h, accepting either TNG field name."""
    if "Coordinates" in gas_data:
        return gas_data["Coordinates"][:].astype("float64")
    if "CenterOfMass" in gas_data:
        return gas_data["CenterOfMass"][:].astype("float64")
    raise KeyError(
        "gas data must contain 'Coordinates' or 'CenterOfMass' to place particles"
    )


def _rotation_to_z(vhat):
    """Rotation matrix taking the unit vector ``vhat`` to +z (shortest arc)."""
    z = np.array([0.0, 0.0, 1.0])
    c = float(np.dot(vhat, z))
    if c > 1.0 - 1e-12:
        return np.eye(3)
    if c < -1.0 + 1e-12:
        # Antiparallel: the axis is not unique, so any 180 deg rotation will do
        return np.diag([1.0, -1.0, -1.0])
    axis = np.cross(vhat, z)
    s = float(np.linalg.norm(axis))
    axis = axis / s
    K = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    return np.eye(3) + s * K + (1.0 - c) * (K @ K)


def compute_align_vector(
    gas, stars, subhalo_pos, subhalo_vel, rhalf_gas_ckpc, rhalf_star_ckpc
):
    """Face-on alignment axis. Returns ``(vhat, success, info)``.

    Conventions, calibrated per galaxy against this repository's official
    catalog ``align_ang_mom_vector`` (over 16 galaxies the signed dot product has
    median 0.9956 and 94% have |dot| > 0.9):

      * origin      = ``SubhaloPos``, in **raw ckpc/h** (the same units as
        ``Coordinates``)
      * velocity reference = ``SubhaloVel``
      * gas radius  = ``ALIGN_RGAS_FACTOR`` x R_half_gas (the paper's "within
        half the 3D gas half-mass radius")
      * weight      = cell mass (equivalent to unweighted for Arepo's
        near-equal-mass cells)
      * ``L = sum m (r x (v - v_sub))``, normalised; the sign is **not** flipped
        (measurement says the positive sign is the right one)

    The radius here must be raw ckpc/h, while ``L_map`` must be physical kpc.
    The same numbers in two unit systems; mixing them costs ~2% of the aligned
    mass in practice.
    """
    info = {"n_gas_align": 0, "n_star_align": 0, "branch": "none"}
    z = np.array([0.0, 0.0, 1.0])

    def _n_in(part, rad, mask=None):
        if part is None or part.get("count", 0) == 0:
            return None, None, 0
        c = np.asarray(part["Coordinates"], dtype="float64") - subhalo_pos
        sel = np.einsum("ij,ij->i", c, c) < float(rad) ** 2
        if mask is not None:
            sel = sel & mask
        return c, sel, int(sel.sum())

    cg, sg, ng = _n_in(gas, ALIGN_RGAS_FACTOR * rhalf_gas_ckpc)
    info["n_gas_align"] = ng

    star_mask = None
    cs = ss = None
    ns = 0
    if stars is not None and stars.get("count", 0) > 0:
        # Exclude wind particles (negative GFM_StellarFormationTime)
        star_mask = np.asarray(stars["GFM_StellarFormationTime"], dtype="float64") > 0.0
        cs, ss, ns = _n_in(stars, ALIGN_RSTR_FACTOR * rhalf_star_ckpc, mask=star_mask)
    info["n_star_align"] = ns

    if ng >= MIN_ALIGN_GAS:
        v = np.asarray(gas["Velocities"], dtype="float64") - subhalo_vel
        m = np.asarray(gas["Masses"], dtype="float64")
        L = np.sum(m[sg, None] * np.cross(cg[sg], v[sg]), axis=0)
        n = float(np.linalg.norm(L))
        if n > 0:
            info["branch"] = "gas"
            return L / n, True, info

    if ns >= MIN_ALIGN_STAR:
        v = np.asarray(stars["Velocities"], dtype="float64") - subhalo_vel
        m = np.asarray(stars["Masses"], dtype="float64")
        L = np.sum(m[ss, None] * np.cross(cs[ss], v[ss]), axis=0)
        n = float(np.linalg.norm(L))
        if n > 0:
            info["branch"] = "star"
            return L / n, True, info

    # Neither branch has enough particles. The official catalog still assigns a
    # projected mass to rows with align_success=0, so this is not fatal: skip
    # the rotation and carry on.
    return z.copy(), False, info


def deposit_gaussian_maps(
    xy_kpc,
    weights,
    sigma_kpc,
    L_map_kpc,
    cfg: HIModelConfig,
    n=None,
    sub=None,
    sigma_bins=None,
    clip_outside=None,
):
    """Deposit weighted points onto an n^2 grid, each with its own Gaussian.

    Every gas cell has its own kernel width ``sigma = L_cell/2`` (each cell has
    its own ``(m/rho)^(1/3)``), so a single ``gaussian_filter`` over the whole
    map is not equivalent: a typical ``sigma`` of ~300 pc is smaller than the
    median pixel of 799 pc, and whole-image smoothing degenerates into
    "one point, one pixel", while the paper is explicit that bare
    particle-in-cell deposition systematically overestimates f_mol. Instead the
    widths are log-binned, supersampled by (sub*n)^2, block-summed back down to
    n^2, and each bin gets its own kernel.

    The width bins are computed once and shared by every quantity in
    ``weights``: quantities that form a ratio must share a kernel or the ratio
    is formed between mismatched maps.

    Parameters
    ----------
    xy_kpc : ndarray, shape (N, 2)
        Physical-kpc coordinates relative to the grid centre, already rotated.
    weights : dict of str -> ndarray
        Several weight sets over the same points.
    sigma_kpc : ndarray
        Per-point kernel width; a value <= 0 deposits exactly (no smoothing).
    L_map_kpc : float
        Map edge length in physical kpc.

    Returns
    -------
    dict of str -> ndarray, shape (n, n)
        The **total** weight per pixel (not divided by the pixel area).
    """
    n = cfg.ngrid_map if n is None else int(n)
    sub = cfg.map_sub if sub is None else int(sub)
    sigma_bins = cfg.map_sigma_bins if sigma_bins is None else int(sigma_bins)
    n_sub = n * sub
    dx_sub = L_map_kpc / n_sub
    names = list(weights)
    empty = {k: np.zeros((n, n), dtype="float64") for k in names}
    if len(xy_kpc) == 0:
        return empty

    half = L_map_kpc / 2.0
    clip_outside = cfg.map_clip_outside if clip_outside is None else bool(clip_outside)
    if clip_outside:
        # Cells outside the box are clamped to the edge pixel (the official map
        # loses no mass -- see HIModelConfig.map_clip_outside)
        eps = 1e-6 * L_map_kpc
        x = np.clip(xy_kpc[:, 0], -half + eps, half - eps) + half
        y = np.clip(xy_kpc[:, 1], -half + eps, half - eps) + half
        keep = np.ones(len(xy_kpc), dtype=bool)
    else:
        keep = (np.abs(xy_kpc[:, 0]) < half) & (np.abs(xy_kpc[:, 1]) < half)
        if not np.any(keep):
            return empty
        x = xy_kpc[keep, 0] + half
        y = xy_kpc[keep, 1] + half
    s = np.maximum(sigma_kpc[keep], 0.0)
    w = {k: np.asarray(weights[k], dtype="float64")[keep] for k in names}

    rng = ((0.0, L_map_kpc), (0.0, L_map_kpc))
    acc = {k: np.zeros((n_sub, n_sub), dtype="float64") for k in names}

    def _add(mask, sig_kpc):
        if not np.any(mask):
            return
        h = {
            k: np.histogram2d(x[mask], y[mask], bins=n_sub, range=rng, weights=w[k][mask])[0]
            for k in names
        }
        if sig_kpc > 0:
            sp = sig_kpc / dx_sub
            for k in names:
                acc[k] += gaussian_filter(h[k], sigma=sp, mode="constant", truncate=4.0)
        else:
            for k in names:
                acc[k] += h[k]

    pos = s > 0
    if np.any(pos):
        lo = max(float(np.percentile(s[pos], 0.5)), 1e-4 * L_map_kpc)
        hi = max(float(np.percentile(s[pos], 99.5)), lo * 1.0001)
        edges = np.logspace(np.log10(lo), np.log10(hi), sigma_bins + 1)
        which = np.clip(np.digitize(s, edges) - 1, 0, sigma_bins - 1)
        for b in range(sigma_bins):
            m = pos & (which == b)
            _add(m, 0.5 * (edges[b] + edges[b + 1]))
    neg = ~pos
    if np.any(neg):
        _add(neg, 0.0)

    return {k: acc[k].reshape(n, sub, n, sub).sum(axis=(1, 3)) for k in names}


def compute_galaxy_fft_uv_map(sigma_sfr, dx_kpc, cfg: HIModelConfig, r_soft_kpc=None):
    """2D FFT convolution: SFR surface-density map -> per-pixel ``U_MW``.

    The paper's §2.4 closes with "In the case of projected modeling, we simply
    apply the FFT to the 2D map of star formation surface density."

    Three differences from the 3D version:

    1. **The source term is ``sigma_sfr * dx^2`` (Msun/yr per pixel), not the
       surface density.** Getting this wrong is a factor ``dx^2`` (4 to 2600)
       and it fails silently -- the map still looks like a UV field, only the
       normalisation is off. The self-test catches this one specifically.
    2. The 1/r^2 kernel is 2D; ``Q_UV_CAL`` is already "flux at 1 kpc", so it
       carries over unchanged.
    3. The normalisation chain is line-for-line the 3D one:
       ``max(f_esc * flux / U_MW0, global UVB)``.
    """
    z = cfg.z
    n = int(sigma_sfr.shape[0])
    sfr_pix = np.asarray(sigma_sfr, dtype="float64") * float(dx_kpc) ** 2

    x_1d = (np.arange(n) - n // 2) * float(dx_kpc)
    XK, YK = np.meshgrid(x_1d, x_1d, indexing="ij")
    RK = np.sqrt(XK**2 + YK**2)
    # The softening of the 2D kernel can be supplied by the caller (cell scale);
    # the default remains 0.25 * dx_pix
    soft = (
        cfg.r_soft_factor_map * float(dx_kpc)
        if r_soft_kpc is None
        else float(r_soft_kpc)
    )
    RK_soft = np.maximum(RK, soft)
    K2 = Q_UV_CAL * (1.0 / RK_soft) ** 2

    F = np.fft.ifftn(np.fft.fftn(sfr_pix) * np.fft.fftn(np.fft.ifftshift(K2))).real
    UMW = cfg.f_esc * F / U_MW0_CAL
    return np.maximum(UMW, calc_UMW_global_Faucher11(z))


def _map_grid(xy_kpc, L_map_kpc, n):
    """Map physical-kpc coordinates to pixel indices (nearest-neighbour lookup)."""
    idx = np.floor((xy_kpc + L_map_kpc / 2.0) / (L_map_kpc / n)).astype("int64")
    return np.clip(idx, 0, n - 1)


__all__ = [
    "compute_galaxy_fft_uv",
    "compute_galaxy_fft_uv_map",
    "compute_align_vector",
    "deposit_gaussian_maps",
    "_rotation_to_z",
    "_map_grid",
]
