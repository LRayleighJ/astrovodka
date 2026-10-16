"""
Group-scale (FOF group / cluster) 3D data products for TNG.

`SubhaloDataProducts` grids *one* subhalo. This module grids a whole
**FOF group**: every subhalo it contains, plus the group's unbound "fuzz" gas,
classified into

===========  ==========================================================
``ISM``      bound to a subhalo and within ``k x R_ISM`` of its centre
``CGM``      bound to a subhalo, outside that radius (no R_200c cut)
``IGM``      in the group but bound to no subhalo (the fuzz)
===========  ==========================================================

Membership comes from the group catalog offsets rather than from a spatial
guess, and it is exact: the gas cells of subhalo ``s`` occupy the contiguous
range ``[Subhalo/SnapByType[s], + SubhaloLenType[s])`` of the snapshot's
``PartType0`` array, and a group's own range is tiled by those blocks followed
by its fuzz.

Two details of the offsets file are worth knowing before touching
`_owner_index`:

1. **Lengths must come from ``SubhaloLenType``, never from differencing the
   offsets array.** ``SnapByType`` is a cumulative start-offset table and, for
   the *last* subhalo of a group, its next row already belongs to the next
   group, so the difference silently swallows that group's fuzz. Measured on
   TNG-100 snapshot 99 group 59: the last subhalo has ``SubhaloLenType == 0``
   yet a 319,249-cell step. Summing block masses with the catalog lengths
   reproduces ``SubhaloMassType`` at a ratio of 1.0000.
2. The cells of a group are ordered ``[subhalo blocks in catalog order][fuzz]``,
   but `_owner_index` places each block by its absolute offset, so leading fuzz
   -- should a group ever have any -- is handled correctly too.

The fuzz is a real component, not a rounding error: 1.5% of the gas mass in the
most massive TNG-100 group, ~10% at the mass of group 59, and ~24% in smaller
groups. It is also *hot and diffuse* (mass-weighted median T ~ 2.4e6 K, n_H ~
2e-5 cm^-3 for group 59), so it contributes essentially no HI at any UV field:
`from_group` models it anyway, for completeness.
"""

from __future__ import annotations

import numpy as np
from astropy.cosmology import Planck15

from .data_products import _cell_volume
from .himodel import (
    DEFAULT_METHOD,
    GAS_FIELDS_FULL,
    HIModelConfig,
    compute_hi_mass,
    is_known,
    resolve,
)
from .himodel.physics import calc_temperature, rho_code_to_nH
from .himodel.uvfield import compute_galaxy_fft_uv
from .sources._illustris_tools import (
    gcPath,
    getSnapOffsets,
    loadHeader,
    loadSingle,
    loadSubset,
    offsetPath,
)

ISM = 0
CGM = 1
IGM = 2

CLASS_LABELS = {ISM: "ISM", CGM: "CGM", IGM: "IGM"}

# Fields the volumetric models read, plus the electron abundance the
# temperature needs; velocities are here for the redshift-space grids.
_GROUP_GAS_FIELDS = list(GAS_FIELDS_FULL) + ["ElectronAbundance"]
_GROUP_GAS_MDI = [0 if f == "GFM_Metals" else None for f in _GROUP_GAS_FIELDS]
_DM_FIELDS = ["Coordinates", "Velocities"]

# (weight key, class filter or None). "dm_mass" weights the dark matter, the
# rest the gas; a class filter selects one of ISM/CGM/IGM.
_REAL_GRIDS = {
    "dm_real": ("dm_mass", None),
    "dm_igm_real": ("dm_mass", IGM),
    "gas_real": ("m_gas", None),
    "gas_ism_real": ("m_gas", ISM),
    "gas_cgm_real": ("m_gas", CGM),
    "gas_igm_real": ("m_gas", IGM),
    "hi_real": ("m_hi", None),
    "hi_ism_real": ("m_hi", ISM),
    "hi_cgm_real": ("m_hi", CGM),
    "hi_igm_real": ("m_hi", IGM),
    "h2_real": ("m_h2", None),
}
_RSD_GRIDS = {
    "dm_rsd": ("dm_mass", None),
    "gas_rsd": ("m_gas", None),
    "hi_rsd": ("m_hi", None),
    "hi_ism_rsd": ("m_hi", ISM),
    "hi_cgm_rsd": ("m_hi", CGM),
    "hi_igm_rsd": ("m_hi", IGM),
}
_MASS_GRIDS = {**_REAL_GRIDS, **_RSD_GRIDS}
_MEAN_GRIDS = ["t_mw", "nh_mw"]
_ALL_GRIDS = list(_MASS_GRIDS) + _MEAN_GRIDS


def ism_radius(rad_gas, rad_star, k_gas=1.0, k_star=2.0):
    """Effective ISM radius of each subhalo [same unit as the inputs].

    ``min(k_gas * R_half,gas, k_star * R_half,star)``, from
    ``SubhaloHalfmassRadType`` columns 0 and 4. The ``min`` is what keeps the
    definition usable at group scale: a group's central subhalo has an
    ``R_half,gas`` of hundreds of ckpc/h because it contains the whole hot halo
    (measured: 263 ckpc/h for TNG-100 group 59, against 19 ckpc/h for
    ``R_half,star``), so a gas-only radius would swallow most of the group.
    A non-positive or non-finite radius is treated as "no constraint"; if
    neither radius is usable the result is 0 and that subhalo contributes no
    ISM.

    Parameters
    ----------
    rad_gas, rad_star : array-like
        Gas and stellar half-mass radii [ckpc/h].
    k_gas, k_star : float, optional
        Multipliers, defaulting to 1 and 2.

    Returns
    -------
    ~numpy.ndarray
        ISM radius per subhalo [ckpc/h].
    """
    r_gas = np.asarray(rad_gas, dtype="float64") * float(k_gas)
    r_star = np.asarray(rad_star, dtype="float64") * float(k_star)
    r_gas = np.where(np.isfinite(r_gas) & (r_gas > 0), r_gas, np.inf)
    r_star = np.where(np.isfinite(r_star) & (r_star > 0), r_star, np.inf)
    out = np.minimum(r_gas, r_star)
    return np.where(np.isfinite(out), out, 0.0)


def classify_ism_cgm_igm(owner, r_gal, r_ism, nH=None, T=None, nH_max=None, T_max=None):
    """Classify gas cells as ISM, CGM or IGM.

    Parameters
    ----------
    owner : ndarray of int
        Index of the subhalo that owns each cell, ``-1`` for the group's fuzz
        (see `_owner_index`).
    r_gal : ndarray
        Distance from each cell to the centre of its owning subhalo [ckpc/h];
        ignored where ``owner < 0``.
    r_ism : ndarray or float
        ISM radius of the owning subhalo, already indexed by cell (or a
        scalar).
    nH : ndarray, optional
        Hydrogen number density [cm^-3]; required with ``nH_max``.
    T : ndarray, optional
        Temperature [K]; required with ``T_max``.
    nH_max, T_max : float, optional
        Optional phase caps: a cell failing either is demoted out of the ISM.
        Both default to None. Note that the cap is a no-op for the population it
        targets -- dense *and* hot shocked gas (n_H > 0.1 cm^-3 with T > 1e6 K)
        is 0.000 of the ISM mass in TNG-100 group 59 -- which is why it is off.

    Returns
    -------
    ndarray of int8
        ``ISM`` (0), ``CGM`` (1) or ``IGM`` (2) for every cell.
    """
    owner = np.asarray(owner)
    ism = owner >= 0
    if r_ism is not None:
        ism = ism & (np.asarray(r_gal) < np.asarray(r_ism))
    if nH_max is not None:
        if nH is None:
            raise ValueError("nH_max was given but nH was not")
        ism = ism & (np.asarray(nH) < float(nH_max))
    if T_max is not None:
        if T is None:
            raise ValueError("T_max was given but T was not")
        ism = ism & (np.asarray(T) < float(T_max))
    out = np.where(owner >= 0, CGM, IGM).astype("int8")
    out[ism] = ISM
    return out


def load_group_members(basePath, snapNum, groupID):
    """Catalog records of every subhalo in a FOF group.

    A group's subhalos occupy the contiguous catalog rows
    ``[GroupFirstSub, + GroupNsubs)``, so this is a single chunk read.

    Parameters
    ----------
    basePath : str
        TNG output directory (the one containing ``groups_<snapNum>/``).
    snapNum : int
        Snapshot number.
    groupID : int
        FOF group ID (``haloID`` in the group catalog).

    Returns
    -------
    dict
        ``group`` (the FOF group's own catalog record), ``first`` (the group's
        first subhalo ID), ``pos`` (n, 3) [ckpc/h], ``len_type`` (n, 6),
        ``mass_type`` (n, 6) [1e10 Msun/h], ``mass`` (n), ``rhalf_type`` (n, 6)
        [ckpc/h], ``sfr`` (n), ``flag`` (n), and ``start`` (n, 6): the start
        offset of each subhalo's particles of every type in the snapshot array.
    """
    import h5py

    group = loadSingle(basePath, snapNum, haloID=int(groupID))
    first = int(group["GroupFirstSub"])
    nsub = int(group["GroupNsubs"])

    with h5py.File(offsetPath(basePath, snapNum), "r") as f:
        file_offsets = np.asarray(f["FileOffsets/Subhalo"][()], dtype="int64")
        start = (
            np.asarray(f["Subhalo/SnapByType"][first : first + nsub, :], dtype="int64")
            if nsub
            else np.zeros((0, 6), dtype="int64")
        )

    rec = {
        "group": group,
        "first": first,
        "start": start,
        "pos": np.zeros((nsub, 3)),
        "len_type": np.zeros((nsub, 6), dtype="int64"),
        "mass_type": np.zeros((nsub, 6)),
        "mass": np.zeros(nsub),
        "rhalf_type": np.zeros((nsub, 6)),
        "sfr": np.zeros(nsub),
        "flag": np.zeros(nsub, dtype="int64"),
    }
    if nsub == 0:
        return rec

    chunk = int(np.max(np.where(first - file_offsets >= 0)))
    local = first - int(file_offsets[chunk])
    with h5py.File(gcPath(basePath, snapNum, chunk), "r") as f:
        s = f["Subhalo"]
        rec["pos"] = np.asarray(s["SubhaloPos"][local : local + nsub], dtype="float64")
        rec["len_type"] = np.asarray(s["SubhaloLenType"][local : local + nsub])
        rec["mass_type"] = np.asarray(s["SubhaloMassType"][local : local + nsub])
        rec["mass"] = np.asarray(s["SubhaloMass"][local : local + nsub], dtype="float64")
        rec["rhalf_type"] = np.asarray(
            s["SubhaloHalfmassRadType"][local : local + nsub], dtype="float64"
        )
        rec["sfr"] = np.asarray(s["SubhaloSFR"][local : local + nsub], dtype="float64")
        rec["flag"] = np.asarray(s["SubhaloFlag"][local : local + nsub])
    return rec


def _load_group_block(basePath, snapNum, groupID, ptype, fields, mdi):
    """Load a group's contiguous particle block of one type.

    Returns the field dict and the block's start offset in the snapshot array,
    which `_owner_index` needs in order to place the subhalo blocks.
    """
    subset = getSnapOffsets(basePath, snapNum, int(groupID), "Group")
    data = loadSubset(
        basePath, snapNum, ptype, fields=fields, subset=subset, mdi=mdi, sq=False
    )
    return data, int(subset["offsetType"][ptype])


def _owner_index(start, lengths, group_start, n_cells, label):
    """Map each cell of a group's block to the subhalo that owns it.

    ``start`` holds absolute snapshot offsets (``Subhalo/SnapByType``) and
    ``lengths`` the catalog counts (``SubhaloLenType``) -- see the module
    docstring for why the two must not be confused.

    Returns
    -------
    owner : ndarray of int32
        Subhalo index per cell, ``-1`` for the group's fuzz.
    blocks : ndarray of int64
        Per-subhalo start index within the loaded block, ``-1`` if it has no
        cells of this type.
    """
    owner = np.full(int(n_cells), -1, dtype="int32")
    blocks = np.full(len(lengths), -1, dtype="int64")
    for s in np.flatnonzero(np.asarray(lengths) > 0):
        a = int(start[s]) - int(group_start)
        b = a + int(lengths[s])
        if a < 0 or b > n_cells:
            raise RuntimeError(
                f"{label}: the block of subhalo {s} spans [{a}, {b}) of the "
                f"group's {n_cells} cells; the offsets file and the group "
                "catalog disagree about this group's layout"
            )
        owner[a:b] = s
        blocks[s] = a
    assigned = int((owner >= 0).sum())
    if assigned != int(np.sum(lengths)):
        raise RuntimeError(
            f"{label}: subhalo blocks overlap ({assigned} cells assigned for "
            f"{int(np.sum(lengths))} catalog cells)"
        )
    return owner, blocks


def _check_ownership(owner, blocks, lengths, masses, catalog_mass, label, groupID):
    """Cross-check block membership against the catalog's per-subhalo masses.

    This is the end-to-end evidence that the offsets convention still holds:
    sum the cell masses of each block and compare with ``SubhaloMassType``.
    """
    bad = []
    for s in np.flatnonzero(np.asarray(lengths) > 0):
        if catalog_mass[s] <= 0:
            continue
        a = int(blocks[s])
        ratio = masses[a : a + int(lengths[s])].sum() / catalog_mass[s]
        if not 0.99 < ratio < 1.01:
            bad.append((int(s), float(ratio)))
    if bad:
        raise RuntimeError(
            f"subhalo membership does not reproduce the catalog masses for group "
            f"{groupID} ({label}): {len(bad)} subhalo(s) off, worst is subhalo "
            f"{bad[0][0]} at ratio {bad[0][1]:.4f}. The offsets file and the "
            "snapshot layout have changed; refusing to continue."
        )


def _grid_sum(coords, weights, edges):
    """Histogram a weighted point set, in Msun / (ckpc/h)^3."""
    hist, _ = np.histogramdd(np.asarray(coords), bins=edges, weights=np.asarray(weights))
    return (hist / _cell_volume(edges)).astype("float32")


def _grid_mean(coords, weights, values, edges):
    """Mass-weighted mean of a per-cell scalar, per grid cell."""
    num, _ = np.histogramdd(
        np.asarray(coords), bins=edges, weights=np.asarray(weights) * np.asarray(values)
    )
    den, _ = np.histogramdd(np.asarray(coords), bins=edges, weights=np.asarray(weights))
    out = np.zeros_like(den)
    np.divide(num, den, out=out, where=den > 0)
    return out.astype("float32")


class GroupDataProducts:
    """3D gridded gas and dark-matter products for one TNG FOF group.

    Grids are comoving mass densities in ``Msun / (ckpc/h)^3`` on an
    ``(ngrid, ngrid, ngrid)`` box centred on the group, with ``x`` along the
    line of sight and ``y-z`` the sky plane, matching
    `~astrovodka.data_products.SubhaloDataProducts`. Mass grids come in
    real-space and redshift-space versions; the mass-weighted temperature and
    density are per-cell scalars and therefore real-space only.

    Attributes
    ----------
    ngrid : int
        Edge length of the grids.
    dm_real, dm_rsd : ndarray
        Total dark-matter mass density.
    dm_igm_real : ndarray
        Dark matter belonging to no subhalo.
    gas_real, gas_rsd : ndarray
        Total gas mass density.
    gas_ism_real, gas_cgm_real, gas_igm_real : ndarray
        Gas mass density of each class.
    hi_real, hi_rsd : ndarray
        Neutral hydrogen mass density.
    hi_ism_real, hi_ism_rsd, hi_cgm_real, hi_cgm_rsd, hi_igm_real, hi_igm_rsd : ndarray
        Neutral hydrogen mass density of each class.
    h2_real : ndarray
        Molecular hydrogen mass density.
    t_mw, nh_mw : ndarray
        Gas-mass-weighted temperature [K] and hydrogen number density [cm^-3].
    grid_edges : tuple of ndarray
        Three 1D edge arrays [ckpc/h].
    metadata : dict
        Group identity, mass budget, classification settings and grid bounds.
    cells : dict
        Per-cell arrays (positions, masses, class, owner, density, temperature,
        HI/H2), so that a re-classification or a custom profile does not need the
        snapshot again.
    subhalos : dict
        The member subhalo catalog table, as returned by `load_group_members`.
    """

    def __init__(self, ngrid=256):
        self.ngrid = int(ngrid)
        for name in _ALL_GRIDS:
            setattr(self, name, None)
        self.grid_edges = None
        self.grid_edges_rsd = None
        self.metadata = {}
        self.cells = {}
        self.subhalos = {}

    # ------------------------------------------------------------------ build
    @classmethod
    def from_group(
        cls,
        basePath,
        snapNum,
        groupID,
        *,
        ngrid=256,
        hi_model=DEFAULT_METHOD,
        hi_config=None,
        uv_mode="adaptive_fft",
        compute_hi=True,
        compute_rsd=False,
        k_gas=1.0,
        k_star=2.0,
        nH_max=None,
        T_max=None,
        igm_uv="group",
        ncpu=1,
        cosmology=None,
        box_halfsize_ckpc=None,
    ):
        """Load a FOF group and compute all of its 3D products.

        Parameters
        ----------
        basePath : str
            TNG output directory: the one holding ``groups_<snapNum>/`` and
            ``snapdir_<snapNum>/``, with the ``postprocessing/offsets`` files
            beside its parent, as `~astrovodka.sources._illustris_tools.offsetPath`
            expects.
        snapNum : int
            Snapshot number.
        groupID : int
            FOF group ID (``haloID``).
        ngrid : int, optional
            Edge length of the 3D grids (default 256). Memory is ``4 *
            ngrid^3`` bytes per grid; use 128 for the most massive groups.
        hi_model : str, optional
            One of the four Diemer+2018 models, but only the volumetric
            ``*_vol`` ones: the projected models build a face-on map about a
            single galaxy, which a group does not have. Default ``K13_vol``.
        hi_config : HIModelConfig or dict, optional
            Model settings; ``h``, ``z`` and the box size always come from the
            snapshot header.
        uv_mode : {'adaptive_fft', 'global'}, optional
            UV field for the volumetric models. Under ``'adaptive_fft'`` the FFT
            is run **per subhalo**, so each galaxy keeps the adaptive box the
            model was calibrated with -- doing it once for the whole group would
            smear the field over ~30 kpc cells at cluster scale.
        compute_hi : bool, optional
            Set False to skip the HI/H2 model (gas and dark matter only). The HI
            grids are then None and the per-cell ``m_hi``/``m_h2`` are zero.
        compute_rsd : bool, optional
            Also build the redshift-space grids. Off by default, and the reason
            is geometric rather than a matter of cost: at z = 0 one km/s of
            line-of-sight velocity is 10 ckpc/h, so a group whose members span
            a few hundred km/s occupies several Mpc along the line of sight in
            redshift space -- measured, group 59 needs a 17 cMpc axis against
            1.9 cMpc in the sky plane, and the grids are stretched to match.
            Pass True when the RSD is the point (the box then grows to fit, and
            ``f_mass_in_rsd_box_*`` reports what a pinned box would keep).
        k_gas, k_star : float, optional
            ISM radius multipliers; see `ism_radius`.
        nH_max, T_max : float, optional
            Optional phase caps on the ISM class; see `classify_ism_cgm_igm`.
        igm_uv : {'group', 'uvb'}, optional
            UV field for the fuzz, which belongs to no subhalo. ``'group'``
            (default) runs one more FFT over the whole group's star formation
            and interpolates it onto the fuzz cells; ``'uvb'`` uses the uniform
            Faucher-Giguere+2011 background. Both give essentially no HI -- the
            fuzz is hot and diffuse -- but ``'group'`` is the better motivated.
        ncpu : int, optional
            Processes for the per-subhalo HI loop. ``ncpu > 1`` needs
            ``multiprocess`` and a ``fork`` start method (the gas arrays reach
            the workers by inheritance, not by pickling).
        cosmology : astropy.cosmology.Cosmology, optional
            Cosmology for the redshift-space shift; defaults to Planck15.
        box_halfsize_ckpc : float, optional
            Half-width of the grid box [ckpc/h]. Defaults to 1.05x the largest
            coordinate offset from the group centre, which for a FOF group can
            be several R_200c: the group is a percolated, filamentary object.

        Returns
        -------
        GroupDataProducts
        """
        if not is_known(hi_model):
            raise ValueError(
                f"unknown hi_model {hi_model!r}; choose one of "
                "['GD14_vol', 'GD14_map', 'K13_vol', 'K13_map']"
            )
        spec = resolve(hi_model)
        if compute_hi and spec["kind"] != "vol":
            raise ValueError(
                f"{spec['label']} is a projected model: it builds a face-on map "
                "about one galaxy (using its position, velocity and half-mass "
                "radii), which a whole group does not have. Use a volumetric "
                "model (K13_vol or GD14_vol) for group products."
            )
        if uv_mode not in ("adaptive_fft", "global"):
            raise ValueError(
                f"uv_mode must be 'adaptive_fft' or 'global', got {uv_mode!r}"
            )
        if igm_uv not in ("group", "uvb"):
            raise ValueError(f"igm_uv must be 'group' or 'uvb', got {igm_uv!r}")

        dp = cls(ngrid=ngrid)
        groupID = int(groupID)

        # --- header, model configuration, group record ---
        header = loadHeader(basePath, snapNum)
        h = float(header["HubbleParam"])
        a = float(header["Time"])
        z = float(header["Redshift"])
        cfg = HIModelConfig.build(
            hi_config, h=h, z=z, boxsize_ckpc_per_h=header.get("BoxSize")
        )
        members = load_group_members(basePath, snapNum, groupID)
        group = members["group"]
        nsub = int(group["GroupNsubs"])
        r200 = float(group["Group_R_Crit200"])
        centre = np.asarray(group["GroupPos"], dtype="float64")

        # --- particles and per-cell ownership ---
        gas, gas_start = _load_group_block(
            basePath, snapNum, groupID, 0, _GROUP_GAS_FIELDS, _GROUP_GAS_MDI
        )
        n_gas = int(gas["count"])
        gas_len = int(np.asarray(group["GroupLenType"], dtype="int64")[0])
        if n_gas != gas_len:
            raise RuntimeError(
                f"loaded {n_gas} gas cells for group {groupID} but the catalog "
                f"reports GroupLenType[0] = {gas_len}"
            )
        owner, blocks_gas = _owner_index(
            members["start"][:, 0], members["len_type"][:, 0], gas_start, n_gas, "gas"
        )
        m_gas = np.asarray(gas["Masses"], dtype="float64") * (1e10 / h)  # Msun
        _check_ownership(
            owner,
            blocks_gas,
            members["len_type"][:, 0],
            m_gas,
            np.asarray(members["mass_type"][:, 0], dtype="float64") * (1e10 / h),
            "gas",
            groupID,
        )

        dm, dm_start = _load_group_block(basePath, snapNum, groupID, 1, _DM_FIELDS, None)
        n_dm = int(dm["count"])
        owner_dm, _ = _owner_index(
            members["start"][:, 1], members["len_type"][:, 1], dm_start, n_dm, "dm"
        )
        dm_mass_pp = _dm_particle_mass(basePath, snapNum, h)

        # --- per-cell physics and classification ---
        coords = np.asarray(gas["Coordinates"], dtype="float64")  # ckpc/h, comoving
        vel = np.asarray(gas["Velocities"], dtype="float64")  # km/s * sqrt(a)
        X_H = np.asarray(gas["GFM_Metals"], dtype="float64")
        xe = np.asarray(gas["ElectronAbundance"], dtype="float64")
        nH = rho_code_to_nH(np.asarray(gas["Density"], dtype="float64"), X_H, cfg=cfg)
        T = calc_temperature(np.asarray(gas["InternalEnergy"], dtype="float64"), X_H, xe)

        r_gal = np.full(n_gas, np.nan)
        owned = owner >= 0
        if owned.any():
            r_gal[owned] = np.linalg.norm(
                coords[owned] - members["pos"][owner[owned]], axis=1
            )
        r_grp = np.linalg.norm(coords - centre, axis=1)

        r_ism_sub = ism_radius(
            members["rhalf_type"][:, 0], members["rhalf_type"][:, 4], k_gas, k_star
        )
        if nsub:
            r_ism_cell = np.where(owned, r_ism_sub[np.clip(owner, 0, nsub - 1)], 0.0)
        else:
            r_ism_cell = np.zeros(n_gas)
        class_id = classify_ism_cgm_igm(
            owner, r_gal, r_ism_cell, nH=nH, T=T, nH_max=nH_max, T_max=T_max
        )

        # --- HI / H2 ---
        m_hi = np.zeros(n_gas)
        m_h2 = np.zeros(n_gas)
        if compute_hi:
            m_hi, m_h2 = _hi_group(
                gas, owner, members, cfg, spec["label"], uv_mode, igm_uv, ncpu
            )

        # --- redshift-space positions ---
        dm_mass = np.full(n_dm, dm_mass_pp)
        dm_coords = np.asarray(dm["Coordinates"], dtype="float64")
        if cosmology is None:
            cosmology = Planck15
        # s = x + v_los / (a H(z)), with s and x comoving and v_los physical.
        # Particle velocities are stored as km/s * sqrt(a) while the catalog's
        # group velocity is physical km/s.
        shift = 1.0e3 * h / (a * cosmology.H(z).to_value("km/s/Mpc"))  # ckpc/h per km/s
        v_bulk = float(np.asarray(group["GroupVel"], dtype="float64")[0])
        # Subtract the group's bulk motion first. Keeping it would leave the
        # shift dominated by the (arbitrary) systemic velocity of the group --
        # at z = 0 one km/s is 10 ckpc/h, so a 600 km/s group velocity displaces
        # every particle by 6 Mpc and the RSD map comes out empty. This is also
        # the frame the mock observations are made in: the source sits at its
        # systemic velocity and only internal motions distort it.
        coords_rsd = dm_coords_rsd = None
        if compute_rsd:
            coords_rsd = coords.copy()
            coords_rsd[:, 0] += (vel[:, 0] * np.sqrt(a) - v_bulk) * shift
            dm_coords_rsd = dm_coords.copy()
            dm_coords_rsd[:, 0] += (
                np.asarray(dm["Velocities"], dtype="float64")[:, 0] * np.sqrt(a) - v_bulk
            ) * shift

        # --- grid boxes ---
        # The real-space box follows the particles. Redshift space needs its own,
        # longer line-of-sight axis, because a velocity range maps to a length
        # (10 ckpc/h per km/s at z = 0) rather than to the geometry. The two are
        # kept separate precisely so that the real-space products do not change
        # when RSD is switched on: the classification and every real-space grid
        # are real-space quantities, and RSD is a transform applied afterwards.
        half_yz = 1.05 * max(
            float(np.max(np.abs(coords - centre))),
            float(np.max(np.abs(dm_coords - centre))),
            1e-3,
        )
        half_los = half_yz
        half_los_rsd = None
        if compute_rsd:
            half_los_rsd = max(
                half_yz,
                1.05 * float(np.max(np.abs(coords_rsd[:, 0] - centre[0]))),
                1.05 * float(np.max(np.abs(dm_coords_rsd[:, 0] - centre[0]))),
            )
        if box_halfsize_ckpc is not None:
            pinned = max(float(box_halfsize_ckpc), 1e-3)
            half_yz = half_los = pinned
            if compute_rsd:
                half_los_rsd = pinned
        edges = [
            np.linspace(centre[0] - half_los, centre[0] + half_los, dp.ngrid + 1),
            np.linspace(centre[1] - half_yz, centre[1] + half_yz, dp.ngrid + 1),
            np.linspace(centre[2] - half_yz, centre[2] + half_yz, dp.ngrid + 1),
        ]
        dp.grid_edges = tuple(edges)
        dp.grid_edges_rsd = (
            tuple(
                [
                    np.linspace(
                        centre[0] - half_los_rsd, centre[0] + half_los_rsd, dp.ngrid + 1
                    ),
                    edges[1].copy(),
                    edges[2].copy(),
                ]
            )
            if compute_rsd
            else None
        )

        weights = {"m_gas": m_gas, "m_hi": m_hi, "m_h2": m_h2}
        for name, (key, klass) in _MASS_GRIDS.items():
            rsd = name.endswith("_rsd")
            if rsd and not compute_rsd:
                continue
            grid_edges = dp.grid_edges_rsd if rsd else edges
            if key == "dm_mass":
                c, w = (dm_coords_rsd, dm_mass) if rsd else (dm_coords, dm_mass)
                if klass is not None:
                    sel = owner_dm == klass
                    c, w = c[sel], w[sel]
            else:
                if not compute_hi and key in ("m_hi", "m_h2"):
                    continue
                c, w = (coords_rsd, weights[key]) if rsd else (coords, weights[key])
                if klass is not None:
                    sel = class_id == klass
                    c, w = c[sel], w[sel]
            setattr(dp, name, _grid_sum(c, w, grid_edges))
        dp.t_mw = _grid_mean(coords, m_gas, T, edges)
        dp.nh_mw = _grid_mean(coords, m_gas, nH, edges)

        # --- per-cell output and mass budget ---
        dp.cells = {
            "coords_ckpc_h": coords.astype("float32"),
            "m_gas_msun": m_gas.astype("float32"),
            "m_hi_msun": m_hi.astype("float32"),
            "m_h2_msun": m_h2.astype("float32"),
            "nH_cm3": nH.astype("float32"),
            "T_K": T.astype("float32"),
            "r_gal_ckpc_h": r_gal.astype("float32"),
            "r_grp_over_r200c": (
                (r_grp / r200).astype("float32") if r200 > 0 else r_grp.astype("float32")
            ),
            "owner": owner,
            "class_id": class_id,
        }
        dp.subhalos = members
        dp.metadata.update(
            groupID=groupID,
            snapNum=int(snapNum),
            redshift=z,
            hubble_h=h,
            scale_factor=a,
            hi_model=spec["label"] if compute_hi else "none",
            hi_reference=spec["reference"] if compute_hi else "none",
            uv_mode=uv_mode if compute_hi else "none",
            igm_uv=igm_uv if compute_hi else "none",
            vol_fine_grid=float(cfg.vol_fine_grid),
            grid_unit="ckpc/h",
            mass_unit="Msun",
            M200c_msun=float(group["Group_M_Crit200"]) * (1e10 / h),
            R200c_ckpc_h=r200,
            R200c_kpc=r200 / h / (1.0 + z),
            group_pos_ckpc_h=[float(v) for v in centre],
            group_first_sub=members["first"],
            Nsubs=nsub,
            N_sub_with_gas=int((members["len_type"][:, 0] > 0).sum()),
            n_gas=n_gas,
            n_dm=n_dm,
            ism_k_gas=float(k_gas),
            ism_k_star=float(k_star),
            ism_nH_max=nH_max,
            ism_T_max=T_max,
            box_halfsize_ckpc=[half_los, half_yz, half_yz],
            box_bounds_ckpc_h=[float(v) for ax in edges for v in (ax[0], ax[-1])],
            rsd_velocity_frame="group_bulk_subtracted",
        )
        if dp.grid_edges_rsd is not None:
            dp.metadata["box_halfsize_ckpc_rsd"] = [half_los_rsd, half_yz, half_yz]
            dp.metadata["box_bounds_ckpc_h_rsd"] = [
                float(v)
                for ax in dp.grid_edges_rsd
                for v in (ax[0], ax[-1])
            ]
        dp.metadata.update(_mass_budget(m_gas, m_hi, m_h2, class_id, n_dm, dm_mass_pp, owner_dm))
        # How much of the mass actually lands in the redshift-space box. A pinned
        # `box_halfsize_ckpc` that is too small along the line of sight shows up
        # here rather than silently. The two boxes have different cell volumes,
        # so compare masses, not densities.
        if dp.gas_rsd is not None and dp.gas_real is not None:
            dp.metadata["f_mass_in_rsd_box_gas"] = float(
                dp.gas_rsd.sum() * _cell_volume(dp.grid_edges_rsd)
            ) / max(float(dp.gas_real.sum() * _cell_volume(edges)), 1e-30)
        return dp

    # ------------------------------------------------------------------- I/O
    def save_hdf5(self, filename):
        """Write the grids, per-cell arrays, subhalo table and metadata.

        Parameters
        ----------
        filename : str
            Output HDF5 path.
        """
        import h5py

        with h5py.File(filename, "w") as f:
            for i, axis in enumerate("xyz"):
                f.create_dataset(
                    f"grid_edges_{axis}", data=self.grid_edges[i].astype("float64")
                )
                if self.grid_edges_rsd is not None:
                    f.create_dataset(
                        f"grid_edges_{axis}_rsd",
                        data=self.grid_edges_rsd[i].astype("float64"),
                    )
            for name in _ALL_GRIDS:
                val = getattr(self, name, None)
                if val is not None:
                    f.create_dataset(
                        name,
                        data=np.asarray(val, dtype="float32"),
                        compression="lzf",
                        chunks=(min(32, self.ngrid),) * 3,
                    )
            for name, val in self.cells.items():
                f.create_dataset(
                    f"cells/{name}", data=np.asarray(val), compression="lzf"
                )
            for name, val in self.subhalos.items():
                if name == "group":
                    continue  # a dict of catalog scalars, kept in the metadata
                # Scalars (``first``) come back as length-1 arrays: HDF5 scalar
                # datasets cannot carry the compression the rest of the file uses.
                f.create_dataset(
                    f"subhalos/{name}",
                    data=np.atleast_1d(np.asarray(val)),
                    compression="lzf",
                )
            for key, val in self.metadata.items():
                f.attrs[key] = "none" if val is None else val
            f.attrs["ngrid"] = self.ngrid

    @classmethod
    def load_hdf5(cls, filename):
        """Read back a file written by `save_hdf5`.

        Parameters
        ----------
        filename : str
            Input HDF5 path.

        Returns
        -------
        GroupDataProducts
        """
        import h5py

        with h5py.File(filename, "r") as f:
            dp = cls(ngrid=int(f.attrs["ngrid"]))
            dp.grid_edges = tuple(
                np.asarray(f[f"grid_edges_{axis}"][:]) for axis in "xyz"
            )
            if "grid_edges_x_rsd" in f:
                dp.grid_edges_rsd = tuple(
                    np.asarray(f[f"grid_edges_{axis}_rsd"][:]) for axis in "xyz"
                )
            for name in _ALL_GRIDS:
                if name in f:
                    setattr(dp, name, np.asarray(f[name][:]))
            if "cells" in f:
                dp.cells = {k: np.asarray(f[f"cells/{k}"][:]) for k in f["cells"]}
            if "subhalos" in f:
                dp.subhalos = {
                    k: np.asarray(f[f"subhalos/{k}"][:]) for k in f["subhalos"]
                }
            dp.metadata = dict(f.attrs.items())
        return dp


def _dm_particle_mass(basePath, snapNum, h):
    """Dark-matter particle mass [Msun] from the snapshot's mass table."""
    import h5py

    from .sources._illustris_tools import snapPath

    with h5py.File(snapPath(basePath, snapNum), "r") as f:
        mass_table = f["Header"].attrs.get("MassTable", [0.0] * 6)
    return float(mass_table[1]) * (1e10 / h)


def _mass_budget(m_gas, m_hi, m_h2, class_id, n_dm, dm_mass_pp, owner_dm):
    """Group-level mass totals, overall and per class."""
    budget = {
        "m_gas_total_msun": float(m_gas.sum()),
        "m_hi_total_msun": float(m_hi.sum()),
        "m_h2_total_msun": float(m_h2.sum()),
        "m_dm_total_msun": float(n_dm * dm_mass_pp),
        "m_dm_igm_total_msun": float((owner_dm < 0).sum() * dm_mass_pp),
    }
    for klass, label in CLASS_LABELS.items():
        sel = class_id == klass
        name = label.lower()
        budget[f"m_gas_{name}_msun"] = float(m_gas[sel].sum())
        budget[f"m_hi_{name}_msun"] = float(m_hi[sel].sum())
        budget[f"n_cells_{name}"] = int(sel.sum())
    return budget


def _hi_group(gas, owner, members, cfg, model, uv_mode, igm_uv, ncpu):
    """Per-cell HI and H2 for a whole group.

    Every subhalo is modelled on its own, in its own adaptive FFT box, which is
    the prescription the volumetric models were calibrated under; the fuzz is
    modelled separately, with the group-wide UV field if ``igm_uv='group'``.

    The models return TNG code units (1e10 Msun/h); the returned arrays are in
    Msun, like every other mass in this module.
    """
    unit = 1e10 / cfg.h  # code -> Msun
    n_gas = len(np.asarray(gas["Masses"]))
    m_hi = np.zeros(n_gas)
    m_h2 = np.zeros(n_gas)

    fields = [f for f in _GROUP_GAS_FIELDS if f in gas]
    views = {k: np.asarray(gas[k]) for k in fields}

    lengths = np.asarray(members["len_type"][:, 0], dtype="int64")
    starts = np.asarray(members["start"][:, 0], dtype="int64")
    # Blocks are placed by absolute offset, relative to the start of the first
    # subhalo that has gas -- the same origin `_owner_index` uses, since a
    # group's cells begin with its first subhalo (leading fuzz would sit before
    # that and belong to no block either way).
    first = int(np.flatnonzero(lengths > 0)[0]) if np.any(lengths > 0) else None
    blocks = []
    if first is not None:
        ref = int(starts[first])
        blocks = [
            (int(starts[s]) - ref, int(lengths[s])) for s in np.flatnonzero(lengths > 0)
        ]

    def run(a, n):
        sub = {k: v[a : a + n] for k, v in views.items()}
        return compute_hi_mass(sub, model, cfg, uv_mode=uv_mode)

    results = (
        _hi_parallel(run, blocks, ncpu)
        if ncpu > 1
        else [run(a, n) for a, n in blocks]
    )
    for (a, n), res in zip(blocks, results):
        m_hi[a : a + n] = res.m_hi * unit
        m_h2[a : a + n] = res.m_h2 * unit

    fuzz = owner < 0
    if fuzz.any():
        igm_gas = {k: v[fuzz] for k, v in views.items()}
        if igm_uv == "group":
            res = compute_hi_mass(igm_gas, model, cfg, uv=compute_galaxy_fft_uv(views, cfg)[fuzz])
        else:
            res = compute_hi_mass(igm_gas, model, cfg, uv_mode="global")
        m_hi[fuzz] = res.m_hi * unit
        m_h2[fuzz] = res.m_h2 * unit
    return m_hi, m_h2


_HI_WORKER = {}


def _hi_worker(block):
    """Run one subhalo's HI model; the gas arrays arrive through the global."""
    run = _HI_WORKER.get("run")
    if run is None:
        raise RuntimeError(
            "the per-subhalo HI worker did not inherit the group's gas arrays; "
            "use ncpu=1 (the parallel path needs a fork start method)"
        )
    return run(*block)


def _hi_parallel(run, blocks, ncpu):
    """Map ``run`` over the subhalo blocks in worker processes."""
    from multiprocess import Pool

    _HI_WORKER["run"] = run
    try:
        with Pool(processes=int(ncpu)) as pool:
            return pool.map(_hi_worker, blocks)
    finally:
        _HI_WORKER.clear()


__all__ = [
    "CGM",
    "CLASS_LABELS",
    "GroupDataProducts",
    "IGM",
    "ISM",
    "classify_ism_cgm_igm",
    "ism_radius",
    "load_group_members",
]
