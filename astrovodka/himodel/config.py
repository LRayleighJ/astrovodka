"""
Runtime configuration for the astrovodka HI/H2 models.

One frozen dataclass, passed explicitly. The verified reference implementation
(and `TNGHIcalc`) instead used a module-level singleton reached through
``_units()`` / ``_flags()`` accessors, because its verification strategy was to
run the reference and the new code side by side in one process. astrovodka does
not need that: it already threads ``h``/``redshift`` through every public
signature, and a second global configuration object inside a library is a
liability — two packages calling their own ``configure()`` silently share the
first caller's settings, workers depend on ``fork`` to see it, and nothing can
be reset.

So the accessors became a parameter. The *expressions* are unchanged, operand
order included, which is what the bit-for-bit agreement rests on. Note in
particular that the unit conversions below are written ``x / self.h`` and never
``x * (1.0 / self.h)``: those differ by 1 ulp.

TNG/Arepo code units, as `illustris_python` returns them unchanged:

    Masses            code = 1e10 Msun / h
    Coordinates       code = ckpc / h                     (comoving)
    Density           code = 1e10 Msun h^2 / (ckpc/h)^3   (comoving)
    InternalEnergy    code = 1e10 (km/s)^2                (no h)
    StarFormationRate      = Msun/yr                      (already physical)
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from . import constants as const

# TNG100-1. Only defaults: production paths read HubbleParam from the snapshot
# header and box size from the group catalog header.
HUBBLE_DEFAULT = 0.6774
BOXSIZE_CKPC_PER_H_DEFAULT = 75000.0  # TNG100-1 box = 110.717 cMpc

# Rule out the grossest errors (a forgotten argument, a unit mix-up) without
# pretending to validate cosmology.
_H_RANGE = (0.1, 1.5)
_Z_RANGE = (-0.01, 20.0)

K13_VARIANTS = ("hydro", "hydro_nostar", "2p")
RHO_DM_MODES = ("none", "nfw", "const")
R_SOFT_MAP_MODES = ("dx", "sigma", "cell")


@dataclass(frozen=True)
class HIModelConfig:
    """Cosmology + numerical knobs for one HI-model evaluation.

    Parameters
    ----------
    h : float
        Dimensionless Hubble parameter of the snapshot. Take it from the
        snapshot header's ``HubbleParam``.
    z : float
        Redshift of the snapshot. This is the **simulation** redshift: it sets
        the unit conversions, never the observing geometry.
    boxsize_ckpc_per_h : float
        Periodic box length [ckpc/h], used for the projection wrap guard.
    ngrid : int
        Edge length of the 3D FFT UV grid.
    f_esc : float
        UV escape fraction, Diemer+2018 §4.2 calibration.
    box_expand : float
        Adaptive UV box edge = ``box_expand`` x the coordinate extent.
    ngrid_map : int
        Edge length of the 2D projection grid (``map_r_npix`` upstream).
    map_sub : int
        Deposition supersampling factor.
    map_sigma_bins : int
        Number of log-spaced bins for the per-cell Gaussian kernel width.
    r_soft_factor_map : float
        2D 1/r^2 kernel softening, in units of the pixel size.
    vol_fine_grid : float
        Half-width of the second, higher-resolution 3D UV box, in units of the
        gas half-mass radius (paper Appendix A.2). ``0.0`` disables it, which is
        the reference implementation's behaviour.
    r_soft_map_mode : {'dx', 'sigma', 'cell'}
        Where the 2D kernel's softening length comes from.
    map_clip_outside : bool
        If True, gas cells outside the projection box are clamped to the edge
        pixel instead of dropped (mass conserving).
    k13_variant : {'hydro', 'hydro_nostar', '2p'}
        Which prescription sets ``n_CNM`` in the projected K13 model.
    rho_dm_mode : {'none', 'nfw', 'const'}
        How the projected K13 model estimates the dark-matter mid-plane density.
    """

    h: float = HUBBLE_DEFAULT
    z: float = 0.0
    boxsize_ckpc_per_h: float = BOXSIZE_CKPC_PER_H_DEFAULT

    ngrid: int = const.NGRID
    f_esc: float = const.F_ESC
    box_expand: float = const.BOX_EXPAND

    ngrid_map: int = const.NGRID_MAP
    map_sub: int = const.MAP_SUB
    map_sigma_bins: int = const.MAP_SIGMA_BINS
    r_soft_factor_map: float = const.R_SOFT_FACTOR_MAP

    vol_fine_grid: float = 0.0
    r_soft_map_mode: str = "dx"
    map_clip_outside: bool = False
    k13_variant: str = "2p"
    rho_dm_mode: str = "none"

    def __post_init__(self):
        if not _H_RANGE[0] < self.h < _H_RANGE[1]:
            raise ValueError(f"h={self.h:.6g} is outside the plausible range {_H_RANGE}")
        if not _Z_RANGE[0] < self.z < _Z_RANGE[1]:
            raise ValueError(f"z={self.z:.6g} is outside the plausible range {_Z_RANGE}")
        for name, allowed in (
            ("k13_variant", K13_VARIANTS),
            ("rho_dm_mode", RHO_DM_MODES),
            ("r_soft_map_mode", R_SOFT_MAP_MODES),
        ):
            value = getattr(self, name)
            if value not in allowed:
                raise ValueError(f"{name} must be one of {list(allowed)}, got {value!r}")
        if self.ngrid < 8 or self.ngrid_map < 8:
            raise ValueError("ngrid and ngrid_map must be at least 8")

    @classmethod
    def build(cls, config=None, *, h, z, boxsize_ckpc_per_h=None) -> "HIModelConfig":
        """Resolve a user-supplied config against a snapshot's own header values.

        ``config`` may be None, a dict of field overrides, or an
        ``HIModelConfig``. Whatever it says, ``h``, ``z`` and (when the header
        has one) the box size are taken from the snapshot being read: those
        three describe the data, not a modelling choice, and mixing a
        snapshot's particles with another snapshot's units scales every mass
        and density silently.
        """
        if config is None:
            base = cls()
        elif isinstance(config, cls):
            base = config
        elif isinstance(config, dict):
            base = cls().with_(**config)
        else:
            raise TypeError(
                "the HI model config must be an HIModelConfig, a dict of "
                f"HIModelConfig fields, or None; got {type(config).__name__}"
            )

        overrides = {"h": float(h), "z": float(z)}
        if boxsize_ckpc_per_h:
            overrides["boxsize_ckpc_per_h"] = float(boxsize_ckpc_per_h)
        return base.with_(**overrides)

    def with_(self, **overrides) -> "HIModelConfig":
        """Return a copy with the given fields replaced.

        Unlike the reference implementation's ``configure()``, this never
        inherits state from a previous call: every field not named here keeps
        *this* config's value, and the caller can always get back to the
        documented defaults with ``HIModelConfig(h=..., z=...)``.
        """
        known = set(self.__dataclass_fields__)
        unknown = set(overrides) - known
        if unknown:
            raise TypeError(
                f"unknown HIModelConfig field(s) {sorted(unknown)}; "
                f"available: {sorted(known)}"
            )
        return replace(self, **overrides)

    # ------------------------------------------------------------------ lengths
    def ckpc_to_cmpc(self, x):
        """Raw ckpc/h -> comoving Mpc."""
        return x / self.h / 1e3

    def ckpc_to_physical_kpc(self, x):
        """Raw ckpc/h -> physical kpc."""
        return x / self.h / (1.0 + self.z)

    @property
    def box_cmpc(self):
        """Box length [comoving Mpc]."""
        return self.boxsize_ckpc_per_h / self.h / 1e3

    # ------------------------------------------------------------------ masses
    def mass_code_to_msun(self, m):
        """Raw code (1e10 Msun/h) -> Msun."""
        return m * 1e10 / self.h

    def mass_code_to_legacy(self, m):
        """Raw code (1e10 Msun/h) -> the 1e10 Msun unit used on disk."""
        return m / self.h

    # ----------------------------------------------------------------- density
    def dens_code_to_msun_kpc3(self, rho):
        """Raw code -> physical Msun/kpc^3."""
        return rho * 1e10 * self.h**2 * (1.0 + self.z) ** 3

    def dens_code_to_g_cm3(self, rho):
        """Raw code -> physical g/cm^3."""
        return (
            rho
            * 1e10
            * const.MSUN
            * self.h**2
            / const.kpc_to_cm**3
            * (1.0 + self.z) ** 3
        )

    def summary(self) -> str:
        return (
            f"h={self.h:.5f}, z={self.z:.5f}, box={self.box_cmpc:.4f} cMpc "
            f"({self.boxsize_ckpc_per_h:.1f} ckpc/h)"
        )


__all__ = [
    "HIModelConfig",
    "HUBBLE_DEFAULT",
    "BOXSIZE_CKPC_PER_H_DEFAULT",
    "K13_VARIANTS",
    "RHO_DM_MODES",
    "R_SOFT_MAP_MODES",
]
