"""Tests for the four selectable HI/H2 models (astrovodka.himodel).

These are analytic/self-contained: no TNG data is required. The stronger
evidence -- bit-for-bit agreement with the reference implementation on real
snapshots -- lives in `docs/en/verification.md`, which records the chain.
"""

import numpy as np
import pytest

from astrovodka.cosmology_utils import rho_code_to_nH, rho_code_to_physical
from astrovodka.himodel import (
    DEFAULT_METHOD,
    MAP_METHODS,
    METHODS,
    VOL_METHODS,
    HIModelConfig,
    calc_UMW_global_Faucher11,
    compute_hi_mass,
    is_known,
    resolve,
)
from astrovodka.himodel.galaxy import compute_HI_mass_particle
from astrovodka.himodel.models import f_mol_GD14_map, f_mol_K13_map
from astrovodka.himodel.physics import f_mol_GD14_vol, f_mol_K13_vol, get_f_neutral_H
from astrovodka.himodel.uvfield import compute_galaxy_fft_uv

CFG = HIModelConfig(h=0.6774, z=0.0)


# --------------------------------------------------------------------------
# synthetic galaxy used by the driver tests
# --------------------------------------------------------------------------
def make_gas(n=400, seed=42, offset=(0.0, 0.0, 0.0)):
    """A small blob of gas cells in raw ckpc/h, deliberately off-centre."""
    rng = np.random.default_rng(seed)
    centre = np.array([12000.0, 8000.0, 20000.0]) + np.asarray(offset)
    return dict(
        count=n,
        Coordinates=centre + rng.normal(0, 250.0, (n, 3)),
        Velocities=rng.normal(0, 120.0, (n, 3)),
        Masses=np.full(n, 1.4e6),
        Density=10 ** rng.uniform(4, 8, n),
        InternalEnergy=10 ** rng.uniform(0.5, 2.5, n),
        StarFormationRate=np.where(rng.random(n) < 0.3, 1e-2, 0.0),
        NeutralHydrogenAbundance=rng.uniform(0, 1, n),
        GFM_Metallicity=10 ** rng.uniform(-4, -1.2, n),
        GFM_Metals=rng.uniform(0.6, 0.9, n),  # 1-D, as mdi=0 gives
    )


def make_stars(n=200, seed=43):
    rng = np.random.default_rng(seed)
    centre = np.array([12000.0, 8000.0, 20000.0])
    return dict(
        count=n,
        Coordinates=centre + rng.normal(0, 150.0, (n, 3)),
        Velocities=rng.normal(0, 150.0, (n, 3)),
        Masses=np.full(n, 5e5),
        GFM_StellarFormationTime=rng.uniform(0.1, 1.0, n),
    )


SUB = dict(
    pos=np.array([12000.0, 8000.0, 20000.0]),
    vel=np.zeros(3),
    rhalf_gas=200.0,
    rhalf_star=150.0,
)


def test_uvb_positive_and_increasing():
    for z in [0.0, 0.5, 1.0, 2.0, 5.0]:
        assert calc_UMW_global_Faucher11(z) > 0.0
    assert calc_UMW_global_Faucher11(1.0) > calc_UMW_global_Faucher11(0.0)
    assert calc_UMW_global_Faucher11(2.0) > calc_UMW_global_Faucher11(0.0)


def test_uvb_at_z0():
    """At z=0 the Faucher-Giguere+2011 field is ~1.4e-3 of the Draine field."""
    assert 1e-4 < calc_UMW_global_Faucher11(0.0) < 1e-2


def test_f_neutral_star_forming_cells_are_fully_neutral():
    nH = np.array([0.5, 0.3, 0.9, 0.1])
    sfr = np.array([0.0, 1e-3, 0.0, 5.0])
    out = get_f_neutral_H(nH, sfr)
    assert out[1] == 1.0
    assert out[3] == 1.0
    assert out[0] == 0.5
    assert out[2] == 0.9


def test_f_neutral_does_not_mutate_input():
    nH = np.array([0.5, 0.3])
    sfr = np.array([0.0, 1.0])
    get_f_neutral_H(nH, sfr)
    assert np.array_equal(nH, np.array([0.5, 0.3]))


def test_density_conversion_redshift_scaling():
    """rho ~ (1+z)^3 for a comoving density, pinned to better than 1%."""
    rho = np.array([1e6, 1e7])
    r0 = rho_code_to_physical(rho, 0.0, 0.6774)
    r1 = rho_code_to_physical(rho, 1.0, 0.6774)
    assert np.allclose(r1 / r0, 8.0, rtol=0.01)
    assert np.all(rho_code_to_nH(rho, 0.0, 0.6774) > 0)


# --------------------------------------------------------------------------
# the registry: keyword selection
# --------------------------------------------------------------------------
def test_four_methods_are_registered():
    assert METHODS == ("GD14_vol", "GD14_map", "K13_vol", "K13_map")
    assert set(VOL_METHODS) == {"GD14_vol", "K13_vol"}
    assert set(MAP_METHODS) == {"GD14_map", "K13_map"}
    assert DEFAULT_METHOD in METHODS
    assert all(is_known(m) for m in METHODS)


@pytest.mark.parametrize(
    "written,expected",
    [
        ("K13_map", "K13_map"),
        ("k13_map", "K13_map"),
        ("  GD14_VOL ", "GD14_vol"),
        ("K13", "K13_vol"),  # bare model name -> volume variant
        ("gd14", "GD14_vol"),
        (("GD14", "map"), "GD14_map"),
    ],
)
def test_method_name_resolution(written, expected):
    assert resolve(written)["label"] == expected


@pytest.mark.parametrize("bogus", ["K13_volume", "", "GD14_map_extra", None, 3])
def test_unknown_method_raises(bogus):
    with pytest.raises((ValueError, TypeError)):
        resolve(bogus)
    assert not is_known(bogus)


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------
def test_build_overrides_h_and_z_from_the_snapshot():
    """h and z describe the data, so they win over the caller's config."""
    cfg = HIModelConfig.build({"ngrid": 64, "h": 0.5, "z": 3.0}, h=0.6774, z=0.0)
    assert cfg.ngrid == 64  # the modelling choice survives
    assert cfg.h == pytest.approx(0.6774)  # the data wins
    assert cfg.z == pytest.approx(0.0)


def test_build_accepts_none_dict_and_config():
    assert HIModelConfig.build(None, h=0.7, z=0.0).h == pytest.approx(0.7)
    assert HIModelConfig.build({"map_sub": 2}, h=0.7, z=0.0).map_sub == 2
    base = HIModelConfig(k13_variant="hydro")
    assert HIModelConfig.build(base, h=0.7, z=0.0).k13_variant == "hydro"


def test_build_rejects_unknown_field():
    with pytest.raises(TypeError):
        HIModelConfig.build({"not_a_field": 1}, h=0.7, z=0.0)
    with pytest.raises(TypeError):
        HIModelConfig.build(42, h=0.7, z=0.0)


def test_config_validates_ranges_and_choices():
    with pytest.raises(ValueError):
        HIModelConfig(h=99.0)
    with pytest.raises(ValueError):
        HIModelConfig(z=1e4)
    with pytest.raises(ValueError):
        HIModelConfig(k13_variant="nonsense")
    with pytest.raises(ValueError):
        HIModelConfig(r_soft_map_mode="nonsense")


# --------------------------------------------------------------------------
# per-model behaviour
# --------------------------------------------------------------------------
@pytest.mark.parametrize("model", METHODS)
def test_every_model_conserves_hydrogen_mass(model):
    """m_HI + m_H2 must never exceed the cell's hydrogen mass."""
    gas = make_gas()
    kw = {"stars": make_stars(), "sub": SUB} if resolve(model)["kind"] == "map" else {}
    res = compute_hi_mass(gas, model, CFG, **kw)

    m_h = gas["Masses"] * gas["GFM_Metals"]
    assert np.all(res.m_neutral >= 0.0)
    assert np.all(res.m_neutral <= m_h * (1.0 + 1e-10))


@pytest.mark.parametrize("model", METHODS)
def test_every_model_stays_in_range(model):
    gas = make_gas()
    kw = {"stars": make_stars(), "sub": SUB} if resolve(model)["kind"] == "map" else {}
    res = compute_hi_mass(gas, model, CFG, **kw)
    f = res.f_mol
    assert np.all(np.isfinite(f))
    assert np.all((f >= 0.0) & (f <= 1.0))


@pytest.mark.parametrize("model", METHODS)
def test_every_model_returns_full_length_arrays(model):
    """The result lines up with the caller's gas arrays, so it can feed SPHSource."""
    gas = make_gas(n=137)
    kw = {"stars": make_stars(), "sub": SUB} if resolve(model)["kind"] == "map" else {}
    res = compute_hi_mass(gas, model, CFG, **kw)
    assert res.m_hi.shape == (137,)
    assert res.m_h2.shape == (137,)
    assert res.f_mol.shape == (137,)


def test_the_four_methods_actually_differ():
    """Selecting a model must select a model, not a no-op."""
    gas, stars = make_gas(), make_stars()
    totals = {}
    for model in METHODS:
        kw = {"stars": stars, "sub": SUB} if resolve(model)["kind"] == "map" else {}
        totals[model] = compute_hi_mass(gas, model, CFG, **kw).m_hi_msun
    assert len(set(round(v, 6) for v in totals.values())) == len(METHODS)


def test_map_models_need_stars_and_sub():
    gas = make_gas()
    for model in MAP_METHODS:
        with pytest.raises(ValueError, match="projects the galaxy"):
            compute_hi_mass(gas, model, CFG)


def test_volume_models_ignore_the_projection_requirement():
    gas = make_gas()
    for model in VOL_METHODS:
        assert compute_hi_mass(gas, model, CFG) is not None


def test_global_uv_mode_is_a_uniform_background():
    """uv_mode='global' must give every cell the same U_MW = the UVB."""
    gas = make_gas()
    res = compute_hi_mass(gas, "K13_vol", CFG, uv_mode="global")
    uvb = calc_UMW_global_Faucher11(CFG.z)
    assert np.allclose(res.u_mw, uvb)


def test_unknown_uv_mode_raises():
    with pytest.raises(ValueError, match="uv_mode"):
        compute_hi_mass(make_gas(), "K13_vol", CFG, uv_mode="adaptive")


def test_adaptive_uv_is_brighter_than_the_background_near_a_starburst():
    gas = make_gas()
    res = compute_hi_mass(gas, "K13_vol", CFG, uv_mode="adaptive_fft")
    assert res.u_mw.max() > calc_UMW_global_Faucher11(CFG.z)


# --------------------------------------------------------------------------
# the regression that motivated replacing the old model
# --------------------------------------------------------------------------
def test_fft_uv_survives_an_off_centre_source():
    """The FFT UV field must not collapse to zero, however the source is placed.

    The model astrovodka shipped before v0.3.0 built its interpolator with one
    axis' edge array reused for all three. Because each axis has its own centre,
    almost every particle then fell outside the interpolator's bounds and
    ``fill_value=0.0`` returned zero -- so the "adaptive FFT" field was silently
    just the uniform background. The self-test in the reference package pins
    this as ``max(U_MW) == 170.078`` before the UVB floor is applied; here we
    check the weaker, placement-independent property: the field responds to
    where the source is.
    """
    base = compute_galaxy_fft_uv(make_gas(), CFG, uvb_only=False)
    uvb = calc_UMW_global_Faucher11(CFG.z)
    assert base.max() > uvb, "FFT UV collapsed to the background"

    again = compute_galaxy_fft_uv(make_gas(), CFG, uvb_only=False)
    assert np.array_equal(base, again), "the FFT UV field is not deterministic"

    # shifting the whole galaxy moves it relative to its own adaptive box: the
    # per-particle values must change (a broken interpolator returns the
    # background regardless of position)
    shifted = make_gas(offset=(500.0, -400.0, 900.0))
    assert not np.array_equal(
        compute_galaxy_fft_uv(shifted, CFG, uvb_only=False), base
    )


def test_fft_uv_grid_matches_the_interpolated_values():
    """return_grid must hand back the field the per-particle values came from."""
    gas = make_gas()
    per_particle, grid, edges = compute_galaxy_fft_uv(gas, CFG, return_grid=True)
    assert grid.shape == (CFG.ngrid,) * 3
    assert len(edges) == 3 and all(len(e) == CFG.ngrid + 1 for e in edges)
    # grid is the raw field, so it need not include the UVB floor the particles get
    assert per_particle.shape == (gas["count"],)
    assert np.all(per_particle >= calc_UMW_global_Faucher11(CFG.z))


# --------------------------------------------------------------------------
# the underlying physics, directly
# --------------------------------------------------------------------------
def _vol_inputs(n=200, seed=1, Z=0.0127):
    rng = np.random.default_rng(seed)
    return dict(
        rho_code=1e6 * np.ones(n),
        u_code=1.0 * np.ones(n),
        Z=np.full(n, Z),
        nH_code=np.full(n, 0.9),
        sfr_code=np.zeros(n),
        X_H=np.full(n, 0.76),
        U_MW=np.full(n, 0.001),
        m_gas=1.4e6 * np.ones(n),
    )


@pytest.mark.parametrize("func", [f_mol_K13_vol, f_mol_GD14_vol])
def test_volume_models_return_fractions(func):
    d = _vol_inputs()
    f_mol, f_HI, f_H2 = func(
        d["rho_code"], d["u_code"], d["Z"], d["nH_code"], d["sfr_code"],
        d["X_H"], d["U_MW"], CFG, mgas=d["m_gas"],
    )
    for f in (f_mol, f_HI, f_H2):
        assert np.all((f >= 0.0) & (f <= 1.0))
    assert np.all(f_HI + f_H2 <= 1.0 + 1e-10)
    assert np.allclose(f_HI + f_H2, 0.9 * np.ones_like(f_HI))  # f_neutral for SFR=0


@pytest.mark.parametrize("func", [f_mol_K13_vol, f_mol_GD14_vol])
def test_more_metallicity_means_more_molecules(func):
    lo = _vol_inputs(Z=0.001)
    hi = _vol_inputs(Z=0.025)
    kw = dict(mgas=lo["m_gas"])
    f_lo = func(lo["rho_code"], lo["u_code"], lo["Z"], lo["nH_code"],
                lo["sfr_code"], lo["X_H"], lo["U_MW"], CFG, **kw)[2]
    f_hi = func(hi["rho_code"], hi["u_code"], hi["Z"], hi["nH_code"],
                hi["sfr_code"], hi["X_H"], hi["U_MW"], CFG, **kw)[2]
    assert f_hi.mean() >= f_lo.mean()


@pytest.mark.parametrize("func", [f_mol_K13_vol, f_mol_GD14_vol])
def test_more_uv_means_fewer_molecules(func):
    dim = _vol_inputs()
    bright = _vol_inputs()
    kw = dict(mgas=dim["m_gas"])
    f_dim = func(dim["rho_code"], dim["u_code"], dim["Z"], dim["nH_code"],
                 dim["sfr_code"], dim["X_H"], dim["U_MW"], CFG, **kw)[2]
    f_bright = func(bright["rho_code"], bright["u_code"], bright["Z"],
                    bright["nH_code"], bright["sfr_code"], bright["X_H"],
                    np.full(bright["X_H"].size, 100.0), CFG, **kw)[2]
    assert f_bright.mean() <= f_dim.mean()


def test_gd14_map_is_the_resolution_free_limit():
    """Large pixels push D_* to its 0.17 floor, so the map answer saturates."""
    kwargs = dict(
        Sigma_n=np.full(4, 1e7), Z_map=np.full(4, 0.0127), U_MW=np.full(4, 0.01)
    )
    coarse = f_mol_GD14_map(pix_pc=1200.0, **kwargs)
    finer = f_mol_GD14_map(pix_pc=5000.0, **kwargs)
    assert np.allclose(coarse, finer, rtol=1e-4)


def test_k13_map_variants_differ_and_validate():
    kwargs = dict(
        Sigma_n=np.full(4, 1e7),
        Z_map=np.full(4, 0.0127),
        U_MW=np.full(4, 0.01),
        Sigma_star=np.full(4, 1e8),
        sigma_star=np.full(4, 30.0),
    )
    two_p = f_mol_K13_map(cfg=CFG.with_(k13_variant="2p"), **kwargs)
    hydro = f_mol_K13_map(cfg=CFG.with_(k13_variant="hydro"), **kwargs)
    assert not np.allclose(two_p, hydro)
    with pytest.raises(ValueError):
        f_mol_K13_map(variant="nonsense", cfg=CFG, **kwargs)


def test_compute_HI_mass_particle_matches_the_driver():
    """The lower-level entry point and the keyword entry point must agree."""
    gas = make_gas()
    cfg = CFG
    uv = compute_galaxy_fft_uv(gas, cfg)
    direct_hi, direct_h2 = compute_HI_mass_particle(gas, uv, "K13", cfg)
    res = compute_hi_mass(gas, "K13_vol", cfg)
    assert np.array_equal(direct_hi, res.m_hi)
    assert np.array_equal(direct_h2, res.m_h2)


def test_2d_gfm_metals_is_accepted():
    """GFM_Metals arrives 1-D via mdi=0, but the full 2-D array is valid too."""
    gas = make_gas()
    a = compute_hi_mass(gas, "K13_vol", CFG)
    gas2 = dict(gas)
    gas2["GFM_Metals"] = np.repeat(gas["GFM_Metals"][:, None], 10, axis=1)
    b = compute_hi_mass(gas2, "K13_vol", CFG)
    assert np.array_equal(a.m_hi, b.m_hi)


def test_result_reports_msun_totals_consistently():
    res = compute_hi_mass(make_gas(), "K13_vol", CFG)
    assert res.m_hi_msun == pytest.approx(np.sum(res.m_hi) * 1e10 / CFG.h)
    assert res.m_neutral_msun == pytest.approx(res.m_hi_msun + res.m_h2_msun)
    assert repr(res).startswith("HiMassResult(K13_vol")
