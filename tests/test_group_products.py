"""Tests for the group-scale (FOF group) products.

Everything here is synthetic: the membership rules, the class definitions and
the file format are testable without a snapshot, and keeping it that way means
the suite stays runnable without TNG data. The end-to-end checks that *do* need
a snapshot are the assertions of `examples/*/group_cluster_viz.py`.
"""

import numpy as np
import pytest

from astrovodka.group_products import (
    _ALL_GRIDS,
    CGM,
    IGM,
    ISM,
    GroupDataProducts,
    _check_ownership,
    _grid_mean,
    _grid_sum,
    _owner_index,
    classify_ism_cgm_igm,
    ism_radius,
)


class TestIsmRadius:
    def test_is_the_smaller_of_the_two_radii(self):
        r = ism_radius([10.0], [3.0], k_gas=1.0, k_star=2.0)
        assert r[0] == pytest.approx(6.0)

    def test_star_radius_limits_the_central_subhalo(self):
        # The case that motivated the min(): a group's central subhalo has a gas
        # half-mass radius of hundreds of ckpc/h because it holds the hot halo,
        # so a gas-only radius would swallow the group.
        r = ism_radius([263.0], [18.98], k_gas=1.0, k_star=2.0)
        assert r[0] == pytest.approx(2 * 18.98)
        assert r[0] < 263.0

    def test_k_multipliers_scale_each_radius(self):
        assert ism_radius([10.0], [10.0], 0.5, 3.0)[0] == pytest.approx(5.0)
        assert ism_radius([10.0], [10.0], 3.0, 0.5)[0] == pytest.approx(5.0)

    @pytest.mark.parametrize("gas,star,expected", [
        (0.0, 5.0, 10.0),      # a zero gas radius is "no constraint"
        (np.nan, 5.0, 10.0),
        (-1.0, 5.0, 10.0),
        (0.0, 0.0, 0.0),       # neither usable -> no ISM for this subhalo
        (np.nan, np.nan, 0.0),
    ])
    def test_unusable_radii(self, gas, star, expected):
        assert ism_radius([gas], [star], 1.0, 2.0)[0] == pytest.approx(expected)


class TestClassify:
    def test_owner_minus_one_is_igm_whatever_else_says(self):
        # IGM cells are outside every subhalo, so their r_gal is meaningless
        # (NaN in the products) and must not affect the outcome.
        out = classify_ism_cgm_igm(
            owner=np.array([-1, -1]), r_gal=np.array([np.nan, 1e9]), r_ism=10.0
        )
        assert list(out) == [IGM, IGM]

    def test_radius_splits_bound_gas_into_ism_and_cgm(self):
        out = classify_ism_cgm_igm(
            owner=np.array([0, 0]), r_gal=np.array([5.0, 50.0]), r_ism=10.0
        )
        assert list(out) == [ISM, CGM]

    def test_r_ism_may_be_per_cell(self):
        out = classify_ism_cgm_igm(
            owner=np.array([0, 1]), r_gal=np.array([5.0, 5.0]),
            r_ism=np.array([10.0, 1.0]),
        )
        assert list(out) == [ISM, CGM]

    def test_phase_caps_demote_out_of_the_ism(self):
        owner = np.array([0, 0, 0])
        r_gal = np.array([1.0, 1.0, 1.0])
        nH = np.array([1e-3, 1.0, 1e-3])
        T = np.array([1e4, 1e4, 1e8])
        assert list(classify_ism_cgm_igm(owner, r_gal, 10.0)) == [ISM] * 3
        assert list(classify_ism_cgm_igm(owner, r_gal, 10.0, nH=nH, nH_max=0.1)) == [
            ISM, CGM, ISM
        ]
        assert list(classify_ism_cgm_igm(owner, r_gal, 10.0, T=T, T_max=1e6)) == [
            ISM, ISM, CGM
        ]

    def test_dtype_is_int8(self):
        out = classify_ism_cgm_igm(np.array([0]), np.array([1.0]), 10.0)
        assert out.dtype == np.int8

    @pytest.mark.parametrize("kwargs,message", [
        ({"nH_max": 1.0}, "nH"),
        ({"T_max": 1.0}, "T"),
    ])
    def test_cap_without_its_array_raises(self, kwargs, message):
        with pytest.raises(ValueError, match=message):
            classify_ism_cgm_igm(np.array([0]), np.array([1.0]), 10.0, **kwargs)


class TestOwnerIndex:
    def test_places_blocks_and_marks_the_fuzz(self):
        # Two subhalos starting at 0 and 4, of length 4 and 3, then 3 fuzz cells.
        owner, blocks = _owner_index(
            start=np.array([100, 104]), lengths=np.array([4, 3]),
            group_start=100, n_cells=10, label="gas",
        )
        assert list(owner) == [0, 0, 0, 0, 1, 1, 1, -1, -1, -1]
        assert list(blocks) == [0, 4]

    def test_gas_free_subhalos_get_no_cells(self):
        owner, blocks = _owner_index(
            start=np.array([0, 0, 5]), lengths=np.array([5, 0, 5]),
            group_start=0, n_cells=10, label="gas",
        )
        assert list(blocks) == [0, -1, 5]
        assert (owner == 1).sum() == 0

    def test_leading_fuzz_is_left_unassigned(self):
        # A group whose first subhalo starts after the group does: the cells
        # before it belong to no subhalo.
        owner, _ = _owner_index(
            start=np.array([103]), lengths=np.array([4]),
            group_start=100, n_cells=10, label="gas",
        )
        assert list(owner) == [-1, -1, -1, 0, 0, 0, 0, -1, -1, -1]

    def test_overlapping_blocks_raise(self):
        with pytest.raises(RuntimeError, match="overlap"):
            _owner_index(
                start=np.array([0, 5]), lengths=np.array([10, 10]),
                group_start=0, n_cells=20, label="gas",
            )

    def test_block_outside_the_group_raises(self):
        with pytest.raises(RuntimeError, match="disagree"):
            _owner_index(
                start=np.array([0, 100]), lengths=np.array([10, 10]),
                group_start=0, n_cells=20, label="gas",
            )

    def test_a_group_with_no_subhalos_is_all_fuzz(self):
        owner, blocks = _owner_index(
            start=np.zeros(0, dtype="int64"), lengths=np.zeros(0, dtype="int64"),
            group_start=0, n_cells=5, label="gas",
        )
        assert list(owner) == [-1] * 5
        assert len(blocks) == 0


class TestOwnershipCheck:
    def test_matching_masses_pass(self):
        owner, blocks = _owner_index(
            np.array([0, 4]), np.array([4, 4]), 0, 8, "gas"
        )
        masses = np.array([1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0])
        _check_ownership(owner, blocks, np.array([4, 4]), masses,
                         np.array([4.0, 8.0]), "gas", 0)

    def test_mismatched_masses_raise(self):
        # This is the guard against the offsets file and the snapshot layout
        # drifting apart -- the failure is silent otherwise.
        owner, blocks = _owner_index(
            np.array([0, 4]), np.array([4, 4]), 0, 8, "gas"
        )
        masses = np.array([1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0])
        with pytest.raises(RuntimeError, match="catalog masses"):
            _check_ownership(owner, blocks, np.array([4, 4]), masses,
                             np.array([4.0, 3.0]), "gas", 0)

    def test_subhalos_without_a_catalog_mass_are_skipped(self):
        owner, blocks = _owner_index(np.array([0]), np.array([4]), 0, 4, "gas")
        _check_ownership(owner, blocks, np.array([4]), np.ones(4),
                         np.array([0.0]), "gas", 0)


class TestGrids:
    def test_grid_sum_is_a_mass_density(self):
        edges = [np.linspace(0, 2, 3)] * 3  # 2 cells of side 1
        coords = np.array([[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]])
        grid = _grid_sum(coords, np.array([3.0, 4.0]), edges)
        assert grid[0, 0, 0] == pytest.approx(7.0)  # 7 Msun in a unit volume
        assert grid.sum() == pytest.approx(7.0)

    def test_grid_mean_is_mass_weighted(self):
        edges = [np.linspace(0, 2, 3)] * 3
        coords = np.array([[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]])
        weights = np.array([1.0, 3.0])
        values = np.array([10.0, 20.0])
        grid = _grid_mean(coords, weights, values, edges)
        assert grid[0, 0, 0] == pytest.approx((1 * 10 + 3 * 20) / 4)

    def test_empty_cells_are_zero_not_nan(self):
        edges = [np.linspace(0, 2, 3)] * 3
        grid = _grid_mean(np.array([[0.5, 0.5, 0.5]]), np.array([1.0]),
                          np.array([5.0]), edges)
        assert grid[1, 1, 1] == 0.0
        assert np.isfinite(grid).all()


def _synthetic_products(ngrid=8):
    dp = GroupDataProducts(ngrid=ngrid)
    rng = np.random.default_rng(7)
    edges = tuple(np.linspace(0, ngrid, ngrid + 1) for _ in range(3))
    dp.grid_edges = edges
    # Redshift space gets its own, longer line-of-sight axis.
    dp.grid_edges_rsd = (
        np.linspace(0, 3 * ngrid, ngrid + 1), edges[1].copy(), edges[2].copy()
    )
    for name in _ALL_GRIDS:
        setattr(dp, name, rng.random((ngrid,) * 3).astype("float32"))
    dp.cells = {
        "coords_ckpc_h": rng.random((10, 3)).astype("float32"),
        "m_gas_msun": rng.random(10).astype("float32"),
        "owner": np.array([-1, 0, 1] * 3 + [2], dtype="int32"),
        "class_id": np.array([IGM, ISM, CGM] * 3 + [IGM], dtype="int8"),
    }
    dp.subhalos = {
        "first": 42,
        "pos": rng.random((4, 3)),
        "len_type": np.array([[5, 5, 0, 0, 0, 0]] * 4),
    }
    dp.metadata = {
        "groupID": 59,
        "Nsubs": 4,
        "ism_nH_max": None,
        "box_halfsize_ckpc": [1.0, 2.0, 3.0],
        "m_hi_total_msun": 1.5e11,
    }
    return dp


class TestGroupDataProducts:
    def test_fresh_instance_has_no_grids(self):
        dp = GroupDataProducts(ngrid=64)
        assert dp.ngrid == 64
        assert all(getattr(dp, name) is None for name in _ALL_GRIDS)
        assert dp.grid_edges is None
        assert dp.grid_edges_rsd is None
        assert dp.metadata == {}
        assert dp.cells == {}

    def test_hdf5_round_trip(self, tmp_path):
        dp = _synthetic_products()
        path = tmp_path / "group_products.hdf5"
        dp.save_hdf5(path)
        back = GroupDataProducts.load_hdf5(path)

        assert back.ngrid == dp.ngrid
        for axis, (a, b) in enumerate(zip(dp.grid_edges, back.grid_edges)):
            assert np.array_equal(a, b), f"grid_edges axis {axis}"
        for axis, (a, b) in enumerate(zip(dp.grid_edges_rsd, back.grid_edges_rsd)):
            assert np.array_equal(a, b), f"grid_edges_rsd axis {axis}"
        for name in _ALL_GRIDS:
            assert np.array_equal(getattr(dp, name), getattr(back, name)), name
        for name, val in dp.cells.items():
            assert np.array_equal(val, back.cells[name]), name
        assert np.array_equal(dp.subhalos["pos"], back.subhalos["pos"])
        assert float(back.metadata["m_hi_total_msun"]) == dp.metadata["m_hi_total_msun"]

    def test_a_file_without_rsd_grids_loads_with_no_rsd_edges(self, tmp_path):
        dp = GroupDataProducts(ngrid=4)
        dp.grid_edges = tuple(np.linspace(0, 4, 5) for _ in range(3))
        path = tmp_path / "no_rsd.hdf5"
        dp.save_hdf5(path)
        back = GroupDataProducts.load_hdf5(path)
        assert back.grid_edges_rsd is None
        assert back.hi_rsd is None

    def test_scalar_subhalo_entries_survive_as_length_one_arrays(self, tmp_path):
        # h5py cannot compress a scalar dataset, so `first` is stored as an array.
        dp = _synthetic_products()
        path = tmp_path / "group_products.hdf5"
        dp.save_hdf5(path)
        back = GroupDataProducts.load_hdf5(path)
        assert list(np.atleast_1d(back.subhalos["first"])) == [42]

    def test_none_metadata_is_written_as_a_sentinel(self, tmp_path):
        dp = _synthetic_products()
        path = tmp_path / "group_products.hdf5"
        dp.save_hdf5(path)
        back = GroupDataProducts.load_hdf5(path)
        assert back.metadata["ism_nH_max"] == "none"
