"""Tests for SubhaloDataProducts, spectral_unit, and UV grid exposure."""

import numpy as np
import astropy.units as U
import os
import pytest
from astrovodka import DataCube, SubhaloDataProducts
from astrovodka.himodel import HIModelConfig, calc_UMW_global_Faucher11
from astrovodka.data_products import _rsd_los_halfwidth
from astrovodka.himodel.uvfield import compute_galaxy_fft_uv


class TestSpectralUnit:
    """Test the spectral_unit parameter in DataCube."""

    def test_default_auto_detects_velocity(self):
        """With no spectral_unit the axis adopts the units of channel_width.

        This matches MARTINI, and is what a velocity channel_width must give.
        """
        dc = DataCube(
            n_px_x=32, n_px_y=32, n_channels=16,
            px_size=10 * U.arcsec,
            channel_width=20 * U.km / U.s,
            spectral_centre=0 * U.km / U.s,
        )
        assert dc._freq_channel_mode is False
        assert dc.wcs.wcs.ctype[2] == "VRAD"
        assert U.get_physical_type(dc.channel_width) == "velocity"

    def test_default_auto_detects_frequency(self):
        """A frequency channel_width with no spectral_unit gives a frequency axis."""
        dc = DataCube(
            n_px_x=32, n_px_y=32, n_channels=16,
            px_size=10 * U.arcsec,
            channel_width=100 * U.kHz,
            spectral_centre=1.42 * U.GHz,
        )
        assert dc._freq_channel_mode is True
        assert dc.wcs.wcs.ctype[2] == "FREQ"

    def test_velocity_to_frequency_width_is_differential(self):
        """A velocity width must convert to a small frequency *offset*, not f_rest.

        Regression test: applying astropy's ``doppler_radio`` equivalency directly to
        a width returns the absolute frequency at that velocity (~1420 MHz for
        20 km/s, i.e. f_rest) rather than the width (~94.76 kHz), which shrank the
        flux entering the cube by ~7 orders of magnitude.
        """
        from astrovodka.cosmology_utils import HIfreq
        from astropy import constants as const

        dc = DataCube(
            n_px_x=32, n_px_y=32, n_channels=16,
            px_size=10 * U.arcsec,
            channel_width=20 * U.km / U.s,
            spectral_centre=0 * U.km / U.s,
            spectral_unit="MHz",
        )
        expected = (HIfreq * 20 * U.km / U.s / const.c).to(U.Hz)
        assert np.isclose(dc.channel_width.to_value(U.Hz), expected.to_value(U.Hz))
        # and nowhere near the rest frequency
        assert dc.channel_width < 1 * U.MHz

    def test_frequency_to_velocity_width_is_differential(self):
        """Inverse of the above: a frequency width maps to a small velocity width."""
        from astrovodka.cosmology_utils import HIfreq
        from astropy import constants as const

        width = 100 * U.kHz
        dc = DataCube(
            n_px_x=32, n_px_y=32, n_channels=16,
            px_size=10 * U.arcsec,
            channel_width=width,
            spectral_centre=1.42 * U.GHz,
            spectral_unit="km/s",
        )
        expected = (const.c * width / HIfreq).to(U.km / U.s)
        assert np.isclose(
            dc.channel_width.to_value(U.km / U.s), expected.to_value(U.km / U.s)
        )
        assert dc.channel_width < 100 * U.km / U.s

    def test_frequency_axis_span_is_sane(self):
        """The whole band must span n_channels * width, not n_channels * f_rest."""
        dc = DataCube(
            n_px_x=16, n_px_y=16, n_channels=32,
            px_size=10 * U.arcsec,
            channel_width=20 * U.km / U.s,
            spectral_centre=0 * U.km / U.s,
            spectral_unit="MHz",
        )
        span = np.abs(
            (dc.frequency_channel_edges[-1] - dc.frequency_channel_edges[0]).to(U.MHz)
        )
        assert span < 5 * U.MHz  # ~0.6 MHz wide in reality; f_rest bug gave ~45 GHz
        assert span > 0.1 * U.MHz

    def test_invalid_channel_width_units(self):
        """channel_width must have frequency or velocity units."""
        for bad in (1.0 * U.dimensionless_unscaled, 1.0 * U.m):
            with pytest.raises(ValueError, match="Channel width"):
                DataCube(
                    n_px_x=8, n_px_y=8, n_channels=8,
                    px_size=10 * U.arcsec,
                    channel_width=bad,
                )

    def test_explicit_kms(self):
        """spectral_unit='km/s' should keep velocity axis."""
        dc = DataCube(
            n_px_x=32, n_px_y=32, n_channels=16,
            px_size=10 * U.arcsec,
            channel_width=20 * U.km / U.s,
            spectral_centre=0 * U.km / U.s,
            spectral_unit="km/s",
        )
        assert dc._freq_channel_mode is False
        assert dc.wcs.wcs.ctype[2] == "VRAD"

    def test_auto_detect_velocity(self):
        """spectral_unit=None with velocity input should keep velocity."""
        dc = DataCube(
            n_px_x=32, n_px_y=32, n_channels=16,
            px_size=10 * U.arcsec,
            channel_width=20 * U.km / U.s,
            spectral_centre=0 * U.km / U.s,
            spectral_unit=None,
        )
        assert dc._freq_channel_mode is False

    def test_auto_detect_frequency(self):
        """spectral_unit=None with frequency input should keep frequency."""
        dc = DataCube(
            n_px_x=32, n_px_y=32, n_channels=16,
            px_size=10 * U.arcsec,
            channel_width=100 * U.kHz,
            spectral_centre=1.42 * U.GHz,
            spectral_unit=None,
        )
        assert dc._freq_channel_mode is True

    def test_frequency_input_to_velocity(self):
        """Frequency input + spectral_unit='km/s' should convert."""
        dc = DataCube(
            n_px_x=32, n_px_y=32, n_channels=16,
            px_size=10 * U.arcsec,
            channel_width=100 * U.kHz,
            spectral_centre=1.42 * U.GHz,
            spectral_unit="km/s",
        )
        assert dc._freq_channel_mode is False
        assert dc.wcs.wcs.ctype[2] == "VRAD"

    def test_invalid_spectral_unit(self):
        """Invalid spectral_unit should raise ValueError."""
        with pytest.raises(ValueError, match="spectral_unit"):
            DataCube(
                n_px_x=32, n_px_y=32, n_channels=16,
                px_size=10 * U.arcsec,
                channel_width=20 * U.km / U.s,
                spectral_unit="invalid",
            )


class TestUVGridExposure:
    """Test that compute_galaxy_fft_uv can return the grid."""

    CFG = HIModelConfig(h=0.6774, z=0.0, ngrid=32)

    def test_return_grid(self):
        """return_grid=True should return (per_particle, grid, edges)."""
        # Create synthetic gas data
        n = 100
        np.random.seed(42)
        gas_data = {
            "Coordinates": np.random.randn(n, 3) * 10 + 50,  # ~50 ckpc/h
            "StarFormationRate": np.zeros(n),
        }
        gas_data["StarFormationRate"][:10] = 1.0  # 10 SF particles

        # Without return_grid
        result = compute_galaxy_fft_uv(gas_data, self.CFG, return_grid=False)
        assert result.shape == (n,)
        assert result.min() > 0  # UVB floor

        # With return_grid
        per_particle, grid, edges = compute_galaxy_fft_uv(
            gas_data, self.CFG, return_grid=True
        )
        assert per_particle.shape == (n,)
        assert grid.shape == (32, 32, 32)
        assert len(edges) == 3
        assert all(len(e) == 33 for e in edges)
        # The three edge arrays are centred on their own axis, not shared: the
        # old model reused one axis' edges for all three, which silently zeroed
        # the interpolated field.
        assert not np.array_equal(edges[0], edges[1])

    def test_uvb_only_matches_the_global_background(self):
        """uv_mode='global' is the uniform UVB, and still exposes a grid."""
        n = 50
        np.random.seed(42)
        gas_data = {
            "Coordinates": np.random.randn(n, 3) * 5 + 300,
            "StarFormationRate": np.ones(n),  # plenty of SF, but ignored
        }
        per_particle, grid, _ = compute_galaxy_fft_uv(
            gas_data, self.CFG, return_grid=True, uvb_only=True
        )
        uvb = calc_UMW_global_Faucher11(self.CFG.z)
        assert np.allclose(per_particle, uvb)
        assert np.allclose(grid, uvb)

    def test_no_star_formation_still_gives_the_uvb_floor(self):
        """With no SF particles the field is just the background."""
        n = 50
        np.random.seed(42)
        gas_data = {
            "Coordinates": np.random.randn(n, 3) * 5,
            "StarFormationRate": np.zeros(n),  # No SF
        }
        per_particle, grid, edges = compute_galaxy_fft_uv(
            gas_data, self.CFG, return_grid=True
        )
        assert grid.shape == (32, 32, 32)
        assert grid.min() >= 0
        assert np.allclose(per_particle, calc_UMW_global_Faucher11(self.CFG.z))


class TestSubhaloDataProducts:
    """Test SubhaloDataProducts HDF5 save/load and structure."""

    def test_init(self):
        """Test basic initialization."""
        dp = SubhaloDataProducts(ngrid=64)
        assert dp.ngrid == 64
        assert dp.dm_real is None
        assert dp.grid_edges_rsd is None
        assert dp.metadata == {}

    def test_hdf5_roundtrip(self):
        """Save and load should preserve all arrays."""
        ngrid = 32
        dp = SubhaloDataProducts(ngrid=ngrid)

        # Fill with synthetic data
        shape = (ngrid, ngrid, ngrid)
        dp.uv_field = np.random.rand(*shape).astype("float32")
        dp.dm_real = np.random.rand(*shape).astype("float32")
        dp.dm_rsd = np.random.rand(*shape).astype("float32")
        dp.hi_real = np.random.rand(*shape).astype("float32")
        dp.hi_rsd = np.random.rand(*shape).astype("float32")
        dp.h2_real = np.random.rand(*shape).astype("float32")
        dp.h2_rsd = np.random.rand(*shape).astype("float32")
        dp.grid_edges = tuple(
            np.linspace(-10, 10, ngrid + 1) for _ in range(3)
        )
        # the redshift-space box has its own, longer line-of-sight axis
        dp.grid_edges_rsd = (
            np.linspace(-100, 100, ngrid + 1),
            np.linspace(-10, 10, ngrid + 1),
            np.linspace(-10, 10, ngrid + 1),
        )
        dp.metadata = {"redshift": 0.0, "hubble_h": 0.6774}

        test_file = "/tmp/test_dp_roundtrip.hdf5"
        try:
            dp.save_hdf5(test_file)
            dp2 = SubhaloDataProducts.load_hdf5(test_file)

            for name in ["uv_field", "dm_real", "dm_rsd", "hi_real",
                         "hi_rsd", "h2_real", "h2_rsd"]:
                assert np.allclose(getattr(dp, name), getattr(dp2, name),
                                   rtol=1e-5), f"{name} mismatch"

            assert dp2.ngrid == ngrid
            assert len(dp2.grid_edges) == 3
            for a, b in zip(dp.grid_edges_rsd, dp2.grid_edges_rsd):
                assert np.array_equal(a, b)
        finally:
            if os.path.exists(test_file):
                os.remove(test_file)

    def test_hdf5_without_rsd_edges_loads_as_none(self):
        """A file written before the RSD box existed still loads."""
        ngrid = 8
        dp = SubhaloDataProducts(ngrid=ngrid)
        dp.grid_edges = tuple(np.linspace(-10, 10, ngrid + 1) for _ in range(3))
        test_file = "/tmp/test_dp_no_rsd.hdf5"
        try:
            dp.save_hdf5(test_file)
            dp2 = SubhaloDataProducts.load_hdf5(test_file)
            assert dp2.grid_edges_rsd is None
        finally:
            if os.path.exists(test_file):
                os.remove(test_file)


class TestRsdLosHalfwidth:
    """The line-of-sight box sizing, which is where the RSD bug lived.

    Before this helper existed, the redshift-space grids shared the real-space
    box: at z = 0 one km/s is 10 ckpc/h, so a galaxy's peculiar velocities put
    almost every particle outside it and the grids kept 0.5% of their mass.
    """

    def test_rsd_box_holds_every_shifted_particle(self):
        real = np.array([-3.0, 4.0, 0.5])
        shifted = np.array([-2000.0, 1500.0, 12.0])
        _, rsd_half = _rsd_los_halfwidth(real, shifted, centre=0.0)
        assert rsd_half >= np.max(np.abs(shifted))
        assert rsd_half >= np.max(np.abs(real))

    def test_rsd_box_honours_the_centre(self):
        # A box centred at 100 must reach the particles at 90 and 140.
        _, rsd_half = _rsd_los_halfwidth(
            np.array([95.0, 105.0]), np.array([90.0, 140.0]), centre=100.0
        )
        assert rsd_half >= 40.0

    def test_real_extent_is_not_stretched_by_the_shift(self):
        # The whole point of returning two values: the real-space grids keep
        # their own, compact box.
        real = np.array([-3.0, 4.0])
        real_half, _ = _rsd_los_halfwidth(real, np.array([-3000.0]), centre=0.0)
        assert real_half == pytest.approx(1.05 * 4.0)

    def test_always_uses_the_larger_of_the_two(self):
        # Shifted inside the real extent: the real extent wins.
        real_half, rsd_half = _rsd_los_halfwidth(
            np.array([-100.0]), np.array([1.0]), centre=0.0
        )
        assert rsd_half == pytest.approx(real_half)
