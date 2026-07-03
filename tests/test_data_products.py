"""Tests for SubhaloDataProducts, spectral_unit, and UV grid exposure."""

import numpy as np
import astropy.units as U
import os
import pytest
from astrovodka import DataCube, SubhaloDataProducts
from astrovodka.sources._k13_model import compute_galaxy_fft_uv


class TestSpectralUnit:
    """Test the spectral_unit parameter in DataCube."""

    def test_default_mhz(self):
        """Default spectral_unit should produce FREQ axis in Hz."""
        dc = DataCube(
            n_px_x=32, n_px_y=32, n_channels=16,
            px_size=10 * U.arcsec,
            channel_width=20 * U.km / U.s,
            spectral_centre=0 * U.km / U.s,
        )
        assert dc._freq_channel_mode is True
        assert dc.wcs.wcs.ctype[2] == "FREQ"
        assert U.get_physical_type(dc.channel_width) == "frequency"

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

        # Without return_grid (backward compatible)
        result = compute_galaxy_fft_uv(gas_data, ngrid=32, return_grid=False)
        assert result.shape == (n,)
        assert result.min() > 0  # UVB floor

        # With return_grid
        per_particle, grid, edges = compute_galaxy_fft_uv(
            gas_data, ngrid=32, return_grid=True
        )
        assert per_particle.shape == (n,)
        assert grid.shape == (32, 32, 32)
        assert len(edges) == 3
        assert all(len(e) == 33 for e in edges)

    def test_global_uvb_fallback(self):
        """When no SF particles, grid should still be created with UVB."""
        n = 50
        np.random.seed(42)
        gas_data = {
            "Coordinates": np.random.randn(n, 3) * 5,
            "StarFormationRate": np.zeros(n),  # No SF
        }
        per_particle, grid, edges = compute_galaxy_fft_uv(
            gas_data, ngrid=32, return_grid=True
        )
        assert grid.shape == (32, 32, 32)
        assert grid.min() >= 0


class TestSubhaloDataProducts:
    """Test SubhaloDataProducts HDF5 save/load and structure."""

    def test_init(self):
        """Test basic initialization."""
        dp = SubhaloDataProducts(ngrid=64)
        assert dp.ngrid == 64
        assert dp.dm_real is None
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
        finally:
            if os.path.exists(test_file):
                os.remove(test_file)
