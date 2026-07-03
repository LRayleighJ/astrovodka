"""Integration tests for SPHSource with cosmology mode."""

import numpy as np
import astropy.units as U
import pytest
from astrovodka import Martini, DataCube
from astrovodka.sources import SPHSource
from astrovodka.beams import GaussianBeam
from astrovodka.spectral_models import GaussianSpectrum
from astrovodka.sph_kernels import CubicSplineKernel


def make_test_source(N=100, redshift=0.01, r_kpc=10.0):
    """Create a simple rotating disk source for testing."""
    np.random.seed(42)
    phi = np.random.rand(N) * 2 * np.pi
    r = np.random.rand(N) * r_kpc * U.kpc

    xyz_g = U.Quantity([
        np.zeros(N) * U.kpc,
        r * np.cos(phi),
        r * np.sin(phi),
    ])
    vphi = 100 * U.km / U.s
    vxyz_g = U.Quantity([
        np.zeros(N) * U.km / U.s,
        -vphi * np.sin(phi),
        vphi * np.cos(phi),
    ])
    T_g = np.ones(N) * 1e4 * U.K
    mHI_g = np.ones(N) * 1e8 * U.Msun
    hsm_g = np.ones(N) * 1.0 * U.kpc

    return SPHSource(
        redshift=redshift,
        vpeculiar=0 * U.km / U.s,
        ra=0.0 * U.deg, dec=0.0 * U.deg,
        T_g=T_g, mHI_g=mHI_g, xyz_g=xyz_g, vxyz_g=vxyz_g, hsm_g=hsm_g,
    )


class TestSPHSourceCosmology:
    """Test SPHSource with cosmological distances."""

    def test_redshift_mode_sets_distance(self):
        """Redshift mode should compute luminosity distance."""
        source = make_test_source(redshift=0.01)
        assert source.distance.value > 0
        assert source.distance.unit == U.Mpc

    def test_vsys_from_redshift(self):
        """vsys should be ~c*z for low z."""
        source = make_test_source(redshift=0.01)
        expected_vsys = 2998  # km/s for z=0.01
        assert abs(source.vsys.to_value(U.km / U.s) - expected_vsys) < 20

    def test_redshift_vs_distance_mode(self):
        """Redshift mode should give similar results to distance mode at z=0.01."""
        source_redshift = make_test_source(redshift=0.01)
        # For z=0.01 with Planck15, D_L ~ 44.6 Mpc, h ~ 0.6774
        # vhubble = h * 100 * D_L ~ 0.6774 * 100 * 44.6 ~ 3021 km/s
        source_distance = SPHSource(
            distance=44.6 * U.Mpc,
            h=0.6774,
            vpeculiar=0 * U.km / U.s,
            ra=0.0 * U.deg, dec=0.0 * U.deg,
            T_g=np.ones(100) * 1e4 * U.K,
            mHI_g=np.ones(100) * 1e8 * U.Msun,
            xyz_g=U.Quantity([
                np.zeros(100) * U.kpc,
                np.random.rand(100) * 10 * U.kpc,
                np.random.rand(100) * 10 * U.kpc,
            ]),
            vxyz_g=U.Quantity([
                np.zeros(100) * U.km / U.s,
                np.random.randn(100) * 50 * U.km / U.s,
                np.random.randn(100) * 50 * U.km / U.s,
            ]),
            hsm_g=np.ones(100) * 1.0 * U.kpc,
        )
        # Distances should be similar (within 5%)
        assert abs(source_redshift.distance.value - source_distance.distance.value) < 5

    def test_error_when_no_distance_or_redshift(self):
        """Should raise ValueError if neither distance nor redshift given."""
        with pytest.raises(ValueError, match="Either 'distance'"):
            SPHSource(
                vpeculiar=0 * U.km / U.s,
                ra=0.0 * U.deg, dec=0.0 * U.deg,
                T_g=np.ones(10) * 1e4 * U.K,
                mHI_g=np.ones(10) * 1e8 * U.Msun,
                xyz_g=U.Quantity([
                    np.zeros(10) * U.kpc,
                    np.random.rand(10) * 10 * U.kpc,
                    np.random.rand(10) * 10 * U.kpc,
                ]),
                vxyz_g=U.Quantity([
                    np.zeros(10) * U.km / U.s,
                    np.random.randn(10) * 50 * U.km / U.s,
                    np.random.randn(10) * 50 * U.km / U.s,
                ]),
                hsm_g=np.ones(10) * 1.0 * U.kpc,
            )

    def test_error_when_both_distance_and_redshift(self):
        """Should raise ValueError if both distance and redshift given."""
        with pytest.raises(ValueError, match="not both"):
            SPHSource(
                distance=10 * U.Mpc,
                redshift=0.01,
                vpeculiar=0 * U.km / U.s,
                ra=0.0 * U.deg, dec=0.0 * U.deg,
                T_g=np.ones(10) * 1e4 * U.K,
                mHI_g=np.ones(10) * 1e8 * U.Msun,
                xyz_g=U.Quantity([
                    np.zeros(10) * U.kpc,
                    np.random.rand(10) * 10 * U.kpc,
                    np.random.rand(10) * 10 * U.kpc,
                ]),
                vxyz_g=U.Quantity([
                    np.zeros(10) * U.km / U.s,
                    np.random.randn(10) * 50 * U.km / U.s,
                    np.random.randn(10) * 50 * U.km / U.s,
                ]),
                hsm_g=np.ones(10) * 1.0 * U.kpc,
            )


class TestPipelineWithCosmology:
    """End-to-end pipeline tests with cosmology mode."""

    def test_full_pipeline_low_z(self):
        """Full pipeline should run without errors at z=0.01."""
        source = make_test_source(N=200, redshift=0.01, r_kpc=10.0)

        datacube = DataCube(
            n_px_x=64, n_px_y=64, n_channels=32,
            px_size=2.0 * U.arcsec,
            channel_width=10.0 * U.km / U.s,
            spectral_centre=source.vsys,
        )

        beam = GaussianBeam(bmaj=10 * U.arcsec, bmin=10 * U.arcsec,
                            bpa=0 * U.deg, truncate=3.0)
        spectral_model = GaussianSpectrum(sigma=7 * U.km / U.s)
        sph_kernel = CubicSplineKernel()

        M = Martini(source=source, datacube=datacube, beam=beam, noise=None,
                    spectral_model=spectral_model, sph_kernel=sph_kernel)
        M.insert_source_in_cube(progressbar=False)
        M.convolve_beam()

        # Check that the cube has non-zero flux
        assert datacube._array.max() > 0

    def test_source_at_moderate_z(self):
        """Source at z=0.1 should have correct distance and vsys."""
        source = make_test_source(N=100, redshift=0.1, r_kpc=10.0)
        # D_L ~ 460 Mpc at z=0.1 with Planck15
        assert 400 < source.distance.value < 520
        # vsys ~ 29500 km/s at z=0.1
        assert 28000 < source.vsys.to_value(U.km / U.s) < 31000
        # angular_diameter_distance should be set
        assert hasattr(source, 'angular_diameter_distance')
