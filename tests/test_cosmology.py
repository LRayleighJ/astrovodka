"""Test cosmological distance computations."""

import numpy as np
import astropy.units as U
from astropy.cosmology import Planck15, FlatLambdaCDM
from astrovodka.cosmology_utils import (
    get_cosmology,
    luminosity_distance,
    angular_diameter_distance,
    comoving_distance,
    hubble_parameter,
    systemic_velocity_from_redshift,
    observed_frequency,
    HIfreq,
)


class TestGetCosmology:
    """Test cosmology factory function."""

    def test_default_planck15(self):
        """Test that default returns Planck15."""
        c = get_cosmology()
        assert c is Planck15

    def test_custom_cosmology(self):
        """Test with a custom cosmology."""
        custom = FlatLambdaCDM(H0=70, Om0=0.3)
        c = get_cosmology(custom)
        assert c is custom
        assert c.H0.value == 70

    def test_invalid_type(self):
        """Test that invalid type raises TypeError."""
        try:
            get_cosmology("not_a_cosmology")
            assert False, "Should have raised TypeError"
        except TypeError:
            pass


class TestDistances:
    """Test cosmological distance computations."""

    def test_luminosity_distance_z0(self):
        """D_L at z=0 should use Hubble-law floor (~0.004 Mpc) to avoid div-by-0."""
        dl = luminosity_distance(0.0)
        # Floor gives cz/H0 for z_eff=1e-6 ~ 0.0044 Mpc
        assert 0.001 < dl.value < 0.01

    def test_luminosity_distance_low_z(self):
        """D_L ~ cz/H0 at low z."""
        z = 0.01
        dl = luminosity_distance(z)
        expected = 3e5 * z / 67.74  # ~44 Mpc for Planck15
        assert 40 < dl.value < 50

    def test_angular_diameter_distance(self):
        """D_A = D_L / (1+z)^2 (Etherington relation)."""
        for z in [0.1, 0.5, 1.0, 2.0]:
            dl = luminosity_distance(z)
            da = angular_diameter_distance(z)
            ratio = (dl / da).to_value(U.dimensionless_unscaled)
            expected = (1 + z) ** 2
            assert abs(ratio - expected) < 1e-6

    def test_comoving_distance(self):
        """Comoving distance = D_A * (1+z)."""
        for z in [0.1, 0.5, 1.0]:
            dc = comoving_distance(z)
            da = angular_diameter_distance(z)
            ratio = (dc / da).to_value(U.dimensionless_unscaled)
            assert abs(ratio - (1 + z)) < 1e-6

    def test_distance_increases_with_z(self):
        """D_L should increase monotonically with z."""
        dl_01 = luminosity_distance(0.1)
        dl_05 = luminosity_distance(0.5)
        dl_10 = luminosity_distance(1.0)
        assert dl_01 < dl_05 < dl_10


class TestHubbleParameter:
    """Test H(z) computation."""

    def test_hubble_z0(self):
        """H(0) should be close to H0."""
        hz = hubble_parameter(0.0)
        assert 60 < hz.value < 80  # Planck15 H0 ~ 67.74

    def test_hubble_increases_with_z(self):
        """H(z) increases with z in Lambda-CDM."""
        h0 = hubble_parameter(0.0)
        h1 = hubble_parameter(1.0)
        assert h1 > h0


class TestSystemicVelocity:
    """Test systemic velocity from redshift."""

    def test_low_z_approx(self):
        """At low z, vsys ~ c*z."""
        z = 0.01
        vsys = systemic_velocity_from_redshift(z)
        expected = 3e5 * z  # km/s
        assert abs(vsys.value - expected) / expected < 0.01

    def test_vsys_less_than_c(self):
        """Systemic velocity should never exceed c."""
        for z in [0.1, 1.0, 5.0, 10.0]:
            vsys = systemic_velocity_from_redshift(z)
            assert vsys.to(U.km / U.s).value < 3e5


class TestObservedFrequency:
    """Test observed frequency computation."""

    def test_z0(self):
        """At z=0, observed frequency = rest frequency."""
        fobs = observed_frequency(0.0)
        assert abs(fobs.value - HIfreq.value) < 1e-6

    def test_redshift_relation(self):
        """f_obs = f_rest / (1+z)."""
        z = 1.0
        fobs = observed_frequency(z)
        assert abs(fobs.value - HIfreq.value / 2) < 1e-6
