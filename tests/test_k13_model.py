"""Test K13 HI model functions."""

import numpy as np
from astrovodka.sources._k13_model import (
    calc_UMW_global_Faucher11,
    get_f_neutral_H,
    f_mol_K13_vol,
    compute_HI_mass_particle,
)
from astrovodka.cosmology_utils import rho_code_to_nH, rho_code_to_physical


class TestUVB:
    """Test global UV background computation."""

    def test_output_is_positive(self):
        """U_MW should be positive at all redshifts."""
        for z in [0.0, 0.5, 1.0, 2.0, 5.0]:
            u = calc_UMW_global_Faucher11(z)
            assert u > 0.0

    def test_increases_with_z(self):
        """UVB should generally increase with redshift (up to ~z=2-3)."""
        u0 = calc_UMW_global_Faucher11(0.0)
        u1 = calc_UMW_global_Faucher11(1.0)
        u2 = calc_UMW_global_Faucher11(2.0)
        assert u1 > u0
        assert u2 > u0

    def test_at_z0(self):
        """At z=0, U_MW should be ~1.4e-3 (Milky Way = 1)."""
        u = calc_UMW_global_Faucher11(0.0)
        assert 1e-4 < u < 1e-2


class TestFNeutral:
    """Test neutral fraction estimation."""

    def test_sf_particles_fully_neutral(self):
        """Star-forming particles should have f_neutral = 1."""
        nH = np.array([0.5, 0.8, 0.3])
        sfr = np.array([1.0, 0.5, 2.0])  # all SF
        f = get_f_neutral_H(nH, sfr)
        assert np.all(f == 1.0)

    def test_non_sf_preserves_fHI(self):
        """Non-SF particles should preserve f_HI value."""
        nH = np.array([0.5, 0.8, 0.3])
        sfr = np.array([0.0, 0.0, 0.0])  # all non-SF
        f = get_f_neutral_H(nH, sfr)
        assert np.allclose(f, nH)

    def test_mixed_population(self):
        """Mixed SF and non-SF population."""
        nH = np.array([0.5, 0.8, 0.3, 0.6])
        sfr = np.array([0.0, 1.0, 0.0, 0.5])
        f = get_f_neutral_H(nH, sfr)
        assert f[0] == 0.5  # non-SF
        assert f[1] == 1.0  # SF
        assert f[2] == 0.3  # non-SF
        assert f[3] == 1.0  # SF


class TestDensityConversion:
    """Test TNG density unit conversions."""

    def test_rho_code_to_physical_positive(self):
        """Output should be positive."""
        rho = np.array([1.0, 10.0, 100.0])
        phys = rho_code_to_physical(rho, 0.0, 0.6774)
        assert np.all(phys > 0)

    def test_rho_code_to_physical_redshift_scaling(self):
        """Physical density scales as (1+z)^3."""
        rho = np.array([1.0])
        phys_z0 = rho_code_to_physical(rho, 0.0, 0.6774)
        phys_z1 = rho_code_to_physical(rho, 1.0, 0.6774)
        ratio = phys_z1 / phys_z0
        assert abs(ratio[0] - 8.0) < 0.01  # (1+1)^3 = 8

    def test_rho_code_to_nH_positive(self):
        """nH should be positive."""
        rho = np.array([1.0, 10.0])
        nh = rho_code_to_nH(rho, 0.0, 0.6774)
        assert np.all(nh > 0)


class TestK13Model:
    """Test the core K13 f_mol model."""

    def make_test_data(self, n=10):
        """Create synthetic test data for K13 model."""
        np.random.seed(42)
        rho_code = np.ones(n) * 1e6  # typical ISM density
        u_code = np.ones(n) * 1.0    # internal energy
        Z = np.ones(n) * 0.0127      # solar metallicity
        nH_code = np.ones(n) * 0.9   # mostly neutral
        sfr_code = np.zeros(n)       # no SF
        X_H = np.ones(n) * 0.76
        U_MW = 0.001  # low UV field
        return rho_code, u_code, Z, nH_code, sfr_code, X_H, U_MW

    def test_output_range(self):
        """f_mol, f_HI, f_H2 should all be in [0, 1]."""
        data = self.make_test_data(10)
        f_mol, f_HI, f_H2 = f_mol_K13_vol(*data)
        for arr, name in [(f_mol, 'f_mol'), (f_HI, 'f_HI'), (f_H2, 'f_H2')]:
            assert np.all(arr >= 0), f"{name} has negative values"
            assert np.all(arr <= 1), f"{name} has values > 1"

    def test_sum_constraint(self):
        """f_HI + f_H2 = f_neutral <= 1."""
        data = self.make_test_data(20)
        f_mol, f_HI, f_H2 = f_mol_K13_vol(*data)
        total = f_HI + f_H2
        assert np.all(total <= 1.0 + 1e-10)

    def test_high_metallicity_more_molecular(self):
        """Higher metallicity should promote H2 formation."""
        n = 10
        rho_code = np.ones(n) * 1e6
        u_code = np.ones(n) * 1.0
        Z_low = np.ones(n) * 0.001  # 0.1 Z_sun
        Z_high = np.ones(n) * 0.025  # 2 Z_sun
        nH_code = np.ones(n) * 0.9
        sfr_code = np.zeros(n)
        X_H = np.ones(n) * 0.76
        U_MW = 0.001

        _, _, f_H2_low = f_mol_K13_vol(
            rho_code, u_code, Z_low, nH_code, sfr_code, X_H, U_MW)
        _, _, f_H2_high = f_mol_K13_vol(
            rho_code, u_code, Z_high, nH_code, sfr_code, X_H, U_MW)

        assert np.mean(f_H2_high) >= np.mean(f_H2_low), \
            "Higher Z should give more H2"

    def test_high_uv_less_molecular(self):
        """Higher UV field should reduce H2 fraction."""
        n = 10
        rho_code = np.ones(n) * 1e6
        u_code = np.ones(n) * 1.0
        Z = np.ones(n) * 0.0127
        nH_code = np.ones(n) * 0.9
        sfr_code = np.zeros(n)
        X_H = np.ones(n) * 0.76

        _, _, f_H2_low_uv = f_mol_K13_vol(
            rho_code, u_code, Z, nH_code, sfr_code, X_H, 0.001)
        _, _, f_H2_high_uv = f_mol_K13_vol(
            rho_code, u_code, Z, nH_code, sfr_code, X_H, 100.0)

        assert np.mean(f_H2_high_uv) <= np.mean(f_H2_low_uv), \
            "Higher UV should give less H2"


class TestComputeHIMass:
    """Test compute_HI_mass_particle."""

    def test_output_shapes(self):
        """Output arrays should match input particle count."""
        n = 50
        gas_data = {
            'Masses': np.ones(n) * 0.01,
            'Density': np.ones(n) * 1e6,
            'InternalEnergy': np.ones(n) * 1.0,
            'NeutralHydrogenAbundance': np.ones(n) * 0.9,
            'StarFormationRate': np.zeros(n),
            'GFM_Metals': np.ones(n) * 0.76,  # 1D: already sliced to H column
            'GFM_Metallicity': np.ones(n) * 0.0127,
        }
        U_MW = 0.001
        m_HI, m_H2 = compute_HI_mass_particle(gas_data, U_MW)
        assert len(m_HI) == n
        assert len(m_H2) == n

    def test_mass_conservation(self):
        """Total H mass should not exceed gas mass * X_H."""
        n = 50
        mgas = 0.01
        X_H_val = 0.76
        gas_data = {
            'Masses': np.ones(n) * mgas,
            'Density': np.ones(n) * 1e6,
            'InternalEnergy': np.ones(n) * 1.0,
            'NeutralHydrogenAbundance': np.ones(n) * 0.9,
            'StarFormationRate': np.zeros(n),
            'GFM_Metals': np.ones((n, 1)) * X_H_val,
            'GFM_Metallicity': np.ones(n) * 0.0127,
        }
        U_MW = 0.001
        m_HI, m_H2 = compute_HI_mass_particle(gas_data, U_MW)
        total_H = (m_HI + m_H2)  # in code units
        max_possible = mgas * X_H_val
        assert np.all(total_H <= max_possible * (1 + 1e-10))
