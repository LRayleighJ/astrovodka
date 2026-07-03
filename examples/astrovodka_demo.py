#!/usr/bin/env python
"""
astrovodka demo script.

Demonstrates the basic usage of astrovodka for creating synthetic HI data cubes
from IllustrisTNG simulations with cosmological distances and K13 HI modeling.

Usage:
    python examples/astrovodka_demo.py

Or run in the TNG environment with actual data:
    python examples/astrovodka_demo.py --basePath /data/hw2/projects/FAST/HD2_high/TNG/TNG-100/ --snapNum 99 --subID 0
"""

import argparse
import numpy as np
import astropy.units as U
from astrovodka import Martini, DataCube
from astrovodka.beams import GaussianBeam
from astrovodka.noise import GaussianNoise
from astrovodka.spectral_models import GaussianSpectrum
from astrovodka.sph_kernels import CubicSplineKernel


def demo_simple_source(output_file="astrovodka_demo.fits"):
    """
    Create a mock observation using a simple analytical source at z=0.01.

    This demonstrates the cosmological distance mode of SPHSource.
    """
    print("=" * 60)
    print("  astrovodka Demo: Simple rotating disk at z=0.01")
    print("=" * 60)

    # Create a simple rotating disk galaxy
    N = 1000
    np.random.seed(42)
    phi = np.random.rand(N) * 2 * np.pi
    r = np.random.rand(N) * 10 * U.kpc  # 10 kpc exponential disk

    # Exponential vertical profile
    z_height = -np.log(np.random.rand(N)) * 0.5 * U.kpc
    z_height *= np.sign(np.random.rand(N) - 0.5)

    xyz_g = U.Quantity([
        np.zeros(N) * U.kpc,       # x: line-of-sight (thin disk)
        r * np.cos(phi),            # y: sky plane
        r * np.sin(phi) + z_height, # z: sky plane + vertical
    ])

    # Flat rotation curve
    vphi = 100 * U.km / U.s
    vxyz_g = U.Quantity([
        np.random.randn(N) * 5 * U.km / U.s,  # small random LOS velocities
        -vphi * np.sin(phi),
        vphi * np.cos(phi),
    ])

    T_g = np.ones(N) * 1e4 * U.K
    mHI_g = np.ones(N) * 1e8 * U.Msun  # 10^8 Msun per particle -> 10^11 Msun total
    hsm_g = np.ones(N) * 1.0 * U.kpc

    # ---- Redshift mode (new feature!) ----
    redshift = 0.01  # ~45 Mpc with Planck15

    from astrovodka.sources import SPHSource
    source = SPHSource(
        redshift=redshift,
        vpeculiar=0 * U.km / U.s,
        ra=30.0 * U.deg, dec=45.0 * U.deg,
        T_g=T_g, mHI_g=mHI_g, xyz_g=xyz_g, vxyz_g=vxyz_g, hsm_g=hsm_g,
    )

    print(f"\nSource properties:")
    print(f"  Redshift: {redshift}")
    print(f"  Luminosity distance: {source.distance:.1f}")
    print(f"  Systemic velocity: {source.vsys.to(U.km/U.s):.0f}")
    print(f"  Particles: {source.npart}")
    print(f"  Total HI mass: {source.input_mass:.2e}")

    # ---- DataCube ----
    datacube = DataCube(
        n_px_x=128, n_px_y=128, n_channels=64,
        px_size=2.0 * U.arcsec,
        channel_width=5.0 * U.km / U.s,
        spectral_centre=source.vsys,
        ra=source.ra, dec=source.dec,
    )
    print(f"\nDataCube: {datacube.n_px_x}x{datacube.n_px_y} pixels, "
          f"{datacube.n_channels} channels")

    # ---- Configure sub-modules ----
    beam = GaussianBeam(bmaj=8 * U.arcsec, bmin=6 * U.arcsec,
                        bpa=30 * U.deg, truncate=3.0)
    noise = GaussianNoise(rms=1e-4 * U.Jy * U.beam**-1, seed=0)
    spectral_model = GaussianSpectrum(sigma=7 * U.km / U.s)
    sph_kernel = CubicSplineKernel()

    # ---- Run pipeline ----
    print("\nRunning pipeline...")
    M = Martini(source=source, datacube=datacube, beam=beam, noise=noise,
                spectral_model=spectral_model, sph_kernel=sph_kernel)

    M.insert_source_in_cube()
    M.add_noise()
    M.convolve_beam()

    # ---- Write output ----
    M.write_fits(output_file)
    print(f"\nOutput written to: {output_file}")
    print("Done!")


def demo_tng_source(basePath, snapNum, subID):
    """Create a mock observation from real TNG data."""
    print("=" * 60)
    print(f"  astrovodka Demo: TNG-100 snap={snapNum} subID={subID}")
    print("=" * 60)

    from astrovodka.sources import TNGSource

    # Load TNG data with K13 HI model
    source = TNGSource(
        basePath=basePath,
        snapNum=snapNum,
        subID=subID,
        uv_mode="global",  # use global UVB (fast)
        ra=0.0 * U.deg,
        dec=0.0 * U.deg,
    )

    print(f"\nSource properties:")
    print(f"  Redshift: {source.source_redshift}")
    print(f"  Cosmology: Planck15")
    print(f"  Luminosity distance: {source.D_L:.1f}")
    print(f"  Angular diameter distance: {source.D_A:.1f}")
    print(f"  Systemic velocity: {source.vsys.to(U.km/U.s):.0f}")
    print(f"  Particles: {source.npart}")

    # DataCube centered on the source
    datacube = DataCube(
        n_px_x=256, n_px_y=256, n_channels=64,
        px_size=2.0 * U.arcsec,
        channel_width=10.0 * U.km / U.s,
        spectral_centre=source.vsys,
        ra=source.ra, dec=source.dec,
    )

    beam = GaussianBeam(bmaj=15 * U.arcsec, bmin=15 * U.arcsec,
                        bpa=0 * U.deg, truncate=3.0)
    noise = GaussianNoise(rms=1e-5 * U.Jy * U.beam**-1, seed=0)
    spectral_model = GaussianSpectrum(sigma="thermal")
    sph_kernel = CubicSplineKernel()

    M = Martini(source=source, datacube=datacube, beam=beam, noise=noise,
                spectral_model=spectral_model, sph_kernel=sph_kernel)

    M.insert_source_in_cube()
    M.add_noise()
    M.convolve_beam()

    outfile = f"astrovodka_TNG_snap{snapNum}_sub{str(subID).zfill(7)}.fits"
    M.write_fits(outfile)
    print(f"\nOutput written to: {outfile}")
    print("Done!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="astrovodka demo - synthetic HI observations"
    )
    parser.add_argument("--basePath", type=str, default=None,
                        help="Path to TNG data directory")
    parser.add_argument("--snapNum", type=int, default=99,
                        help="TNG snapshot number")
    parser.add_argument("--subID", type=int, default=0,
                        help="Subhalo ID")
    args = parser.parse_args()

    if args.basePath is not None:
        demo_tng_source(args.basePath, args.snapNum, args.subID)
    else:
        demo_simple_source()
