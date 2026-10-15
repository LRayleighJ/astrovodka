#!/usr/bin/env python
"""astrovodka demo script.

Creates mock resolved HI (21 cm) data cubes with astrovodka, either from a
simple analytic galaxy or from a real IllustrisTNG subhalo loaded off local
disk. With no arguments it runs the analytic galaxy, which needs no
simulation data; pass ``--basePath`` (or set ``TNG_BASE_PATH``) to run the TNG
path instead.

Usage
-----
    # Analytic rotating disc at z = 0.01, no simulation data required
    python examples/en/astrovodka_demo.py

    # A real TNG-100 subhalo (path can also come from $TNG_BASE_PATH)
    python examples/en/astrovodka_demo.py \\
        --basePath /path/to/TNG/TNG-100/ --snapNum 99 --subID 113 \\
        --hiModel K13_vol

Outputs are written to ``examples/outputs/`` by default (``--outdir`` to
change that); the directory is in .gitignore.
"""

import argparse
import os
from pathlib import Path

import numpy as np
import astropy.units as U

from astrovodka import Martini, DataCube
from astrovodka.beams import GaussianBeam
from astrovodka.noise import GaussianNoise
from astrovodka.spectral_models import GaussianSpectrum
from astrovodka.sph_kernels import CubicSplineKernel


def default_outdir() -> str:
    """Return ``examples/outputs/``, resolved relative to this file.

    Resolving against ``__file__`` rather than the working directory means the
    script can be run from anywhere in the repository and still write to the
    same place.
    """
    here = Path(__file__).resolve()
    outdir = here.parent.parent / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    return str(outdir)


def demo_simple_source(output_file):
    """Run the pipeline on an analytic rotating disc at z = 0.01.

    This is the whole mock-observation workflow in miniature and needs no
    simulation data: build a source, define a cube, pick sub-modules, run.

    Parameters
    ----------
    output_file : str
        Path of the FITS cube to write.
    """
    print("=" * 68)
    print("  astrovodka demo 1: analytic rotating disc at z = 0.01")
    print("=" * 68)

    # ---- An analytic disc galaxy: no simulation needed ----
    n_part = 1000
    rng = np.random.default_rng(42)
    phi = rng.random(n_part) * 2 * np.pi
    radius = rng.random(n_part) * 10 * U.kpc

    # Exponential vertical profile, symmetric about the mid-plane.
    z_height = -np.log(rng.random(n_part)) * 0.5 * U.kpc
    z_height *= np.sign(rng.random(n_part) - 0.5)

    xyz_g = U.Quantity(
        [
            np.zeros(n_part) * U.kpc,  # x is the line of sight
            radius * np.cos(phi),  # y-z is the sky plane
            radius * np.sin(phi) + z_height,
        ]
    )

    # Flat rotation curve of 100 km/s.
    v_phi = 100 * U.km / U.s
    vxyz_g = U.Quantity(
        [
            rng.standard_normal(n_part) * 5 * U.km / U.s,
            -v_phi * np.sin(phi),
            v_phi * np.cos(phi),
        ]
    )

    T_g = np.full(n_part, 1e4) * U.K
    mHI_g = np.full(n_part, 1e8) * U.Msun  # 10^8 Msun each, 10^11 Msun total
    hsm_g = np.full(n_part, 1.0) * U.kpc

    # ---- The source, placed at a cosmological redshift ----
    # Passing `redshift` (rather than `distance`) puts SPHSource in cosmology
    # mode: D_L comes from Planck15 and drives the flux normalisation, and D_A
    # drives every angular scale.
    from astrovodka.sources import SPHSource

    redshift = 0.01
    source = SPHSource(
        redshift=redshift,
        vpeculiar=0 * U.km / U.s,
        ra=30.0 * U.deg,
        dec=45.0 * U.deg,
        T_g=T_g,
        mHI_g=mHI_g,
        xyz_g=xyz_g,
        vxyz_g=vxyz_g,
        hsm_g=hsm_g,
    )

    print("\nSource")
    print(f"  redshift                 : {redshift}")
    print(f"  luminosity distance D_L  : {source.distance:.1f}")
    print(f"  systemic velocity        : {source.vsys.to(U.km / U.s):.0f}")
    print(f"  particles                : {source.npart}")
    print(f"  total HI mass            : {source.input_mass:.2e}")

    # ---- The data cube ----
    # Leaving spectral_unit at its default makes the spectral axis adopt the
    # units of channel_width (a velocity axis here), as in MARTINI.
    datacube = DataCube(
        n_px_x=128,
        n_px_y=128,
        n_channels=64,
        px_size=2.0 * U.arcsec,
        channel_width=5.0 * U.km / U.s,
        spectral_centre=source.vsys,
        ra=source.ra,
        dec=source.dec,
    )
    print(
        f"\nDataCube: {datacube.n_px_x} x {datacube.n_px_y} pixels, "
        f"{datacube.n_channels} channels"
    )

    # ---- Instrument and line models ----
    beam = GaussianBeam(bmaj=8 * U.arcsec, bmin=6 * U.arcsec, bpa=30 * U.deg, truncate=3.0)
    noise = GaussianNoise(rms=1e-4 * U.Jy * U.beam**-1, seed=0)
    spectral_model = GaussianSpectrum(sigma=7 * U.km / U.s)
    sph_kernel = CubicSplineKernel()

    # ---- Run the pipeline ----
    print("\nRunning pipeline...")
    martini = Martini(
        source=source,
        datacube=datacube,
        beam=beam,
        noise=noise,
        spectral_model=spectral_model,
        sph_kernel=sph_kernel,
    )
    martini.insert_source_in_cube()
    martini.add_noise()
    martini.convolve_beam()
    martini.write_fits(output_file)

    print(f"\nWrote {output_file}")
    print("Done.")


def demo_tng_source(basePath, snapNum, subID, observed_redshift, hi_model, output_file):
    """Run the pipeline on a real IllustrisTNG subhalo read from local disk.

    Parameters
    ----------
    basePath : str
        TNG simulation output directory, containing ``groups_<snap>/`` and
        ``snapdir_<snap>/``.
    snapNum : int
        Snapshot number.
    subID : int
        Subhalo ID within that snapshot.
    observed_redshift : float
        Redshift to place the mock observation at.
    hi_model : str
        One of ``K13_vol``, ``K13_map``, ``GD14_vol``, ``GD14_map``.
    output_file : str
        Path of the FITS cube to write.
    """
    print("=" * 68)
    print(f"  astrovodka demo 2: TNG-100 snap={snapNum} subID={subID}")
    print("=" * 68)

    from astrovodka.sources import TNGSource

    # `observed_redshift` separates the two redshifts that matter here. The
    # snapshot's own redshift sets the physics (unit conversions, the UV
    # background); `observed_redshift` sets the geometry (D_L, D_A, and so the
    # angular size and the flux). For a z~0 snapshot such as snap 99 the
    # angular diameter distance is essentially zero, so the galaxy would be an
    # unresolvable point without this.
    #
    # `hi_model` picks one of four HI/H2 models. The projected (*_map) ones
    # also need the stellar particles, and assign zero HI to gas that falls
    # outside the projection box.
    source = TNGSource(
        basePath=basePath,
        snapNum=snapNum,
        subID=subID,
        hi_model=hi_model,
        observed_redshift=observed_redshift,
        ra=0.0 * U.deg,
        dec=0.0 * U.deg,
    )

    print("\nSource")
    print(f"  HI model                 : {source.hi_model}")
    print(f"  snapshot redshift        : {source.source_redshift}")
    print(f"  observed redshift        : {observed_redshift}")
    print(f"  cosmology                : Planck15")
    print(f"  luminosity distance D_L  : {source.D_L:.1f}")
    print(f"  angular diameter dist D_A: {source.D_A:.1f}")
    print(f"  systemic velocity        : {source.vsys.to(U.km / U.s):.0f}")
    print(f"  particles                : {source.npart}")
    print(f"  total HI mass            : {source.m_HI_total:.4e}")
    print(f"  total H2 mass            : {source.m_H2_g.sum():.4e}")

    datacube = DataCube(
        n_px_x=256,
        n_px_y=256,
        n_channels=64,
        px_size=2.0 * U.arcsec,
        channel_width=10.0 * U.km / U.s,
        spectral_centre=source.vsys,
        ra=source.ra,
        dec=source.dec,
    )

    beam = GaussianBeam(bmaj=15 * U.arcsec, bmin=15 * U.arcsec, bpa=0 * U.deg, truncate=3.0)
    noise = GaussianNoise(rms=1e-5 * U.Jy * U.beam**-1, seed=0)
    spectral_model = GaussianSpectrum(sigma="thermal")  # kB*T/m_p per particle
    sph_kernel = CubicSplineKernel()

    martini = Martini(
        source=source,
        datacube=datacube,
        beam=beam,
        noise=noise,
        spectral_model=spectral_model,
        sph_kernel=sph_kernel,
    )
    martini.insert_source_in_cube()
    martini.add_noise()
    martini.convolve_beam()
    martini.write_fits(output_file)

    print(f"\nWrote {output_file}")
    print("Done.")


def main():
    parser = argparse.ArgumentParser(
        description="astrovodka demo -- synthetic HI observations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--basePath",
        type=str,
        default=os.environ.get("TNG_BASE_PATH"),
        help="TNG simulation output directory. Falls back to $TNG_BASE_PATH. "
        "If neither is set, the analytic-galaxy demo runs instead.",
    )
    parser.add_argument("--snapNum", type=int, default=99, help="TNG snapshot number.")
    parser.add_argument(
        "--subID",
        type=int,
        default=113,
        help="Subhalo ID. The default, 113, is a gas-rich, star-forming subhalo "
        "of TNG-100 snap 99 (M_HI ~ 4e9 Msun). Avoid low IDs: the catalog is "
        "sorted by descending mass, so ID 0 is the cluster core (~38 million "
        "gas cells, >10 GB to load).",
    )
    parser.add_argument(
        "--observedRedshift",
        type=float,
        default=0.02,
        help="Redshift to place the mock observation at, decoupled from the "
        "snapshot redshift.",
    )
    parser.add_argument(
        "--hiModel",
        type=str,
        default="K13_vol",
        choices=["K13_vol", "K13_map", "GD14_vol", "GD14_map"],
        help="Which HI/H2 model to use.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help="Directory for the output FITS file. Defaults to examples/outputs/.",
    )
    args = parser.parse_args()

    outdir = args.outdir or default_outdir()
    os.makedirs(outdir, exist_ok=True)

    if args.basePath:
        basePath = args.basePath.rstrip("/") + "/"
        output_file = os.path.join(
            outdir,
            f"astrovodka_TNG_snap{args.snapNum}_sub{args.subID:07d}"
            f"_z{args.observedRedshift:.2f}.fits",
        )
        demo_tng_source(
            basePath,
            args.snapNum,
            args.subID,
            observed_redshift=args.observedRedshift,
            hi_model=args.hiModel,
            output_file=output_file,
        )
    else:
        output_file = os.path.join(outdir, "astrovodka_demo.fits")
        demo_simple_source(output_file)


if __name__ == "__main__":
    main()
