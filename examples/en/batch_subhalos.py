#!/usr/bin/env python
"""Batch-process subhalos from an IllustrisTNG group catalog.

Reads the subhalo catalog, selects objects by gas-cell count, and builds one
mock HI data cube per object. This is the pattern to reach for when you want a
sample rather than a single galaxy.

Usage
-----
    export TNG_BASE_PATH=/path/to/TNG/TNG-100/
    python examples/en/batch_subhalos.py --limit 20
    python examples/en/batch_subhalos.py --limit 10 --hiModel K13_map --snapNum 99

Outputs land in ``examples/outputs/batch/`` unless ``--outdir`` says otherwise.

A note on subhalo selection
---------------------------
The group catalog is sorted by **descending mass**, so the first entries are
cluster cores: subhalo 0 of TNG-100 snapshot 99 has roughly 38 million gas
cells and needs more than 10 GB to load. Taking ``ids[:N]`` is therefore a good
way to hang a login node.

Filtering on the gas-cell count first — which is what ``--minGas``/``--maxGas``
do — selects the low-mass, gas-rich satellites that are actually interesting
for HI work, and bounds the memory per object.
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import astropy.units as U
import h5py

from astrovodka import Martini, DataCube
from astrovodka.beams import GaussianBeam
from astrovodka.spectral_models import GaussianSpectrum
from astrovodka.sph_kernels import CubicSplineKernel
from astrovodka.sources import TNGSource


def default_outdir() -> str:
    """Return ``examples/outputs/batch/``, resolved relative to this file."""
    outdir = Path(__file__).resolve().parent.parent / "outputs" / "batch"
    outdir.mkdir(parents=True, exist_ok=True)
    return str(outdir)


def _read_catalog(basePath, snapNum):
    """Read the gas-cell count and SFR of every subhalo in the catalog.

    The catalog is split across several hundred chunk files, and the chunks are
    sequential in subhalo ID, so concatenating them in numeric order reproduces
    the global ID ordering. Reading only the first chunk is a trap worth
    avoiding: because the catalog is sorted by descending mass, chunk 0 holds
    the most massive objects in the box and nothing else.

    Returns
    -------
    n_gas, sfr : numpy.ndarray
        Per-subhalo gas-cell count and star formation rate [Msun/yr], indexed
        by subhalo ID.
    """
    catalog_dir = Path(basePath) / f"groups_{snapNum:03d}"
    chunks = sorted(
        catalog_dir.glob(f"fof_subhalo_tab_{snapNum:03d}.*.hdf5"),
        key=lambda p: int(p.name.rsplit(".", 2)[-2]),
    )
    if not chunks:
        raise FileNotFoundError(
            f"No group catalog files in {catalog_dir}.\n"
            "basePath must be the directory that *contains* groups_<snap>/ and "
            "snapdir_<snap>/, and that snapshot must have been downloaded."
        )

    n_gas_parts, sfr_parts = [], []
    for chunk in chunks:
        with h5py.File(chunk, "r") as handle:
            subhalo = handle.get("Subhalo")
            # Some chunks are empty placeholders: zero subhalos in that ID
            # range, and the datasets simply do not exist in the file. They
            # contribute nothing, so skip rather than special-case them.
            if subhalo is None or "SubhaloLenType" not in subhalo:
                continue
            # SubhaloLenType[:, 0] is the gas-cell count, part type 0 being gas.
            n_gas_parts.append(np.asarray(subhalo["SubhaloLenType"][:, 0]))
            sfr_parts.append(np.asarray(subhalo["SubhaloSFR"]))
    if not n_gas_parts:
        raise RuntimeError(f"no subhalos found in any chunk under {catalog_dir}")
    return np.concatenate(n_gas_parts), np.concatenate(sfr_parts)


def select_subhalos(basePath, snapNum, min_gas, max_gas, limit, min_sfr=0.0):
    """Return subhalo IDs with a gas-cell count in ``[min_gas, max_gas]``.

    Parameters
    ----------
    basePath : str
        TNG simulation directory, containing ``groups_<snap>/``.
    snapNum : int
        Snapshot number.
    min_gas, max_gas : int
        Inclusive bounds on the number of gas cells.
    limit : int or None
        Keep at most this many, taking the gas-richest first.
    min_sfr : float, optional
        Minimum star formation rate [Msun/yr]. Star-forming objects are the
        interesting ones for HI work: the adaptive FFT UV field only does
        visible work where there are young stars.

    Returns
    -------
    list of tuple
        ``(subhalo_id, n_gas, sfr)`` for each selected object.
    """
    n_gas, sfr = _read_catalog(basePath, snapNum)
    print(f"Catalog: {len(n_gas)} subhalos in {snapNum:03d}")

    inside = np.flatnonzero(
        (n_gas >= min_gas) & (n_gas <= max_gas) & (sfr >= min_sfr)
    )
    # Gas-richest first, so --limit keeps the most interesting objects.
    inside = inside[np.argsort(n_gas[inside])[::-1]]
    if limit is not None:
        inside = inside[:limit]

    print(
        f"Selected {len(inside)} with {min_gas} <= N_gas <= {max_gas} "
        f"and SFR >= {min_sfr} Msun/yr"
    )
    return [(int(i), int(n_gas[i]), float(sfr[i])) for i in inside]


def build_cube(basePath, snapNum, subID, hi_model, observed_redshift, outdir):
    """Build and write one mock observation.

    Returns
    -------
    dict or None
        A summary record, or None if the subhalo produced no HI to observe.
    """
    source = TNGSource(
        basePath=basePath,
        snapNum=snapNum,
        subID=subID,
        hi_model=hi_model,
        observed_redshift=observed_redshift,
        ra=0.0 * U.deg,
        dec=0.0 * U.deg,
    )

    if source.m_HI_total <= 0:
        return None

    # Size the cube to the source: a fixed angular size would either clip a
    # nearby galaxy or waste pixels on a distant one. The velocity axis needs
    # the gas kinematics plus room for the line profile, since particles whose
    # line falls outside the channel range are pruned silently (see the guide).
    arcsec_per_kpc = (U.rad / source.D_A.to(U.kpc)).to(U.arcsec / U.kpc)
    extent_arcsec = 2.5 * 20.0 * arcsec_per_kpc.value  # ~2.5x a 20 kpc galaxy
    beam_fwhm = 15.0 * U.arcsec
    n_px = int(np.clip(np.ceil(extent_arcsec / (beam_fwhm.value / 3.0)), 32, 256))

    datacube = DataCube(
        n_px_x=n_px,
        n_px_y=n_px,
        n_channels=64,
        px_size=(extent_arcsec / n_px) * U.arcsec,
        channel_width=20.0 * U.km / U.s,
        spectral_centre=source.vsys,
        ra=source.ra,
        dec=source.dec,
    )

    martini = Martini(
        source=source,
        datacube=datacube,
        beam=GaussianBeam(bmaj=beam_fwhm),
        # A fresh spectral model per cube: one instance caches its spectra
        # array and cannot be reused across cubes of different channel count.
        spectral_model=GaussianSpectrum(sigma="thermal"),
        sph_kernel=CubicSplineKernel(),
        quiet=True,
    )
    martini.insert_source_in_cube(progressbar=False)
    martini.convolve_beam()

    outfile = os.path.join(
        outdir,
        f"TNG{snapNum}_sub{subID:07d}_{hi_model}_z{observed_redshift:.2f}.fits",
    )
    martini.write_fits(outfile)

    return {
        "subID": subID,
        "m_hi": float(source.m_HI_total.to_value(U.Msun)),
        "m_h2": float(source.m_H2_g.sum().to_value(U.Msun)),
        "n_px": n_px,
        "npart": source.npart,
        "outfile": outfile,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Build mock HI cubes for a set of TNG subhalos.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--basePath",
        default=os.environ.get("TNG_BASE_PATH"),
        help="TNG simulation directory. Falls back to $TNG_BASE_PATH.",
    )
    parser.add_argument("--snapNum", type=int, default=99, help="Snapshot number.")
    parser.add_argument("--limit", type=int, default=10, help="How many subhalos to process.")
    parser.add_argument("--minGas", type=int, default=64, help="Minimum gas cells per subhalo.")
    parser.add_argument(
        "--maxGas",
        type=int,
        default=5000,
        help="Maximum gas cells per subhalo. Keeps cluster cores, which are "
        "slow and memory-hungry, out of the sample.",
    )
    parser.add_argument(
        "--hiModel",
        default="K13_vol",
        choices=["K13_vol", "K13_map", "GD14_vol", "GD14_map"],
        help="HI/H2 model.",
    )
    parser.add_argument("--observedRedshift", type=float, default=0.02, help="Observing redshift.")
    parser.add_argument(
        "--minSfr",
        type=float,
        default=0.0,
        help="Minimum star formation rate [Msun/yr]. Raising it above zero "
        "restricts the sample to star-forming objects, where the adaptive FFT "
        "UV field does visible work.",
    )
    parser.add_argument("--outdir", default=None, help="Output directory.")
    args = parser.parse_args()

    if not args.basePath:
        parser.error(
            "no TNG data directory. Pass --basePath, or set TNG_BASE_PATH:\n"
            "    export TNG_BASE_PATH=/path/to/TNG/TNG-100/"
        )
    basePath = args.basePath.rstrip("/") + "/"
    outdir = args.outdir or default_outdir()
    os.makedirs(outdir, exist_ok=True)

    subhalos = select_subhalos(
        basePath, args.snapNum, args.minGas, args.maxGas, args.limit, args.minSfr
    )
    if not subhalos:
        print("Nothing selected; loosen --minGas/--maxGas.")
        return 1

    print(f"\nModel {args.hiModel}, observing at z = {args.observedRedshift}")
    print(f"Writing to {outdir}\n")

    results, skipped, failures = [], [], []
    started = time.time()
    for index, (subID, n_gas, sfr) in enumerate(subhalos, start=1):
        print(
            f"[{index}/{len(subhalos)}] subID {subID} "
            f"({n_gas} gas cells, SFR {sfr:.2f} Msun/yr)...",
            flush=True,
        )
        try:
            record = build_cube(
                basePath, args.snapNum, subID, args.hiModel, args.observedRedshift, outdir
            )
        except Exception as exc:  # keep going: one bad object should not stop the batch
            failures.append((subID, repr(exc)))
            print(f"    failed: {exc}")
            continue

        if record is None:
            skipped.append(subID)
            print("    no HI in this model; skipped")
            continue

        results.append(record)
        print(f"    M_HI = {record['m_hi']:.3e} Msun, {record['npart']} particles used")

    elapsed = time.time() - started
    print("\n" + "=" * 68)
    print(f"  {len(results)} cubes written, {len(skipped)} skipped, "
          f"{len(failures)} failed  ({elapsed:.0f} s)")
    print("=" * 68)

    if results:
        m_hi = np.array([r["m_hi"] for r in results])
        m_h2 = np.array([r["m_h2"] for r in results])
        print(f"  M_HI: median {np.median(m_hi):.3e}, max {m_hi.max():.3e} Msun")
        print(f"  M_H2: median {np.median(m_h2):.3e}, max {m_h2.max():.3e} Msun")
        print(f"\n  These are {args.hiModel} masses. The four models disagree with")
        print("  each other by design — on the same galaxy the H2 mass can differ by")
        print("  a factor of several — so any number you publish needs the model named.")
    for subID, error in failures:
        print(f"  failed: subID {subID}: {error}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
