#!/usr/bin/env python
"""
Group/cluster-scale view of a TNG FOF group: gas, dark matter and ISM/CGM/IGM.

`SubhaloDataProducts` (and the tutorial notebook) show *one* subhalo. This
script shows the whole FOF group that subhalo belongs to, using the group
catalog's own membership: every gas cell and dark-matter particle is assigned
to the subhalo that owns it, or to the group's unbound fuzz if none does, and
that assignment is what splits the gas into

    ISM   bound to a subhalo AND within k x R_ISM of its centre
    CGM   bound to a subhalo but outside that radius
    IGM   in the group but bound to no subhalo (the fuzz)

The products come from `astrovodka.GroupDataProducts`; this script only plots
them. Three figures: the raw distributions, the three classes, and a set of
one-dimensional diagnostics.

Usage
-----
    export TNG_BASE_PATH=/path/to/TNG/TNG-100
    python examples/en/group_cluster_viz.py --groupID 59 --ncpu 8

A cheap group to try first is gid 14258 (15 subhalos, ~1e4 gas cells, seconds);
gid 59 is the sensible showcase (1710 subhalos, 3.0e6 gas cells, ~20 s with
--ncpu 8). The most massive TNG-100 groups have ~4.5e7 gas cells: pass
--ngrid 128 and expect several GB of reads.
"""

import argparse
import os
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

from astrovodka import GroupDataProducts  # noqa: E402
from astrovodka.group_products import ism_radius  # noqa: E402

ISM, CGM, IGM = 0, 1, 2
CLASS_COLORS = {ISM: "#c0392b", CGM: "#2471a3", IGM: "#7f8c8d"}
CLASS_NAMES = {ISM: "ISM", CGM: "CGM", IGM: "IGM"}


def default_outdir():
    """`examples/outputs/`, resolved against this file rather than the cwd."""
    outdir = Path(__file__).resolve().parent.parent / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def kpc_per_ckpc(dp):
    """Physical kpc per comoving ckpc/h, for this snapshot."""
    return float(dp.metadata["hubble_h"]) * float(dp.metadata["scale_factor"])


def physical_edges(dp):
    """Grid edges in physical kpc (the product stores comoving ckpc/h)."""
    s = kpc_per_ckpc(dp)
    return [np.asarray(e) / s for e in dp.grid_edges]


def project(arr, los_axis):
    """Sum a 3D grid along one axis, ready for imshow.

    Summing leaves the two sky axes in ascending order, so the result's rows are
    the lower-numbered one and its columns the higher-numbered one -- which is
    what `sky_axes` returns as ``(vertical, horizontal)``.
    """
    return arr.sum(axis=los_axis)


def sky_axes(los_axis):
    """The ``(vertical, horizontal)`` sky axes for a projection along ``los``."""
    p, q = [i for i in range(3) if i != los_axis]
    return p, q


def extent_of(edges, los_axis, centre):
    """imshow extent ``[x0, x1, y0, y1]``, as offsets from the group centre."""
    p, q = sky_axes(los_axis)
    return [
        edges[q][0] - centre[q],
        edges[q][-1] - centre[q],
        edges[p][0] - centre[p],
        edges[p][-1] - centre[p],
    ]


def zoom_maps(dp, centre, halfsize, npix=256):
    """Surface-density maps of a sub-volume, gridded from the per-cell arrays.

    A single group-wide grid cannot show galaxies: group 59 spans 5.7 Mpc, so
    even at ngrid = 512 a cell is 11 kpc and a galactic disc is a few pixels.
    Re-gridding just the cells inside the zoom box is cheap and does not
    require the snapshot again, which is why the products keep per-cell arrays.
    """
    cells = dp.cells
    xyz = np.asarray(cells["coords_ckpc_h"], dtype="float64")
    mass = np.asarray(cells["m_gas_msun"], dtype="float64")
    klass = np.asarray(cells["class_id"])
    inside = np.all(np.abs(xyz - centre) < halfsize, axis=1)
    edges = [
        np.linspace(centre[i] - halfsize, centre[i] + halfsize, npix + 1) for i in range(2)
    ]
    area = (edges[0][1] - edges[0][0]) * (edges[1][1] - edges[1][0])
    masks = {"ism": klass == ISM, "cgm": klass == CGM, "all": np.ones(len(xyz), bool)}
    out = {}
    for name, sel in masks.items():
        sel = sel & inside
        h, _, _ = np.histogram2d(xyz[sel, 0], xyz[sel, 1], bins=edges, weights=mass[sel])
        out[name] = h / area
    return out, edges


def label_axes(ax, los_axis):
    """Label the sky axes; the maps are drawn as offsets from the centre."""
    p, q = sky_axes(los_axis)
    ax.set_xlabel(rf"$\Delta${'xyz'[q]} [kpc]")
    ax.set_ylabel(rf"$\Delta${'xyz'[p]} [kpc]")


def overlay_group(ax, dp, los_axis):
    """R200c and the gas-bearing subhalos, projected like the map under them."""
    s = kpc_per_ckpc(dp)
    centre = np.asarray(dp.metadata["group_pos_ckpc_h"]) / s
    a, b = sky_axes(los_axis)
    ax.add_patch(
        Circle(
            (0.0, 0.0),  # the maps are in offsets from the group centre
            dp.metadata["R200c_ckpc_h"] / s,
            fill=False,
            ec="cyan",
            lw=1.3,
            ls="--",
            label=r"$R_{200c}$",
        )
    )
    sel = np.asarray(dp.subhalos["len_type"])[:, 0] > 0
    pos = np.asarray(dp.subhalos["pos"])[sel] / s - centre
    mgas = np.asarray(dp.subhalos["mass_type"])[sel, 0]
    size = 8.0 + 30.0 * np.log10(np.maximum(mgas, 1.0)) / max(np.log10(mgas.max()), 1.0)
    ax.scatter(pos[:, a], pos[:, b], s=size, facecolor="none", edgecolor="cyan", lw=0.6)


def plot_projections(dp, outdir, tag):
    """Dark matter, total gas and HI, projected along each of the three axes."""
    e = physical_edges(dp)
    s_ckpc = kpc_per_ckpc(dp)
    centre = np.asarray(dp.metadata["group_pos_ckpc_h"]) / s_ckpc
    rows = [
        ("dm_real", "Dark matter", "bone"),
        ("gas_real", "Gas (total)", "viridis"),
        ("hi_real", "HI", "inferno"),
    ]
    fig, axes = plt.subplots(len(rows), 3, figsize=(15.5, 4.6 * len(rows)))
    for r, (name, title, cmap) in enumerate(rows):
        arr = getattr(dp, name)
        if arr is None:
            continue
        vmax_hi = float(np.percentile(arr[arr > 0], 99.5)) if (arr > 0).any() else 1.0
        for c, los in enumerate((2, 1, 0)):
            ax = axes[r, c]
            proj = project(arr, los)
            pos = proj[proj > 0]
            im = ax.imshow(
                proj,
                origin="lower",
                extent=extent_of(e, los, centre),
                cmap=cmap,
                norm=LogNorm(vmin=pos.min() if pos.size else 1e-6, vmax=vmax_hi),
            )
            ax.set_title(f"{title} — {'xyz'[los]} projection", fontsize=11)
            label_axes(ax, los)
            if name == "gas_real" and c == 0:
                overlay_group(ax, dp, los)
                ax.legend(loc="upper right", fontsize=8)
            plt.colorbar(im, ax=ax, shrink=0.85, label=r"$M_\odot/\mathrm{kpc}^3$")
    fig.suptitle(
        f"TNG group {dp.metadata['groupID']} (snap {dp.metadata['snapNum']}): "
        f"{dp.metadata['Nsubs']} subhalos, "
        f"R200c = {dp.metadata['R200c_kpc']:.0f} kpc",
    )
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outdir / f"{tag}_projections.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_classes(dp, outdir, tag):
    """The three classes projected, plus a zoom on the central galaxies.

    Two scales are needed: a FOF group reaches several Mpc (it percolates along
    filaments), while its galaxies are tens of kpc across, so the ISM is
    invisible in the wide view.
    """
    e = physical_edges(dp)
    s = kpc_per_ckpc(dp)
    centre = np.asarray(dp.metadata["group_pos_ckpc_h"]) / s
    r200 = dp.metadata["R200c_ckpc_h"] / s

    r_ism = ism_radius(
        np.asarray(dp.subhalos["rhalf_type"])[:, 0],
        np.asarray(dp.subhalos["rhalf_type"])[:, 4],
        dp.metadata["ism_k_gas"],
        dp.metadata["ism_k_star"],
    )
    fig, axes = plt.subplots(4, 3, figsize=(15.5, 4.4 * 4))
    for r, (name, klass) in enumerate(
        [("gas_ism_real", ISM), ("gas_cgm_real", CGM), ("gas_igm_real", IGM)]
    ):
        arr = getattr(dp, name)
        vmax = float(np.percentile(arr[arr > 0], 99.5)) if (arr > 0).any() else 1.0
        for c, los in enumerate((2, 1, 0)):
            ax = axes[r, c]
            proj = project(arr, los)
            pos = proj[proj > 0]
            im = ax.imshow(
                proj,
                origin="lower",
                extent=extent_of(e, los, centre),
                cmap="magma",
                norm=LogNorm(vmin=pos.min() if pos.size else 1e-6, vmax=vmax),
            )
            ax.set_title(f"{CLASS_NAMES[klass]} — {'xyz'[los]} projection", fontsize=11)
            label_axes(ax, los)
            if c == 0:
                overlay_group(ax, dp, los)
            plt.colorbar(im, ax=ax, shrink=0.85, label=r"$M_\odot/\mathrm{kpc}^3$")

    # Bottom row: the brightest galaxy, on its own high-resolution grid. The
    # group-wide grid is ~30 kpc per cell, so zooming into it would show blocks;
    # the per-cell arrays let this be re-gridded at a galactic scale instead.
    gal = int(np.argmax(np.asarray(dp.subhalos["mass_type"])[:, 0]))
    gcentre = np.asarray(dp.subhalos["pos"])[gal] / s
    halfsize = float(np.clip(4.0 * r_ism[gal] / s, 60.0, 0.3 * r200))
    maps, zedges = zoom_maps(dp, gcentre * s, halfsize * s, npix=256)
    zext = [
        zedges[0][0] / s - gcentre[0],
        zedges[0][-1] / s - gcentre[0],
        zedges[1][0] / s - gcentre[1],
        zedges[1][-1] / s - gcentre[1],
    ]
    normed = {}
    for key in ("ism", "cgm"):
        d = maps[key]
        pos = d[d > 0]
        lo = float(np.percentile(pos, 25)) if pos.size else 1e-6
        hi = float(np.percentile(pos, 99.5)) if pos.size else 1.0
        normed[key] = np.clip(
            (np.log10(np.maximum(d, lo)) - np.log10(lo))
            / max(np.log10(max(hi, lo * 1.01)) - np.log10(lo), 1e-9),
            0.0,
            1.0,
        )
        normed[key + "_range"] = (lo, hi)
    for c, (key, title) in enumerate(
        [
            ("ism", "ISM of the brightest galaxy"),
            ("cgm", "CGM around the same galaxy"),
            (None, "both (red = ISM, blue = CGM)"),
        ]
    ):
        ax = axes[3, c]
        if key is None:
            img = np.zeros(normed["ism"].T.shape + (3,))
            img[..., 0] = normed["ism"].T
            img[..., 2] = normed["cgm"].T
            ax.imshow(img, origin="lower", extent=zext)
        else:
            lo, hi = normed[key + "_range"]
            im = ax.imshow(
                maps[key].T, origin="lower", extent=zext, cmap="magma",
                norm=LogNorm(vmin=lo, vmax=hi),
            )
            plt.colorbar(im, ax=ax, shrink=0.85, label=r"$M_\odot/\mathrm{kpc}^2$")
        ax.set_title(f"{title} (±{halfsize:.0f} kpc)", fontsize=11)
        ax.set_xlabel(r"$\Delta$x [kpc]")
        ax.set_ylabel(r"$\Delta$y [kpc]")
    fig.suptitle(
        f"ISM (bound to a subhalo and within {dp.metadata['ism_k_gas']:g}xR_half,gas "
        f"or {dp.metadata['ism_k_star']:g}xR_half,star of it) / CGM / IGM",
    )
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outdir / f"{tag}_ism_cgm_igm.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_rsd(dp, outdir, tag):
    """Real space against redshift space, for the gas and for the HI.

    The classes are assigned in real space and only then transformed, so what
    this figure shows is purely a change of coordinates: the same cells, moved
    along the line of sight by ``v_los / (a H(z))`` about the group's bulk
    velocity. Both maps are drawn on the redshift-space box, which is why the
    real-space one is a thin slab -- that difference *is* the effect, and it is
    large: at z=0 one km/s is 10 ckpc/h, so a few hundred km/s of dispersion
    stretches the group over several Mpc.
    """
    if dp.grid_edges_rsd is None:
        print("  (no redshift-space grids: pass compute_rsd=True)")
        return False
    s = kpc_per_ckpc(dp)
    centre = np.asarray(dp.metadata["group_pos_ckpc_h"]) / s
    e_real = physical_edges(dp)
    e_rsd = [np.asarray(a) / s for a in dp.grid_edges_rsd]
    los = 0  # x is the line of sight
    p, q = sky_axes(los)

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.0))
    for r, (name, title) in enumerate([("gas", "Gas"), ("hi", "HI")]):
        for c, (space, edges) in enumerate([("real", e_real), ("rsd", e_rsd)]):
            space_label = "real space" if space == "real" else "redshift space"
            ax = axes[r, c]
            arr = getattr(dp, f"{name}_{space}")
            proj = project(arr, los)
            pos = proj[proj > 0]
            im = ax.imshow(
                proj,
                origin="lower",
                extent=extent_of(edges, los, centre),
                cmap="inferno",
                aspect="auto",
                norm=LogNorm(vmin=pos.min() if pos.size else 1e-6,
                             vmax=float(np.percentile(pos, 99.9)) if pos.size else 1.0),
            )
            ax.set_title(f"{title} — {space_label}", fontsize=11)
            ax.set_xlabel(r"$\Delta$x [kpc] (line of sight)")
            ax.set_ylabel(rf"$\Delta${'xyz'[p]} [kpc]")
            plt.colorbar(im, ax=ax, shrink=0.8)

        # Line-of-sight profile: the elongation quantified.
        ax = axes[r, 2]
        for space, edges, colour in (
            ("real", e_real, "#2471a3"),
            ("rsd", e_rsd, "#c0392b"),
        ):
            arr = getattr(dp, f"{name}_{space}")
            prof = arr.sum(axis=(1, 2))
            c_ax = 0.5 * (edges[0][1:] + edges[0][:-1]) - centre[0]
            ax.step(c_ax, prof, where="mid", color=colour, lw=1.6,
                    label=space_label)
        ax.set_yscale("log")
        ax.set_xlabel(r"$\Delta$x [kpc] (line of sight)")
        ax.set_ylabel(r"surface density [$M_\odot$/kpc$^2$]")
        ax.set_title(f"{title} along the line of sight", fontsize=11)
        ax.legend()
    fig.suptitle(
        f"group {dp.metadata['groupID']}: classification is done in real space, "
        "then the cells are moved along the line of sight",
    )
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outdir / f"{tag}_rsd.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return True


def plot_diagnostics(dp, outdir, tag):
    """Radial profiles, the phase diagram, the budget and the HI inventory."""
    cells = dp.cells
    md = dp.metadata
    r = np.asarray(cells["r_grp_over_r200c"], dtype="float64")
    m_gas = np.asarray(cells["m_gas_msun"], dtype="float64")
    m_hi = np.asarray(cells["m_hi_msun"], dtype="float64")
    class_id = np.asarray(cells["class_id"])
    nH = np.asarray(cells["nH_cm3"], dtype="float64")
    T = np.asarray(cells["T_K"], dtype="float64")

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.5))

    # (a) enclosed gas mass by class
    ax = axes[0, 0]
    bins = np.linspace(0, max(float(np.percentile(r, 99)), 1.2), 60)
    for klass in (ISM, CGM, IGM):
        sel = class_id == klass
        if not sel.any():
            continue
        h, edges = np.histogram(r[sel], bins=bins, weights=m_gas[sel])
        ax.plot(0.5 * (edges[1:] + edges[:-1]), np.cumsum(h),
                color=CLASS_COLORS[klass], label=CLASS_NAMES[klass], lw=2)
    ax.set_yscale("log")
    ax.axvline(1.0, color="k", ls=":", lw=1)
    ax.set_xlabel(r"$r / R_{200c}$")
    ax.set_ylabel(r"enclosed gas mass [$M_\odot$]")
    ax.set_title("(a) gas mass by class, from the group centre")
    ax.legend()

    # (b) phase diagram, one point per gas cell
    ax = axes[0, 1]
    for klass in (IGM, CGM, ISM):
        sel = class_id == klass
        if not sel.any():
            continue
        ax.scatter(nH[sel], T[sel], s=0.4, alpha=0.2, lw=0,
                   color=CLASS_COLORS[klass], label=CLASS_NAMES[klass], rasterized=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axvline(0.1, color="k", ls="--", lw=1)
    ax.set_xlabel(r"$n_H$ [cm$^{-3}$]")
    ax.set_ylabel("T [K]")
    ax.set_title("(b) phase diagram, one point per gas cell")
    # Proxy handles rather than the scatter artists: the cells are drawn nearly
    # transparent, so their own legend entries would be invisible.
    ax.legend(
        handles=[
            Line2D([], [], marker="o", ls="", color=CLASS_COLORS[k],
                   label=CLASS_NAMES[k], markersize=8)
            for k in (ISM, CGM, IGM)
        ],
        loc="lower right",
    )

    # (c) mass budget: what fraction of the gas and of the HI is in each class
    ax = axes[1, 0]
    x = np.arange(3)
    gas_frac = [md[f"m_gas_{n}_msun"] / md["m_gas_total_msun"] for n in ("ism", "cgm", "igm")]
    hi_frac = [md[f"m_hi_{n}_msun"] / max(md["m_hi_total_msun"], 1e-30)
               for n in ("ism", "cgm", "igm")]
    ax.bar(x - 0.2, gas_frac, 0.4, label="gas mass", color="#5d6d7e")
    ax.bar(x + 0.2, hi_frac, 0.4, label="HI mass", color="#d68910")
    for i, (g, h) in enumerate(zip(gas_frac, hi_frac)):
        ax.text(i - 0.2, g + 0.02, f"{100 * g:.1f}%", ha="center", fontsize=9)
        ax.text(i + 0.2, h + 0.02, f"{100 * h:.1f}%", ha="center", fontsize=9)
    ax.set_xticks(x, ["ISM", "CGM", "IGM"])
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("fraction of the group total")
    ax.set_title("(c) mass budget")
    ax.legend()

    # (d) HI against gas mass, per subhalo
    ax = axes[1, 1]
    owner = np.asarray(cells["owner"])
    sel = owner >= 0
    n_sub = len(dp.subhalos["len_type"])
    m_gas_sub = np.bincount(owner[sel], weights=m_gas[sel], minlength=n_sub)
    m_hi_sub = np.bincount(owner[sel], weights=m_hi[sel], minlength=n_sub)
    ok = (m_gas_sub > 0) & (m_hi_sub > 0)
    ax.scatter(m_gas_sub[ok], m_hi_sub[ok], s=24, color="#2471a3", alpha=0.85,
               edgecolor="k", lw=0.3)
    lim = [max(m_gas_sub[ok].min(), 1e6) * 0.5, m_gas_sub[ok].max() * 2]
    ax.plot(lim, lim, "k--", lw=1, label="HI = gas")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_xlabel(r"gas mass in the subhalo [$M_\odot$]")
    ax.set_ylabel(r"HI mass in the subhalo [$M_\odot$]")
    ax.set_title(f"(d) HI inventory ({int(ok.sum())} gas-bearing subhalos)")
    ax.legend()

    fig.suptitle(
        f"group {md['groupID']}: M200c = {md['M200c_msun']:.2e} Msun, "
        f"M_HI = {md['m_hi_total_msun']:.2e} Msun, model {md['hi_model']}",
    )
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outdir / f"{tag}_diagnostics.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def summarise(dp):
    """Print the numbers a reader of the figures would want."""
    md = dp.metadata
    print(f"\nGroup {md['groupID']} (snapshot {md['snapNum']}, z={md['redshift']:.4f})")
    print(f"  subhalos           {md['Nsubs']} ({md['N_sub_with_gas']} with gas)")
    print(f"  M200c              {md['M200c_msun']:.3e} Msun")
    print(f"  R200c              {md['R200c_ckpc_h']:.1f} ckpc/h "
          f"= {md['R200c_kpc']:.0f} physical kpc")
    print(f"  gas cells          {md['n_gas']:,}")
    print(f"  DM particles       {md['n_dm']:,}")
    print(f"  box half-width     {[round(float(v), 1) for v in md['box_halfsize_ckpc']]}")
    print(f"  grid               {dp.ngrid}^3")
    print("\n  class        cells        gas mass [Msun]     HI mass [Msun]    f_HI")
    for name in ("ism", "cgm", "igm"):
        g = md[f"m_gas_{name}_msun"]
        hi = md[f"m_hi_{name}_msun"]
        print(f"  {name.upper():<6} {md['n_cells_' + name]:>12,}  {g:>19.4e}  "
              f"{hi:>18.4e}   {hi / max(g, 1):.4f}")
    print(f"  {'total':<6} {md['n_gas']:>12,}  {md['m_gas_total_msun']:>19.4e}  "
          f"{md['m_hi_total_msun']:>18.4e}")
    print(f"\n  M_H2               {md['m_h2_total_msun']:.4e} Msun")
    print(f"  M_DM               {md['m_dm_total_msun']:.4e} Msun "
          f"({100 * md['m_dm_igm_total_msun'] / md['m_dm_total_msun']:.2f}% in the IGM)")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--basePath", default=os.environ.get("TNG_BASE_PATH"),
                        help="TNG output directory (default: $TNG_BASE_PATH)")
    parser.add_argument("--snapNum", type=int, default=99)
    parser.add_argument("--groupID", type=int, default=59)
    parser.add_argument("--ngrid", type=int, default=256)
    parser.add_argument("--hiModel", default="K13_vol",
                        choices=["K13_vol", "GD14_vol", "none"])
    parser.add_argument("--kGas", type=float, default=1.0,
                        help="ISM radius multiplier on R_half,gas")
    parser.add_argument("--kStar", type=float, default=2.0,
                        help="ISM radius multiplier on R_half,star")
    parser.add_argument("--ncpu", type=int, default=1,
                        help="processes for the per-subhalo HI loop")
    parser.add_argument("--noRsd", action="store_true",
                        help="skip the redshift-space grids and their figure")
    parser.add_argument("--outdir", default=None)
    args = parser.parse_args()
    if not args.basePath:
        parser.error("no TNG data directory: pass --basePath or set $TNG_BASE_PATH")
    basePath = args.basePath.rstrip("/") + "/"
    outdir = Path(args.outdir) if args.outdir else default_outdir()
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading group {args.groupID} from snapshot {args.snapNum} ...")
    dp = GroupDataProducts.from_group(
        basePath,
        args.snapNum,
        args.groupID,
        ngrid=args.ngrid,
        hi_model=args.hiModel,
        compute_hi=args.hiModel != "none",
        k_gas=args.kGas,
        k_star=args.kStar,
        ncpu=args.ncpu,
        compute_rsd=not args.noRsd,
    )
    summarise(dp)

    tag = f"group{args.groupID:07d}_snap{args.snapNum:03d}"
    h5 = outdir / f"{tag}_data_products.hdf5"
    dp.save_hdf5(h5)
    print(f"\nSaved {h5}")

    plot_projections(dp, outdir, tag)
    plot_classes(dp, outdir, tag)
    has_rsd = plot_rsd(dp, outdir, tag)
    plot_diagnostics(dp, outdir, tag)
    print(f"\nFigures written to {outdir}:")
    written = ["projections", "ism_cgm_igm"]
    if has_rsd:
        written.append("rsd")
    written.append("diagnostics")
    for name in written:
        print(f"  {tag}_{name}.png")


if __name__ == "__main__":
    main()
