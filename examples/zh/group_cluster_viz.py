#!/usr/bin/env python
"""
TNG 星系群/星系团尺度（FOF group）的气体、暗物质与 ISM/CGM/IGM 分类图。

`SubhaloDataProducts`（以及教程 notebook）画的是*单个* subhalo。这个脚本画它所属的
整个 FOF group，用的是 group catalog 自带的归属关系：每个气体胞元和暗物质粒子都归给
拥有它的 subhalo，若不属于任何 subhalo 则算作 group 的 fuzz。正是这个归属把气体分成

    ISM   归属于某个 subhalo，且距该星系中心 < k x R_ISM
    CGM   归属于某个 subhalo，但在该半径之外
    IGM   在 group 内但不属于任何 subhalo（即 fuzz）

数据来自 `astrovodka.GroupDataProducts`，本脚本只负责画图。三张图：原始分布、
三个类别、以及一组一维诊断。

用法
----
    export TNG_BASE_PATH=/path/to/TNG/TNG-100
    python examples/zh/group_cluster_viz.py --groupID 59 --ncpu 8

想快速试跑可以用 gid 14258（15 个 subhalo、约 1e4 个气体胞元，秒级）；gid 59 适合
正式出图（1710 个 subhalo、3.0e6 个气体胞元，--ncpu 8 约 20 秒）。最大的那些
TNG-100 group 有约 4.5e7 个气体胞元：请用 --ngrid 128，并预留数 GB 内存。
"""

import argparse
import os
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

from astrovodka import GroupDataProducts  # noqa: E402
from astrovodka.group_products import ism_radius  # noqa: E402

ISM, CGM, IGM = 0, 1, 2
CLASS_COLORS = {ISM: "#c0392b", CGM: "#2471a3", IGM: "#7f8c8d"}
CLASS_NAMES = {ISM: "ISM", CGM: "CGM", IGM: "IGM"}


def setup_cjk_font():
    """图里的中文标签需要 CJK 字体，否则 matplotlib 会画成方框（tofu）。"""
    available = {f.name for f in matplotlib.font_manager.fontManager.ttflist}
    for candidate in ("Noto Sans CJK JP", "Noto Sans CJK SC", "Source Han Sans SC",
                      "Droid Sans Fallback", "WenQuanYi Zen Hei"):
        if candidate in available:
            plt.rcParams["font.sans-serif"] = [candidate] + list(
                plt.rcParams["font.sans-serif"]
            )
            print(f"CJK 字体: {candidate}")
            break
    else:
        print("警告：没找到 CJK 字体，图上的中文会显示成方框")
    plt.rcParams["axes.unicode_minus"] = False   # 负号也走同一套字体


def default_outdir():
    """`examples/outputs/`，按本文件的位置解析，而不是当前工作目录。"""
    outdir = Path(__file__).resolve().parent.parent / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def kpc_per_ckpc(dp):
    """该快照下 1 ckpc/h 对应多少物理 kpc。"""
    return float(dp.metadata["hubble_h"]) * float(dp.metadata["scale_factor"])


def physical_edges(dp):
    """网格边界换算到物理 kpc（产品里存的是共动 ckpc/h）。"""
    s = kpc_per_ckpc(dp)
    return [np.asarray(e) / s for e in dp.grid_edges]


def project(arr, los_axis):
    """沿某个轴求和，得到可直接交给 imshow 的二维数组。

    求和后剩下的两个天球轴按序号升序排列，因此结果的「行」是序号较小的那个、
    「列」是较大的那个——这正是 `sky_axes` 返回的 ``(竖直, 水平)``。
    """
    return arr.sum(axis=los_axis)


def sky_axes(los_axis):
    """沿 ``los_axis`` 投影后剩下的 ``(竖直, 水平)`` 两轴。"""
    p, q = [i for i in range(3) if i != los_axis]
    return p, q


def extent_of(edges, los_axis, centre):
    """imshow 的 extent ``[x0, x1, y0, y1]``，用相对 group 中心的偏移。"""
    p, q = sky_axes(los_axis)
    return [
        edges[q][0] - centre[q],
        edges[q][-1] - centre[q],
        edges[p][0] - centre[p],
        edges[p][-1] - centre[p],
    ]


def zoom_maps(dp, centre, halfsize, npix=256):
    """用小体积内的胞元重新做一个高分辨率的投影面密度图。

    单一 group 尺度的网格画不出星系：gid 59 跨度 5.7 Mpc，即使 ngrid = 512
    一格也有 11 kpc，星系盘只剩几个像素。把缩放框内的胞元重新分箱既便宜，也不需要
    再读快照——产品里保留逐胞元数组正是为了这个。
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
    """标注天球面两轴；图是以 group 中心为零点的偏移。"""
    p, q = sky_axes(los_axis)
    ax.set_xlabel(rf"$\Delta${'xyz'[q]} [kpc]")
    ax.set_ylabel(rf"$\Delta${'xyz'[p]} [kpc]")


def overlay_group(ax, dp, los_axis):
    """叠加 R200c 与含气 subhalo 的位置，投影方式与底图一致。"""
    s = kpc_per_ckpc(dp)
    centre = np.asarray(dp.metadata["group_pos_ckpc_h"]) / s
    a, b = sky_axes(los_axis)
    ax.add_patch(
        Circle(
            (0.0, 0.0),  # 图本身就以 group 中心为原点
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
    """暗物质、总气体与 HI，沿三个轴各做一次投影。"""
    e = physical_edges(dp)
    s_ckpc = kpc_per_ckpc(dp)
    centre = np.asarray(dp.metadata["group_pos_ckpc_h"]) / s_ckpc
    rows = [
        ("dm_real", "暗物质", "bone"),
        ("gas_real", "气体（全部）", "viridis"),
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
            ax.set_title(f"{title} — 沿 {'xyz'[los]} 轴投影", fontsize=11)
            label_axes(ax, los)
            if name == "gas_real" and c == 0:
                overlay_group(ax, dp, los)
                ax.legend(loc="upper right", fontsize=8)
            plt.colorbar(im, ax=ax, shrink=0.85, label=r"$M_\odot/\mathrm{kpc}^3$")
    fig.suptitle(
        f"TNG group {dp.metadata['groupID']}（snap {dp.metadata['snapNum']}）："
        f"{dp.metadata['Nsubs']} 个 subhalo，"
        f"R200c = {dp.metadata['R200c_kpc']:.0f} kpc"
    )
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outdir / f"{tag}_projections.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_classes(dp, outdir, tag):
    """三个类别各自的投影，外加最亮星系的高分辨率缩放。

    必须给两个尺度：FOF group 沿纤维结构延伸可达数 Mpc，而它里面的星系只有几十 kpc，
    宽视场里根本看不见 ISM。
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
            ax.set_title(f"{CLASS_NAMES[klass]} — 沿 {'xyz'[los]} 轴投影", fontsize=11)
            label_axes(ax, los)
            if c == 0:
                overlay_group(ax, dp, los)
            plt.colorbar(im, ax=ax, shrink=0.85, label=r"$M_\odot/\mathrm{kpc}^3$")

    # 最下面一行：最亮星系，用它自己的高分辨率网格。整组的网格一格约 30 kpc，
    # 直接放大只会看到方块；逐胞元数组可以按星系尺度重新分箱。
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
            ("ism", "最亮星系的 ISM"),
            ("cgm", "同一个星系周围的 CGM"),
            (None, "两者叠加（红=ISM，蓝=CGM）"),
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
        ax.set_title(f"{title}（±{halfsize:.0f} kpc）", fontsize=11)
        ax.set_xlabel(r"$\Delta$x [kpc]")
        ax.set_ylabel(r"$\Delta$y [kpc]")
    fig.suptitle(
        f"ISM（归属于 subhalo 且在其 {dp.metadata['ism_k_gas']:g}xR_half,gas 或 "
        f"{dp.metadata['ism_k_star']:g}xR_half,star 之内）/ CGM / IGM"
    )
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outdir / f"{tag}_ism_cgm_igm.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_rsd(dp, outdir, tag):
    """实空间与红移空间对照，气体与 HI 各一行。

    类别是在实空间中判定、之后才做变换的，所以这张图展示的纯粹是坐标变换：
    同一批胞元，沿视线方向按 ``v_los / (a H(z))``（相对 group 整体运动）平移。
    两张图都画在红移空间的盒子上，因此实空间那张是一条薄片——这个差别本身就是效应，
    而且很大：z=0 时 1 km/s 就是 10 ckpc/h，几百 km/s 的速度弥散会把 group 拉长
    好几个 Mpc。
    """
    if dp.grid_edges_rsd is None:
        print("  （没有红移空间网格：请传 compute_rsd=True）")
        return False
    s = kpc_per_ckpc(dp)
    centre = np.asarray(dp.metadata["group_pos_ckpc_h"]) / s
    e_real = physical_edges(dp)
    e_rsd = [np.asarray(a) / s for a in dp.grid_edges_rsd]
    los = 0  # x 是视线方向
    p, q = sky_axes(los)

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.0))
    for r, (name, title) in enumerate([("gas", "气体"), ("hi", "HI")]):
        for c, (space, edges) in enumerate([("real", e_real), ("rsd", e_rsd)]):
            space_label = "实空间" if space == "real" else "红移空间"
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
            ax.set_xlabel(r"$\Delta$x [kpc]（视线方向）")
            ax.set_ylabel(rf"$\Delta${'xyz'[p]} [kpc]")
            plt.colorbar(im, ax=ax, shrink=0.8)

        # 视线方向剖面：把拉长量化出来。
        ax = axes[r, 2]
        for space, edges, colour in (("real", e_real, "#2471a3"),
                                     ("rsd", e_rsd, "#c0392b")):
            arr = getattr(dp, f"{name}_{space}")
            prof = arr.sum(axis=(1, 2))
            c_ax = 0.5 * (edges[0][1:] + edges[0][:-1]) - centre[0]
            ax.step(c_ax, prof, where="mid", color=colour, lw=1.6,
                    label=space_label)
        ax.set_yscale("log")
        ax.set_xlabel(r"$\Delta$x [kpc]（视线方向）")
        ax.set_ylabel(r"面密度 [$M_\odot$/kpc$^2$]")
        ax.set_title(f"{title} 沿视线方向的分布", fontsize=11)
        ax.legend()
    fig.suptitle(
        f"group {dp.metadata['groupID']}：分类在实空间中完成，之后才把胞元沿视线方向平移"
    )
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outdir / f"{tag}_rsd.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return True


def plot_diagnostics(dp, outdir, tag):
    """径向剖面、相图、质量收支与 HI 清单。"""
    cells = dp.cells
    md = dp.metadata
    r = np.asarray(cells["r_grp_over_r200c"], dtype="float64")
    m_gas = np.asarray(cells["m_gas_msun"], dtype="float64")
    m_hi = np.asarray(cells["m_hi_msun"], dtype="float64")
    class_id = np.asarray(cells["class_id"])
    nH = np.asarray(cells["nH_cm3"], dtype="float64")
    T = np.asarray(cells["T_K"], dtype="float64")

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.5))

    # (a) 按类别累计的气体质量
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
    ax.set_ylabel(r"累计气体质量 [$M_\odot$]")
    ax.set_title("(a) 自 group 中心起、按类别累计的气体质量")
    ax.legend()

    # (b) 相图，一个气体胞元一个点
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
    ax.set_title("(b) 相图（每个气体胞元一个点）")
    # 用代理图例而不是散点本身：胞元画得几乎透明，它们自己的图例项会看不见。
    ax.legend(
        handles=[
            Line2D([], [], marker="o", ls="", color=CLASS_COLORS[k],
                   label=CLASS_NAMES[k], markersize=8)
            for k in (ISM, CGM, IGM)
        ],
        loc="lower right",
    )

    # (c) 质量收支：每一类占气体质量和 HI 质量的比例
    ax = axes[1, 0]
    x = np.arange(3)
    gas_frac = [md[f"m_gas_{n}_msun"] / md["m_gas_total_msun"] for n in ("ism", "cgm", "igm")]
    hi_frac = [md[f"m_hi_{n}_msun"] / max(md["m_hi_total_msun"], 1e-30)
               for n in ("ism", "cgm", "igm")]
    ax.bar(x - 0.2, gas_frac, 0.4, label="气体质量", color="#5d6d7e")
    ax.bar(x + 0.2, hi_frac, 0.4, label="HI 质量", color="#d68910")
    for i, (g, h) in enumerate(zip(gas_frac, hi_frac)):
        ax.text(i - 0.2, g + 0.02, f"{100 * g:.1f}%", ha="center", fontsize=9)
        ax.text(i + 0.2, h + 0.02, f"{100 * h:.1f}%", ha="center", fontsize=9)
    ax.set_xticks(x, ["ISM", "CGM", "IGM"])
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("占 group 总量的比例")
    ax.set_title("(c) 质量收支")
    ax.legend()

    # (d) 每个 subhalo 的 HI 与其气体质量
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
    ax.plot(lim, lim, "k--", lw=1, label="HI = 气体质量")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_xlabel(r"subhalo 内的气体质量 [$M_\odot$]")
    ax.set_ylabel(r"subhalo 内的 HI 质量 [$M_\odot$]")
    ax.set_title(f"(d) HI 清单（{int(ok.sum())} 个含气 subhalo）")
    ax.legend()

    fig.suptitle(
        f"group {md['groupID']}：M200c = {md['M200c_msun']:.2e} Msun，"
        f"M_HI = {md['m_hi_total_msun']:.2e} Msun，模型 {md['hi_model']}"
    )
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outdir / f"{tag}_diagnostics.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def summarise(dp):
    """打印读图时想知道的那些数。"""
    md = dp.metadata
    print(f"\nGroup {md['groupID']}（快照 {md['snapNum']}，z={md['redshift']:.4f}）")
    print(f"  subhalo 数目       {md['Nsubs']}（其中 {md['N_sub_with_gas']} 个含气体）")
    print(f"  M200c              {md['M200c_msun']:.3e} Msun")
    print(f"  R200c              {md['R200c_ckpc_h']:.1f} ckpc/h "
          f"= {md['R200c_kpc']:.0f} 物理 kpc")
    print(f"  气体胞元           {md['n_gas']:,}")
    print(f"  暗物质粒子         {md['n_dm']:,}")
    print(f"  盒子半宽           {[round(float(v), 1) for v in md['box_halfsize_ckpc']]}")
    print(f"  网格               {dp.ngrid}^3")
    print("\n  类别        胞元数        气体质量 [Msun]      HI 质量 [Msun]    f_HI")
    for name in ("ism", "cgm", "igm"):
        g = md[f"m_gas_{name}_msun"]
        hi = md[f"m_hi_{name}_msun"]
        print(f"  {name.upper():<7} {md['n_cells_' + name]:>12,}  {g:>19.4e}  "
              f"{hi:>18.4e}   {hi / max(g, 1):.4f}")
    print(f"  {'合计':<6} {md['n_gas']:>12,}  {md['m_gas_total_msun']:>19.4e}  "
          f"{md['m_hi_total_msun']:>18.4e}")
    print(f"\n  M_H2               {md['m_h2_total_msun']:.4e} Msun")
    print(f"  M_DM               {md['m_dm_total_msun']:.4e} Msun "
          f"（其中 {100 * md['m_dm_igm_total_msun'] / md['m_dm_total_msun']:.2f}% 在 IGM）")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--basePath", default=os.environ.get("TNG_BASE_PATH"),
                        help="TNG 输出目录（默认取 $TNG_BASE_PATH）")
    parser.add_argument("--snapNum", type=int, default=99)
    parser.add_argument("--groupID", type=int, default=59)
    parser.add_argument("--ngrid", type=int, default=256)
    parser.add_argument("--hiModel", default="K13_vol",
                        choices=["K13_vol", "GD14_vol", "none"])
    parser.add_argument("--kGas", type=float, default=1.0,
                        help="ISM 半径中 R_half,gas 的系数")
    parser.add_argument("--kStar", type=float, default=2.0,
                        help="ISM 半径中 R_half,star 的系数")
    parser.add_argument("--ncpu", type=int, default=1,
                        help="逐 subhalo 算 HI 时使用的进程数")
    parser.add_argument("--noRsd", action="store_true",
                        help="跳过红移空间网格及其配图")
    parser.add_argument("--outdir", default=None)
    args = parser.parse_args()
    if not args.basePath:
        parser.error("没有 TNG 数据目录：请传 --basePath 或设置 $TNG_BASE_PATH")
    basePath = args.basePath.rstrip("/") + "/"
    outdir = Path(args.outdir) if args.outdir else default_outdir()
    outdir.mkdir(parents=True, exist_ok=True)
    setup_cjk_font()

    print(f"正在读取 group {args.groupID}（快照 {args.snapNum}）...")
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
    print(f"\n已保存 {h5}")

    plot_projections(dp, outdir, tag)
    plot_classes(dp, outdir, tag)
    has_rsd = plot_rsd(dp, outdir, tag)
    plot_diagnostics(dp, outdir, tag)
    print(f"\n图已写入 {outdir}：")
    written = ["projections", "ism_cgm_igm"]
    if has_rsd:
        written.append("rsd")
    written.append("diagnostics")
    for name in written:
        print(f"  {tag}_{name}.png")


if __name__ == "__main__":
    main()
