#!/usr/bin/env python
"""批量处理 IllustrisTNG 群目录中的子晕。

读取群目录里的子晕，按气体单元数挑选目标，再为每个目标构建一个模拟 HI
数据立方体。当你需要的是一个样本而不是单个星系时，就该用这个模式。

用法
----
    export TNG_BASE_PATH=/path/to/TNG/TNG-100/
    python examples/zh/batch_subhalos.py --limit 20
    python examples/zh/batch_subhalos.py --limit 10 --hiModel K13_map --snapNum 99

输出默认落在 ``examples/outputs/batch/``，除非用 ``--outdir`` 另行指定。

关于子晕选择的说明
------------------
群目录按**质量降序**排列，所以排在最前面的就是团核心：TNG-100 快照 99 的
子晕 0 约有 3800 万个气体单元，读入需要超过 10 GB。因此 ``ids[:N]`` 是个
把登录节点挂死的好办法。

先按气体单元数筛选 —— 这正是 ``--minGas``/``--maxGas`` 所做的事 ——
挑出来的才是真正值得做 HI 工作的低质量、气体丰富卫星星系，同时也限制了
每个目标的内存占用。
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
    """返回 ``examples/outputs/batch/``，相对于本文件解析。"""
    outdir = Path(__file__).resolve().parent.parent / "outputs" / "batch"
    outdir.mkdir(parents=True, exist_ok=True)
    return str(outdir)


def _read_catalog(basePath, snapNum):
    """读取群目录中每个子晕的气体单元数与恒星形成率。

    群目录被切分成数百个分块文件，各分块在子晕 ID 上是连续的，所以按数字顺序
    把它们拼接起来，就还原出全局的 ID 顺序。只读第一个分块是个值得避开的陷阱：
    群目录按质量降序排列，因此第 0 块里只有盒子里质量最大的那些天体，别的什么
    都没有。

    返回值
    -------
    n_gas, sfr : numpy.ndarray
        逐子晕的气体单元数与恒星形成率 [Msun/yr]，以子晕 ID 为索引。
    """
    catalog_dir = Path(basePath) / f"groups_{snapNum:03d}"
    chunks = sorted(
        catalog_dir.glob(f"fof_subhalo_tab_{snapNum:03d}.*.hdf5"),
        key=lambda p: int(p.name.rsplit(".", 2)[-2]),
    )
    if not chunks:
        raise FileNotFoundError(
            f"{catalog_dir} 中没有群目录文件。\n"
            "basePath 必须指向 *包含* groups_<snap>/ 与 snapdir_<snap>/ 的"
            "那一级目录，并且该快照必须已经下载。"
        )

    n_gas_parts, sfr_parts = [], []
    for chunk in chunks:
        with h5py.File(chunk, "r") as handle:
            subhalo = handle.get("Subhalo")
            # 有些分块是占位用的空块：那个 ID 区间里没有子晕，对应的数据集在
            # 文件里根本不存在。它们什么也不贡献，所以直接跳过，不必特殊处理。
            if subhalo is None or "SubhaloLenType" not in subhalo:
                continue
            # SubhaloLenType[:, 0] 就是气体单元数，粒子类型 0 即气体。
            n_gas_parts.append(np.asarray(subhalo["SubhaloLenType"][:, 0]))
            sfr_parts.append(np.asarray(subhalo["SubhaloSFR"]))
    if not n_gas_parts:
        raise RuntimeError(f"在 {catalog_dir} 下的所有分块中都没有找到子晕")
    return np.concatenate(n_gas_parts), np.concatenate(sfr_parts)


def select_subhalos(basePath, snapNum, min_gas, max_gas, limit, min_sfr=0.0):
    """返回气体单元数落在 ``[min_gas, max_gas]`` 区间内的子晕 ID。

    参数
    ----
    basePath : str
        TNG 模拟目录，其中包含 ``groups_<snap>/``。
    snapNum : int
        快照号。
    min_gas, max_gas : int
        气体单元数的闭区间上下限。
    limit : int or None
        最多保留这么多个，优先取气体最丰富的。
    min_sfr : float, optional
        最低恒星形成率 [Msun/yr]。对 HI 工作而言，正在形成恒星的天体才值得
        关注：自适应 FFT UV 场只在有年轻恒星的地方才真正起作用。

    返回值
    -------
    list of tuple
        每个入选目标的 ``(subhalo_id, n_gas, sfr)``。
    """
    n_gas, sfr = _read_catalog(basePath, snapNum)
    print(f"群目录：快照 {snapNum:03d} 中有 {len(n_gas)} 个子晕")

    inside = np.flatnonzero(
        (n_gas >= min_gas) & (n_gas <= max_gas) & (sfr >= min_sfr)
    )
    # 气体最丰富的排在前面，这样 --limit 保留的就是最值得关注的目标。
    inside = inside[np.argsort(n_gas[inside])[::-1]]
    if limit is not None:
        inside = inside[:limit]

    print(
        f"选中 {len(inside)} 个：{min_gas} <= N_gas <= {max_gas}，"
        f"且 SFR >= {min_sfr} Msun/yr"
    )
    return [(int(i), int(n_gas[i]), float(sfr[i])) for i in inside]


def build_cube(basePath, snapNum, subID, hi_model, observed_redshift, outdir):
    """构建并写出一次模拟观测。

    返回值
    -------
    dict or None
        一条汇总记录；若该子晕没有可供观测的 HI，则为 None。
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

    # 让数据立方体的尺寸跟随源：固定的角尺寸要么把近处的星系裁掉，要么在远处
    # 的星系上浪费像素。速度轴要覆盖气体运动学，再加谱线轮廓的余量，因为谱线落到
    # 通道范围之外的粒子会被静默裁剪掉（见用户指南）。
    arcsec_per_kpc = (U.rad / source.D_A.to(U.kpc)).to(U.arcsec / U.kpc)
    extent_arcsec = 2.5 * 20.0 * arcsec_per_kpc.value  # 约合一个 20 kpc 星系的 2.5 倍
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
        # 每个立方体都用一份新的谱线模型：一个实例会缓存自己的谱数组，
        # 无法跨通道数不同的立方体复用。
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
        description="为一组 TNG 子晕构建模拟 HI 立方体。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--basePath",
        default=os.environ.get("TNG_BASE_PATH"),
        help="TNG 模拟目录。未给出时回退到环境变量 $TNG_BASE_PATH。",
    )
    parser.add_argument("--snapNum", type=int, default=99, help="快照号。")
    parser.add_argument("--limit", type=int, default=10, help="处理多少个子晕。")
    parser.add_argument("--minGas", type=int, default=64, help="每个子晕的最少气体单元数。")
    parser.add_argument(
        "--maxGas",
        type=int,
        default=5000,
        help="每个子晕的最多气体单元数。把又慢又吃内存的团核心排除在样本之外。",
    )
    parser.add_argument(
        "--hiModel",
        default="K13_vol",
        choices=["K13_vol", "K13_map", "GD14_vol", "GD14_map"],
        help="HI/H2 模型。",
    )
    parser.add_argument("--observedRedshift", type=float, default=0.02, help="观测红移。")
    parser.add_argument(
        "--minSfr",
        type=float,
        default=0.0,
        help="最低恒星形成率 [Msun/yr]。把它调到零以上，样本就只剩正在形成恒星的"
        "目标，自适应 FFT UV 场在这些目标上才真正起作用。",
    )
    parser.add_argument("--outdir", default=None, help="输出目录。")
    args = parser.parse_args()

    if not args.basePath:
        parser.error(
            "没有 TNG 数据目录。请传入 --basePath，或设置 TNG_BASE_PATH：\n"
            "    export TNG_BASE_PATH=/path/to/TNG/TNG-100/"
        )
    basePath = args.basePath.rstrip("/") + "/"
    outdir = args.outdir or default_outdir()
    os.makedirs(outdir, exist_ok=True)

    subhalos = select_subhalos(
        basePath, args.snapNum, args.minGas, args.maxGas, args.limit, args.minSfr
    )
    if not subhalos:
        print("没有选中任何目标；请放宽 --minGas/--maxGas。")
        return 1

    print(f"\n模型 {args.hiModel}，观测红移 z = {args.observedRedshift}")
    print(f"输出到 {outdir}\n")

    results, skipped, failures = [], [], []
    started = time.time()
    for index, (subID, n_gas, sfr) in enumerate(subhalos, start=1):
        print(
            f"[{index}/{len(subhalos)}] subID {subID} "
            f"（{n_gas} 个气体单元，SFR {sfr:.2f} Msun/yr）...",
            flush=True,
        )
        try:
            record = build_cube(
                basePath, args.snapNum, subID, args.hiModel, args.observedRedshift, outdir
            )
        except Exception as exc:  # 继续跑下去：一个坏目标不该中断整批任务
            failures.append((subID, repr(exc)))
            print(f"    失败：{exc}")
            continue

        if record is None:
            skipped.append(subID)
            print("    该模型下没有 HI；已跳过")
            continue

        results.append(record)
        print(f"    M_HI = {record['m_hi']:.3e} Msun，用到 {record['npart']} 个粒子")

    elapsed = time.time() - started
    print("\n" + "=" * 68)
    print(f"  写出 {len(results)} 个立方体，跳过 {len(skipped)} 个，"
          f"失败 {len(failures)} 个  ({elapsed:.0f} s)")
    print("=" * 68)

    if results:
        m_hi = np.array([r["m_hi"] for r in results])
        m_h2 = np.array([r["m_h2"] for r in results])
        print(f"  M_HI：中位数 {np.median(m_hi):.3e}，最大 {m_hi.max():.3e} Msun")
        print(f"  M_H2：中位数 {np.median(m_h2):.3e}，最大 {m_h2.max():.3e} Msun")
        print(f"\n  这些是 {args.hiModel} 模型给出的质量。四种模型彼此不同是刻意的：")
        print("  同一个星系上，H2 质量可以相差数倍 —— 所以任何要发表的数值")
        print("  都必须注明用的是哪个模型。")
    for subID, error in failures:
        print(f"  失败：subID {subID}，{error}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
