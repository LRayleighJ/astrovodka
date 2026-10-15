#!/usr/bin/env python
"""astrovodka 演示脚本。

用 astrovodka 生成模拟的、空间可分辨的 HI（21 cm）数据立方体：数据要么来自
一个简单的解析星系，要么来自本地磁盘上的真实 IllustrisTNG 子晕。不带参数
运行时跑的是解析星系，不需要任何模拟数据；传入 ``--basePath``（或设置
``TNG_BASE_PATH``）则改走 TNG 这条路。

用法
----
    # z = 0.01 的解析旋转盘，不需要模拟数据
    python examples/zh/astrovodka_demo.py

    # 真实的 TNG-100 子晕（路径也可来自环境变量 $TNG_BASE_PATH）
    python examples/zh/astrovodka_demo.py \\
        --basePath /path/to/TNG/TNG-100/ --snapNum 99 --subID 113 \\
        --hiModel K13_vol

输出默认写到 ``examples/outputs/``（用 ``--outdir`` 可改），该目录在 .gitignore 中。
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
    """返回 ``examples/outputs/``，相对于本文件解析。

    以 ``__file__`` 而不是当前工作目录为基准来解析，意味着脚本在仓库里的
    任何位置运行，写出的位置都相同。
    """
    here = Path(__file__).resolve()
    outdir = here.parent.parent / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    return str(outdir)


def demo_simple_source(output_file):
    """在 z = 0.01 的解析旋转盘上跑完整条流水线。

    这是整个模拟观测流程的微缩版，不需要任何模拟数据：构造源、定义数据立方体、
    选好各个子模块、运行。

    参数
    ----
    output_file : str
        要写出的 FITS 立方体的路径。
    """
    print("=" * 68)
    print("  astrovodka 演示 1：z = 0.01 的解析旋转盘")
    print("=" * 68)

    # ---- 一个解析盘星系：不需要模拟数据 ----
    n_part = 1000
    rng = np.random.default_rng(42)
    phi = rng.random(n_part) * 2 * np.pi
    radius = rng.random(n_part) * 10 * U.kpc

    # 垂直方向为指数分布，关于中面对称。
    z_height = -np.log(rng.random(n_part)) * 0.5 * U.kpc
    z_height *= np.sign(rng.random(n_part) - 0.5)

    xyz_g = U.Quantity(
        [
            np.zeros(n_part) * U.kpc,  # x 是视线方向
            radius * np.cos(phi),  # y-z 是天空平面
            radius * np.sin(phi) + z_height,
        ]
    )

    # 100 km/s 的平坦旋转曲线。
    v_phi = 100 * U.km / U.s
    vxyz_g = U.Quantity(
        [
            rng.standard_normal(n_part) * 5 * U.km / U.s,
            -v_phi * np.sin(phi),
            v_phi * np.cos(phi),
        ]
    )

    T_g = np.full(n_part, 1e4) * U.K
    mHI_g = np.full(n_part, 1e8) * U.Msun  # 每个粒子 10^8 Msun，总计 10^11 Msun
    hsm_g = np.full(n_part, 1.0) * U.kpc

    # ---- 源：放在一个宇宙学红移上 ----
    # 传 `redshift`（而不是 `distance`）会让 SPHSource 进入宇宙学模式：
    # D_L 由 Planck15 给出，负责流量归一化；D_A 决定所有角尺度。
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

    print("\n源")
    print(f"  红移                     : {redshift}")
    print(f"  光度距离 D_L             : {source.distance:.1f}")
    print(f"  系统速度                 : {source.vsys.to(U.km / U.s):.0f}")
    print(f"  粒子数                   : {source.npart}")
    print(f"  HI 总质量                : {source.input_mass:.2e}")

    # ---- 数据立方体 ----
    # spectral_unit 保持默认时，谱轴会采用 channel_width 的单位（这里是速度轴），
    # 与 MARTINI 的行为一致。
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
        f"\nDataCube：{datacube.n_px_x} x {datacube.n_px_y} 个像素，"
        f"{datacube.n_channels} 个通道"
    )

    # ---- 仪器与谱线模型 ----
    beam = GaussianBeam(bmaj=8 * U.arcsec, bmin=6 * U.arcsec, bpa=30 * U.deg, truncate=3.0)
    noise = GaussianNoise(rms=1e-4 * U.Jy * U.beam**-1, seed=0)
    spectral_model = GaussianSpectrum(sigma=7 * U.km / U.s)
    sph_kernel = CubicSplineKernel()

    # ---- 运行流水线 ----
    print("\n正在运行流水线...")
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

    print(f"\n已写出 {output_file}")
    print("完成。")


def demo_tng_source(basePath, snapNum, subID, observed_redshift, hi_model, output_file):
    """在从本地磁盘读取的真实 IllustrisTNG 子晕上跑完整条流水线。

    参数
    ----
    basePath : str
        TNG 模拟输出目录，其中包含 ``groups_<snap>/`` 和 ``snapdir_<snap>/``。
    snapNum : int
        快照号。
    subID : int
        该快照中的子晕 ID。
    observed_redshift : float
        放置模拟观测所用的红移。
    hi_model : str
        ``K13_vol``、``K13_map``、``GD14_vol``、``GD14_map`` 之一。
    output_file : str
        要写出的 FITS 立方体的路径。
    """
    print("=" * 68)
    print(f"  astrovodka 演示 2：TNG-100 snap={snapNum} subID={subID}")
    print("=" * 68)

    from astrovodka.sources import TNGSource

    # `observed_redshift` 把这里真正起作用的两个红移分开。快照自身的红移决定物理
    # （单位换算、UV 背景）；`observed_redshift` 决定几何（D_L、D_A，进而决定角尺度
    # 与流量）。对 snap 99 这种 z~0 的快照，角直径距离几乎为零，没有它星系
    # 就只是一个无法分辨的点。
    #
    # `hi_model` 从四种 HI/H2 模型中选一种。投影版（*_map）还需要恒星粒子，
    # 并且把落在投影盒外的气体记为 HI 为零。
    source = TNGSource(
        basePath=basePath,
        snapNum=snapNum,
        subID=subID,
        hi_model=hi_model,
        observed_redshift=observed_redshift,
        ra=0.0 * U.deg,
        dec=0.0 * U.deg,
    )

    print("\n源")
    print(f"  HI 模型                  : {source.hi_model}")
    print(f"  快照红移                 : {source.source_redshift}")
    print(f"  观测红移                 : {observed_redshift}")
    print(f"  宇宙学                   : Planck15")
    print(f"  光度距离 D_L             : {source.D_L:.1f}")
    print(f"  角直径距离 D_A           : {source.D_A:.1f}")
    print(f"  系统速度                 : {source.vsys.to(U.km / U.s):.0f}")
    print(f"  粒子数                   : {source.npart}")
    print(f"  HI 总质量                : {source.m_HI_total:.4e}")
    print(f"  H2 总质量                : {source.m_H2_g.sum():.4e}")

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
    spectral_model = GaussianSpectrum(sigma="thermal")  # 逐粒子 kB*T/m_p
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

    print(f"\n已写出 {output_file}")
    print("完成。")


def main():
    parser = argparse.ArgumentParser(
        description="astrovodka 演示 —— 合成 HI 观测",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--basePath",
        type=str,
        default=os.environ.get("TNG_BASE_PATH"),
        help="TNG 模拟输出目录。未给出时回退到环境变量 $TNG_BASE_PATH。"
        "两者都没有设置时，改为运行解析星系演示。",
    )
    parser.add_argument("--snapNum", type=int, default=99, help="TNG 快照号。")
    parser.add_argument(
        "--subID",
        type=int,
        default=113,
        help="子晕 ID。默认值 113 是 TNG-100 快照 99 中一个气体丰富、正在形成恒星"
        "的子晕（M_HI ~ 4e9 Msun）。避免使用太小的 ID：群目录按质量降序排列，"
        "所以 ID 0 是团核心（约 3800 万个气体单元，读入需要超过 10 GB）。",
    )
    parser.add_argument(
        "--observedRedshift",
        type=float,
        default=0.02,
        help="放置模拟观测所用的红移，与快照自身的红移解耦。",
    )
    parser.add_argument(
        "--hiModel",
        type=str,
        default="K13_vol",
        choices=["K13_vol", "K13_map", "GD14_vol", "GD14_map"],
        help="使用哪一种 HI/H2 模型。",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help="输出 FITS 文件的目录。默认为 examples/outputs/。",
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
