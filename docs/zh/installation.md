# 安装与环境配置

本页涵盖让 astrovodka 运行起来所需的全部内容：支持的 Python 版本、三条受
支持的安装路径、可选依赖以及各自解锁的功能、如何在集群上离线工作、如何组
织你的 IllustrisTNG 数据，以及如何确认安装是否完好。

如果你只想要最短的路径，那就是：

```bash
git clone https://github.com/LRayleighJ/astrovodka.git
cd astrovodka
conda env create -f environment.yml
conda activate astrovodka
```

然后跳到[验证安装](#验证安装)。

- [环境要求](#环境要求)
- [安装](#安装)
  - [方式 1：conda（推荐）](#方式-1conda推荐)
  - [方式 2：virtualenv 与 pip](#方式-2virtualenv-与-pip)
  - [方式 3：仅核心，不装可选依赖](#方式-3仅核心不装可选依赖)
- [可选依赖](#可选依赖)
- [离线与隔离网络安装](#离线与隔离网络安装)
- [IllustrisTNG 数据](#illustristng-数据)
- [验证安装](#验证安装)
- [疑难排解](#疑难排解)
- [Python 版本支持](#python-版本支持)

---

## 环境要求

### 软件

| 要求 | 最低版本 | 说明 |
|---|---|---|
| Python | 3.9 | 见 [Python 版本支持](#python-版本支持) |
| pip | 21.3 | 只有 `pyproject.toml` 的项目需要它 |
| numpy | 1.20 | 核心依赖 |
| scipy | 1.6 | FFT 卷积、插值、特殊函数 |
| astropy | 5.0 | 宇宙学、单位、坐标、WCS、FITS 读写 |
| tqdm | 4.40 | 进度条 |
| h5py | 3.0 | *可选，但处理 TNG 数据时必需* — 见下文 |
| matplotlib | 3.4 | 可选：绘图与教程 |
| multiprocess | 0.70 | 可选：`ncpu > 1` 时的源插入 |

`numpy`、`scipy`、`astropy` 和 `tqdm` 会随本包自动安装。其余都是按需选择
的可选依赖组（extra），只不过**任何 TNG 工作流都需要 `h5py`** — 它负责读
取快照文件，而不只是写出 HDF5 结果。

注意 astrovodka 使用的是 **`multiprocess`**，而不是标准库的
`multiprocessing`。MARTINI 用基于 `dill` 的序列化来做并行，标准库模块并
不提供这一能力。只有当你传入 `ncpu > 1` 时才需要安装 `multiprocess`。

### 硬件

本包本身很轻量，但有两个功能很吃内存：

- **加载一个 TNG 子晕。** 开销随该子晕的粒子数增长。教程的目标 — TNG-100
  快照 99 的子晕 113 — 有 8,169 个气体单元，几秒内就能加载完。星系团核心
  则是另一回事：同一快照的子晕 0 约有 3,800 万个气体单元，需要远超 10 GB
  内存。请避开低子晕 ID，那些是最庞大的天体（目录按质量降序排序）。
- **三维数据产品与自适应 FFT UV 场** 会构建完整的 `ngrid³` 网格（默认
  128³），每个数组的每个网格约 10 MB。气体单元多于约 10⁷ 的子晕本可从分
  块处理中受益，但分块尚未实现。

软件本身的磁盘占用可以忽略不计（几 MB）。TNG 数据的磁盘占用则不然 — 见
[IllustrisTNG 数据](#illustristng-数据)。

---

## 安装

astrovodka 不在 PyPI 上。请从源码检出安装。下面三条路径都从这一组命令开
始：

```bash
git clone https://github.com/LRayleighJ/astrovodka.git
cd astrovodka
```

### 方式 1：conda（推荐）

仓库根目录的 `environment.yml` 固定了项目开发与验证所依据的整套依赖，并
以可编辑（editable）模式安装 astrovodka 本身：

```bash
conda env create -f environment.yml
conda activate astrovodka
```

若要使用别的环境名，请先修改 `name:` 一行，或在命令行上覆盖：

```bash
conda env create -f environment.yml -n my-env
```

拉取新提交后，更新已有环境：

```bash
conda activate astrovodka
conda env update -f environment.yml --prune
```

### 方式 2：virtualenv 与 pip

任何 Python 3.9+ 的环境都可以：

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

`[dev]` 会带入 h5py、matplotlib、multiprocess、pytest 和 pytest-cov。若
想要更精简的安装，只点名你需要的可选依赖组：

```bash
python -m pip install -e ".[hdf5]"                 # 用于 TNG 工作，不装其他
python -m pip install -e ".[hdf5,parallel]"        # 另加多核插入
python -m pip install -e ".[hdf5,parallel,plotting]"   # 另加绘图
```

如果你更想用普通的 requirements 文件而不是可选依赖组：

```bash
python -m pip install -r requirements-dev.txt
python -m pip install -e . --no-deps
```

代码还在变动时推荐用可编辑（`-e`）方式。若要固定安装，请改用
`pip install .`，但注意示例和教程都假定仓库仍然存在。

### 方式 3：仅核心，不装可选依赖

物理层 — `astrovodka.himodel`，从 TNG 气体单元数组计算逐单元的 HI 与 H2
质量 — 只依赖 numpy、scipy 和 astropy。它可以单独使用，没有 h5py，没有
matplotlib，也完全见不到数据立方体：

```bash
python -m pip install -e .
```

```python
from astrovodka.himodel import HIModelConfig, compute_hi_mass

cfg = HIModelConfig(h=0.6774, z=0.0, boxsize_ckpc_per_h=75000.0, ngrid=64)
result = compute_hi_mass(gas, model="K13_vol", cfg=cfg)
print(result.m_hi_msun, result.m_h2_msun)
```

包的布局正是为了让这一点成立：`astrovodka.himodel` 不从观测层导入任何东
西，`import astrovodka` 也不会导入 h5py 或 matplotlib。有一项测试同时保
证这两点（`tests/test_selfcontained.py`）。

---

## 可选依赖

| 可选依赖组 | 包 | 解锁的功能 | 适用场景 |
|---|---|---|---|
| `hdf5` | h5py | 读取 TNG 快照与群目录；`write_hdf5`、`save_state`、`SubhaloDataProducts.save_hdf5` | **任何 TNG 工作流** |
| `parallel` | multiprocess | `insert_source_in_cube(ncpu=N)`、多核光谱初始化 | 大型立方体 |
| `plotting` | matplotlib | `Martini.preview()`、`GlobalProfile.plot_spectrum()`、教程插图 | 出图 |
| `testing` | pytest, pytest-cov | 测试套件 | 开发 |
| `all` | 全部三个运行时可选依赖组 | — | 图方便 |
| `dev` | 以上全部 + testing | — | 开发 |

`h5py` 之所以列为可选，是因为物理层确实不需要它。但在实践中，凡是要接触
真实 TNG 数据的路径都需要它，所以除非你只用 `astrovodka.himodel`，否则请
把它当作必需项。

---

## 离线与隔离网络安装

集群登录节点往往没有对外网络。有两件事需要注意。

**1. 构建隔离。** 默认情况下 pip 会创建一个临时虚拟环境来构建本包，并把
`setuptools` 下载进去。如果你已经有较新的 setuptools，就跳过这一步：

```bash
python -m pip install -e . --no-build-isolation
```

`pyproject.toml` 要求 `setuptools>=64`；用
`python -c "import setuptools; print(setuptools.__version__)"` 查看你的
版本。

**2. 依赖。** 要么在联网机器上预先下载 wheel 包，

```bash
# 在联网机器上执行，Python 版本需与目标一致
pip download astrovodka[hdf5,parallel,plotting] -d wheels/ --no-deps
pip download numpy scipy astropy tqdm h5py matplotlib multiprocess pytest -d wheels/
```

然后把 `wheels/` 复制过去并从它安装：

```bash
python -m pip install --no-index --find-links wheels/ -e . --no-build-isolation
```

要么在联网机器上生成一个 conda 环境文件再传过去
（`conda env export > environment.lock.yml`，然后
`conda env create -f environment.lock.yml`）。

**3. conda 离线。** 如果包已经在本地的包缓存里，
`conda env create --offline` 就能用。如果你的集群提供了 TNG 数据却没有科
学 Python 栈，那么先 `conda create -n astrovodka --clone base`、再
`pip install -e . --no-deps` 往往是最快的入门方式，前提是 numpy、scipy
和 astropy 已经就位。

---

## IllustrisTNG 数据

astrovodka **从本地磁盘**读取 TNG 数据。它不使用 TNG 的网页 API，运行时
也不需要网络访问，但它需要那些文件。

### 从哪里获取

- **TNG-100 / TNG-300 公开数据发布：**
  <https://www.tng-project.org/data/>（需要注册）。下载你所需那些快照的
  `groups_<snap>` 和 `snapdir_<snap>` 目录。
- **TNG-50：** 同一站点，体积更小，适合测试。
- 如果你只想要群目录（子晕的位置与质量）而不想要粒子数据，那么仅有
  `groups_<snap>` 就足以*挑选*子晕，但不足以构建立方体。

### 预期的布局

把 `basePath` 指向模拟目录，也就是*包含* `groups_<snap>` 和
`snapdir_<snap>` 的那一层：

```
/path/to/TNG/TNG-100/
├── groups_099/
│   ├── fof_subhalo_tab_099.0.hdf5      # 群/子晕目录
│   └── ...
└── snapdir_099/
    ├── snap_099.0.hdf5                 # 粒子数据，分块存储
    └── ...
```

快照号会补零到三位数字。目录名必须符合这一模式；`TNGSource` 会从
`basePath` 和 `snapNum` 构造路径。

### 大小

对 TNG-100 而言，一个快照的粒子数据约为 **1–2 TB**，群目录约 **13 GB**。
你并不需要每个快照 — 对多数工作来说一两个就够。如果数据放在带清理策略的
暂存文件系统上，可以考虑搭一组符号链接（symlink farm）。

### 告诉 astrovodka 数据在哪里

没有全局配置文件。请显式传入路径：

```python
source = TNGSource(basePath="/path/to/TNG/TNG-100/", snapNum=99, subID=113)
```

示例脚本和两个教程笔记本则改从环境变量中读取它，这样仓库里不会写死任何路
径：

```bash
export TNG_BASE_PATH=/path/to/TNG/TNG-100/
python examples/en/astrovodka_demo.py --snapNum 99 --subID 113
```

---

## 验证安装

### 1. 能导入

```bash
python -c "import astrovodka; print(astrovodka.__version__)"
python -c "from astrovodka.himodel import METHODS; print(METHODS)"
```

第二条应当输出 `('GD14_vol', 'GD14_map', 'K13_vol', 'K13_map')`。

### 2. 可选组件都在

```bash
python -c "import h5py, matplotlib, multiprocess; print('optional deps OK')"
```

缺少 `h5py` 不会让 `import astrovodka` 失败 — 这是有意为之 — 但你一旦加
载 TNG 子晕，它就会立刻报错。

### 3. 测试套件通过

```bash
python -m pytest tests/ -q
```

该套件是合成的：**不需要 TNG 数据，不需要网络，也不需要 GPU**，约 15 秒
就能跑完。它对照解析预期检验四种 HI 模型，还覆盖宇宙学辅助函数、数据立方
体的谱轴换算、三维数据产品，以及上文描述的自包含规则。

### 4. 示例可以运行

不需要任何 TNG 数据时：

```bash
python examples/en/astrovodka_demo.py
```

这会构建一个 z = 0.01 的解析旋转盘，并把 FITS 立方体写到
`examples/outputs/`。它会打印一项质量守恒检查；下面这一行

```
Mass in cube (assuming distance ... ): 1.00e+11 solMass
    [100% of initial source mass]
```

应当显示 100%。

有 TNG 数据时：

```bash
export TNG_BASE_PATH=/path/to/TNG/TNG-100/
python examples/en/astrovodka_demo.py --snapNum 99 --subID 113 --hiModel K13_vol
```

---

## 疑难排解

| 症状 | 原因 | 修复 |
|---|---|---|
| `ModuleNotFoundError: No module named 'astrovodka'` | 没有安装，或装进了另一个解释器 | 用与运行脚本*同一个* `python` 执行 `python -m pip install -e .` |
| 构造 `TNGSource` 时出现 `ModuleNotFoundError: No module named 'h5py'` | h5py 属于可选依赖组 | `pip install -e ".[hdf5]"` |
| `ModuleNotFoundError: No module named 'six'` | 较早版本的随附 TNG I/O 使用了 `six` | 已在 1.0.0 中更新；请更新你的检出 |
| 在 Python 3.9/3.10 上出现 `ModuleNotFoundError: No module named 'typing_extensions'` | `typing.Self` 到 3.11 才进入标准库 | `pip install typing_extensions`，或重新安装本包，让条件依赖被解析出来 |
| 使用 `ncpu > 1` 时出现 `ModuleNotFoundError: No module named 'multiprocess'` | 标准库的 `multiprocessing` 不能替代它 | `pip install -e ".[parallel]"` |
| `ImportError: numpy.core.multiarray failed to import` | numpy/scipy 是针对不匹配的 ABI 构建的 | 在全新环境里执行 `pip install --force-reinstall numpy scipy` |
| `ValueError: datacube padding insufficient for beam convolution` | `write_fits` 会调用 `drop_pad()`，因此在流水线中途写出会移除填充区 | 在下一次 `convolve_beam()` 之前执行 `M.datacube.add_pad(M.beam.needs_pad())` |
| `FileNotFoundError: .../snapdir_099/...` | `basePath` 指向了 `snapdir_099/` 本身，或指向 `TNG-100/` 的上一级 | `basePath` 必须是*包含* `groups_<snap>` 和 `snapdir_<snap>` 的那个目录 |
| 加载一个子晕要花几分钟、几百 GB | 你正在加载星系团核心 | 使用 100 及以上的子晕 ID；见 [硬件](#硬件) |
| `RuntimeError: unknown hi_model` | 拼写错误，或用了旧版本中的名称 | 取 `K13_vol`、`K13_map`、`GD14_vol`、`GD14_map` 之一 |

---

## Python 版本支持

**Python 3.9 是最低版本，也是本次发布所测试的版本**，CI 中连同 3.10 到
3.13 一起覆盖。`pyproject.toml` 声明了 `requires-python = ">=3.9"`。

有两个细节值得了解：

- `astrovodka/datacube.py` 使用了 `typing.Self`，它到 Python 3.11 才进入
  标准库。在 3.9 和 3.10 上它会回退到 `typing_extensions`，后者被声明为
  条件依赖，因此在这些版本上会自动安装。
- `types.EllipsisType`（3.10+）采用防御式导入，并带有 `type(Ellipsis)`
  回退，因此该标注在 3.9 上也能正常工作。

由于 astropy 7 要求 Python 3.11 及更高版本，在 Python 3.9 上 pip 会解析
到较旧的 astropy（6.x），而它是受支持的。astropy 声明的最低版本是 5.0。

随仓库附带的 MARTINI 观测机制在这里没有改写，但它是干净的 3.9 代码，CI
会在每个受支持的版本上运行整套测试。
