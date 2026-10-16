# astrovodka API 参考

Version 1.0.0 — 对 IllustrisTNG 星系生成模拟 HI 观测，具备宇宙学距离，
并提供四种可选的 Diemer+2018 HI/H2 模型。

这里是逐模块的接口参考。背后的概念请先看[用户指南](../guide.md)；
HI/H2 物理见 [HI/H2 模型](../hi-models.md)；安装见
[安装与环境配置](../installation.md)。

## 模块

| 模块 | 说明 |
|--------|-------------|
| [martini](martini.md) | 核心流水线 — Martini 与 GlobalProfile 两个类 |
| [datacube](datacube.md) | DataCube — WCS 网格定义 |
| [sph_source](sph_source.md) | SPHSource — 通用 SPH 源（红移模式 + 距离模式） |
| [tng_source](tng_source.md) | TNGSource — IllustrisTNG 源 + 可选 HI 模型 |
| [beams](beams.md) | GaussianBeam — 望远镜波束模型 |
| [noise](noise.md) | GaussianNoise — 噪声模型 |
| [spectral_models](spectral_models.md) | GaussianSpectrum、DiracDeltaSpectrum — 21-cm 谱线模型 |
| [sph_kernels](sph_kernels.md) | SPH 平滑核（6 种，带自适应回退） |
| [cosmology_utils](cosmology_utils.md) | 宇宙学距离与单位换算辅助函数 |
| [himodel](himodel.md) | 四种可选 HI/H2 模型（GD14/K13，体密度/投影） |
| [data_products](data_products.md) | SubhaloDataProducts — 三维网格化 UV/DM/HI/H2 |
| [group_products](group_products.md) | GroupDataProducts — FOF group 的气体/暗物质，含 ISM/CGM/IGM 分类 |
| [L_coords](L_coords.md) | 角动量对齐坐标 |
| [_illustris_tools](_illustris_tools.md) | TNG 快照数据 I/O 工具 |
