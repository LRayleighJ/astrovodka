# astrovodka 文档

**对 IllustrisTNG 星系的模拟分辨 HI（21 cm）观测。**

astrovodka 把 SPH 粒子数据变成合成的射电数据立方体。它是
[MARTINI](https://github.com/kyleaoman/martini) 的一个分支，把距离计算换
成了真正的宇宙学距离，把中性气体的划分换成了有物理依据的 HI/H₂ 模型。

英文文档见 [`docs/en/`](../en/index.md)。

---

## 从哪里开始

| 如果你想…… | 请读 |
|---|---|
| **把它跑起来** | [安装与环境配置](installation.md) — 三条安装路径、可选依赖、离线集群、TNG 数据布局、疑难排解 |
| **看它实际运行** | [`examples/zh/`](../../examples/zh/) — 一个基于真实 TNG-100 子晕的教程笔记本，外加可运行的脚本 |
| **理解 API** | [用户指南](guide.md) — 流水线、源、距离、数据立方体、常见任务，以及四个耗掉真实调试时间的陷阱 |
| **理解物理** | [HI/H₂ 模型](hi-models.md) — 四种模型各自计算什么、彼此如何不同、它们的假设又会在哪里带来麻烦 |
| **查阅签名** | [API 参考](api/index.md) |
| **核对论断** | [验证记录](verification.md) — 与参考实现逐位（bit-for-bit）一致，以及端到端检查 |

---

## 让它有别于 MARTINI 的两件事

**宇宙学距离。** `D_L` 归一化流量，`D_A` 决定每一个角度尺度，二者都来自
`Planck15`。它们相差 `(1+z)²`，用错一个是在这里能犯的最昂贵的错误。见
[两种距离](guide.md#两种距离)。

**四种可选的 HI/H₂ 模型。** Diemer+2018 把 Gnedin & Draine 2014 与
Krumholz+ 2013 两两配对，每一种都分体密度变体与投影变体，用一个关键字即
可选定。它们彼此不一致是刻意为之 — 在某个子晕上，`K13_vol` 与 `GD14_vol`
的 H₂ 质量相差约 6 倍 — 所以你用的是哪一种，本身就是结果的一部分。见
[HI/H₂ 模型](hi-models.md)。

---

## 文档地图

```
docs/
├── en/                      本目录树
│   ├── index.md             你在这里
│   ├── installation.md      环境、依赖、TNG 数据
│   ├── guide.md             概念与完整示例
│   ├── hi-models.md         四种模型的深入讲解
│   ├── verification.md      物理是如何验证的
│   └── api/                 逐模块参考
├── zh/                      同一目录树，中文版
└── html/                    由上述 Markdown 生成
```

Markdown 才是事实来源。`docs/html/` 由 `python docs/build_html.py`
生成，一旦两者出现偏差，测试就会失败 — 此前手工维护的 HTML 与 Markdown
两棵树就出现过这种情况，最终在版本号和某个参数的默认值上都不一致。
