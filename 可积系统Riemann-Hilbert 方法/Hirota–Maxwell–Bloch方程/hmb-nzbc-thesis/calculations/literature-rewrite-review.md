# HMB 结构重写的文献与证明范围核查

核查日期：2026-11-14。基准提交：8690a1b。最终目标为毕业论文，组会只承担进度汇报。本轮不增加摘要和附录。

## 检索方法

通过期刊官网、arXiv、作者公开全文、MathNet 与 SciSpace 交叉检索，不限本地论文。检索包括 HMB / NLS–MB 的 RHP、sharp-line MB 的固定谱极点、非零光场散射、谱奇点、数据传播、介质恢复、多孤子与适定性，并追踪核心论文的引用和后续工作。

索引中的作者、年份及重复版本以原刊为准。GPT 对话只提供线索。下表“关键段落”不代表逐页认证全篇；只核对摘要的论文用于背景介绍，不作为未核查定理的证明来源。

## 文献对照

| 文献 | 核查材料 | 与本文的关系 |
| --- | --- | --- |
| Li–He–Porsezian, PRE 87 (2013), 012913 | 本地全文与 arXiv v3；Lax、Darboux、三场背景解 | 采用已有 Lax 表示并重新验证零曲率；已有呼吸子和怪波，本稿不把该类波形本身称为新解 |
| Li–Xia–Wei, IJMPB 35 (2021), 2150153 | 本地原文的模型、无反射系统及重构 | 零光场端态；有限维计算可比较，不能替代非零光场谱曲面 |
| Li–Tian–Peng–Yang, TMP 203 (2020), 709–725 | 原刊记录与摘要 | 高阶 NLS–MB 的 RHP 与孤子分类；边界和参数需分别核对 |
| Liu–Zhang–Xu, MPLB 35 (2021), 2150420 | 原刊记录与摘要 | NLS–MB 有理演化结构的线索；不从标题认定光场端态与本稿相同 |
| Guo–Zhang–Liu, JMP 63 (2022), 051505 | AIP 原刊作者、摘要、卷页与 DOI | 高阶 NLS–MB 的无反射公式和数值图；不从摘要推断全局非退化性 |
| Liu–Liu, Nonlinear Dynamics 114 (2026), 821 | Springer 摘要与记录 | 摘要明确仅介质变量非零，光场与极化趋于零；不同于本文光场背景 |
| Biondini–Gabitov–Kovačič–Li, JMP 60 (2019), 073510 | 作者全文：传播、逆问题、三场恢复及无反射系统 | 核心非零光场参考；具有频率相关介质和展宽，传播公式不能直接照搬 |
| Kotlyarov–Moskovchenko, JMFAG 10 (2014), 328–349 | MathNet 原刊全文的 RHP 与重构段落 | 无展宽混合问题，强调两个 Lax 方程；端态和定义域不同 |
| Filipkovska–Kotlyarov–Melamedova, JMFAG 13 (2017), 119–153 | 原刊摘要与元数据；部分全文入口受限 | 规范变换、变换算子和奇点附近 Jost 控制的方法线索；未借用未核查的具体定理 |
| Li–Miller, CPAM 77 (2024), 457–542 | arXiv / 作者全文：假设、因果 RHP 与衰减讨论 | 奇异演化影响边界与衰减；不能自动宣称演化保持所选空间 |
| Wu–He–Huang, arXiv:2505.06986v2 (2025) | 原文第2节的 Lax、反射演化、修改曲线与介质恢复 | 实共振附近快速振荡；绕行曲线不能直接移植到本文实轴加圆周 |
| Biondini–Prinari–Zhang, ANONA 13 (2024), 20240054 | 作者全文：模型、初边值条件、空间及主结果 | 展宽 MB 的适定性；用于区分构造与适定性目标，不声称覆盖局部 HMB |
| Bilman–Miller, CPAM 72 (2019), 1722–1805 | 作者全文的谱奇点与 robust IST 部分 | 非零背景谱奇点参考；需区别 NLS 谱奇点与 HMB 固定时间极点 |
| Abeya–Biondini–Kovačič–Prinari, CMP 405 (2024), 192 | 作者全文的模型、边界与摘要 | 单侧非零背景的反射性质不同，不能直接使用其数据类 |
| Sitai Li, CMP 406 (2025), 89 | arXiv:2402.02166 当前版本与作者发表记录 | 零背景 sharp-line MB 的多孤子与局域性，作为进一步分析线索 |
| Li–Biondini–Kovačič, SAM 154 (2025), e70055 | 作者记录与 arXiv:2405.10117 摘要 | 非零背景耦合 MB 的传播和介质状态，模型与矩阵维数不同 |
| Zhang–Fan, ZAMP 71 (2020), 149 | 本地全文，重点第3、8、9节 | 参考空间分析、迹公式和有限维构造，不再逐章模仿 |

## 对正文的影响

1. 主定义采用空间 Jost 解及一个随时间演化的散射矩阵。完整相位作为辅助，在共振处不作点值代入。
2. 共振反射消失由介质方程与共同端态证明；点值结论不自动给出谱导数可积性。
3. 一般 RHP 的特殊点边界、Cauchy 密度可积性及迹公式对数条件单列。积分表示不替代解存在证明。
4. 一般多轨道公式保留矩阵可逆条件，指定纯虚单轨道则证明全局严格正行列式。
5. 保留两共振原像、Hermitian 约化、Bloch 约束及完整 Lax 对恢复的证明。
6. 单轨道矩阵有一次 Darboux 因子的代数结构；本稿定位为 RHP 内的自足构造与严格逆验证，不作“首次发现该类解”的声明。

## 原文入口

- [LHP2013](https://arxiv.org/abs/1205.1191)
- [LXW2021](https://doi.org/10.1142/S0217979221501538)
- [LTPY2020](https://www.mathnet.ru/eng/tmf9780)
- [LZX2021](https://doi.org/10.1142/S0217984921504200)
- [Guo2022](https://doi.org/10.1063/5.0076361)
- [LL2026](https://doi.org/10.1007/s11071-026-12724-0)
- [BGKL2019](https://www.math.buffalo.edu/~biondini/papers/jmp2019v60p073510.pdf)
- [KM2014](https://www.mathnet.ru/eng/jmag598)
- [FKM2017](https://www.mathnet.ru/eng/jmag666)
- [LM2024](https://arxiv.org/abs/2105.13293)
- [WHH2025](https://arxiv.org/abs/2505.06986)
- [BPZ2024](https://www.math.buffalo.edu/~biondini/papers/anona2024v13p20240054.pdf)
- [BM2019](https://arxiv.org/abs/1710.06568)
- [单侧非零背景 MB](https://www.math.buffalo.edu/~biondini/papers/cmp2024v405p192.pdf)
- [Li 的多孤子论文](https://arxiv.org/abs/2402.02166)
- [Li 的作者发表记录](https://www.sitaili.net/publications)
