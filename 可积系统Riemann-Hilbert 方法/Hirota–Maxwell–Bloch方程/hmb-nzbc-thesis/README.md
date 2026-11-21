# HMB 非零边界论文

本项目以完成毕业论文为目标，研究非零光场背景下的 HMB 方程。当前正文为八节，暂不包含摘要和附录；组会用于汇报进度。

1. 绪论
2. HMB 模型、基本条件与背景谱问题
3. 非零背景的直接散射
4. 散射演化与共振相容性
5. 矩阵 RHP 与三场重构
6. 无反射 RHP 与一般重构
7. 孤子解的显式构造
8. 数值结果与动力学分析

入口为 fduthesis.tex，正文位于 section，成品为 fduthesis.pdf。相邻 hmb-nzbc-thesis-explanation 保留逐式推导，与当前正文的公式和定理编号同步。

2026-11-15 的论文式修订将小节合并为12个。公式整体居中，式2.5、2.21等多行条件在公式内部左端对齐；数学对象采用定义环境，辅助讨论使用“注”。命题和定理不另加名称，证明采用编号及普通字体的简短步骤说明。规范变换列为命题2.2，规范后三场基本条件列为注2.1；式2.17列出Lax方程，式2.18、2.19分别定义两个Lax矩阵。背景矩阵与传播函数集中于2.3节，正文省去其分量计算。

2.3节从背景特征方程引出谱曲面、单值化及谱边界；特征向量矩阵、渐近基本解和Jost积分方程集中于3.1节，不再单列Jost定义环境。介质背景矩阵在第4节共振相容性中首次引入。对应推导与文献核对见 calculations/background-jost-rewrite.md，精确矩阵核验脚本为 scripts/verify_background_jost.py。

## 数学范围

正向散射以统一场类为基础，时间演化另加 HMB 方程条件。主要散射对象为空间 Jost 解和演化散射矩阵，完整相位用于辅助分析与逆验证。

一般含反射表示保留解存在、局部正则性与可积性条件；一般 $N$ 轨道无反射情形通过 Cauchy 矩阵与 Hermitian 分块合同变换证明全时空 $\det L_N(x,t)>0$，并统一完成空间端态、完整 Lax 对恢复、原 HMB 方程、Bloch 约束与直接散射数据复原。第 7 节在此基础上显式构造单孤子与一般复谱双孤子解并导出传播与弹性碰撞位移公式；第 8 节给出单、双孤子的数值图像、动力学分析与 70 位高精度残差核验。

## 编译与检查

在本目录运行 PowerShell：

    .\build.ps1

脚本在 ASCII 临时目录完成 XeLaTeX、BibTeX、XeLaTeX、XeLaTeX，再复制 PDF 与 SyncTeX 回项目。可用 OutputDirectory 参数另选目录。

scripts 保留符号计算、留数系统、单轨道三场、一般多轨道样点及 Cauchy 表示的核验。有限样点不替代一般数据的解析证明。

calculations/literature-rewrite-review.md 记录跨来源文献与读取范围；calculations/structural-rewrite-review.md 记录重写结构及证明范围；PDF 验收见 calculations/pdf-acceptance.json。
