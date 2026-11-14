# HMB 非零边界论文

本项目以完成毕业论文为目标，研究非零光场背景下的 HMB 方程。当前正文为八节，暂不包含摘要和附录；组会用于汇报进度。

1. 绪论
2. HMB 模型、基本条件与背景谱问题
3. 非零背景的直接散射
4. 散射演化与共振相容性
5. 矩阵 RHP 与三场重构
6. 无反射 RHP 与轨道构造
7. 单轨道三场解的全局性与方程验证
8. 数值结果与动力学分析

入口为 fduthesis.tex，正文位于 section，成品为 fduthesis.pdf。相邻 hmb-nzbc-thesis-explanation 保留逐式推导，与当前正文的公式和定理编号同步。

## 数学范围

正向散射以统一场类为基础，时间演化另加 HMB 方程条件。主要散射对象为空间 Jost 解和演化散射矩阵，完整相位用于辅助分析与逆验证。

一般含反射表示保留解存在、局部正则性与可积性条件；一般多轨道公式保留系数矩阵可逆条件。文中指定纯虚单轨道则证明全局矩阵可解、三场边界、基本空间归属与原方程。数值部分采用重构公式求值、参数分析和残差核验。

## 编译与检查

在本目录运行 PowerShell：

    .\build.ps1

脚本在 ASCII 临时目录完成 XeLaTeX、BibTeX、XeLaTeX、XeLaTeX，再复制 PDF 与 SyncTeX 回项目。可用 OutputDirectory 参数另选目录。

scripts 保留符号计算、留数系统、单轨道三场、一般多轨道样点及 Cauchy 表示的核验。有限样点不替代一般数据的解析证明。

calculations/literature-rewrite-review.md 记录跨来源文献与读取范围；calculations/structural-rewrite-review.md 记录重写结构及证明范围；PDF 验收见 calculations/pdf-acceptance.json。
