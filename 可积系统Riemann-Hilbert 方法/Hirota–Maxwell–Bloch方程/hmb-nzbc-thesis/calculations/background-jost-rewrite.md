# 背景谱与 Jost 解的推导重排

日期：2026-11-15。修改直接位于原 main 分支，以修改前工作区快照为基准，保留此前未提交的修改。

## 参考原文与读取范围

Zechuan Zhang、Engui Fan，Inverse scattering transform for the Gerdjikov–Ivanov equation with nonzero boundary conditions，ZAMP 71:149 (2020)，DOI：10.1007/s00033-020-01371-z。

核对了[作者提供的原文](https://zechuanzhang.github.io/files/Zhang-Fan2020_Article_InverseScatteringTransformForT.pdf)第2、3节，PDF第2—7页；本地副本为参考论文目录中的 Zhang-Fan2020_Article_InverseScatteringTransformForT.pdf。
对应原文位置为 paper:doi:10.1007/s00033-020-01371-z#Sec-2 与 paper:doi:10.1007/s00033-020-01371-z#Sec-3。

参考的组织是：背景矩阵与特征值、双叶曲面、单值化；随后由特征向量、对角化、背景基本解、修正方程及全微分形式引出 Volterra 方程，最后证明解析性。对应原文式(2.8)—(2.14)、(3.1)—(3.21)。

## 正文安排

| 位置 | 当前内容 |
| --- | --- |
| 2.3开头 | 仅列出 $U_\infty,Q_\infty$ 与 $V_\infty,F$ 两组背景公式 |
| 2.3后部 | 从特征方程依次得到谱根、双叶曲面、单值化、振荡谱、解析区域、换叶映射及共振原像 |
| 3.1前部 | 分别求两个特征向量，组成 $Y$，计算对角化关系、行列式与逆矩阵 |
| 3.1中部 | 求解渐近谱问题；对完整Lax解推导全微分恒等式；由固定时刻的空间方程引出Jost归一化与Volterra方程 |
| 定理3.1 | 有界空间Jost解的存在唯一性、解析性及边界连续性 |
| 第4节共振相容性之前 | 在首次需要时引入介质背景矩阵 $R_\infty$ |

空间Jost解不再单列为定义环境。讲解版的特征向量、行列式和背景解计算随正文移到第3节；介质背景的四元计算移到第4节。相应标签沿用，公式索引按新顺序生成。

## HMB 中保留的数学结构

- HMB 的背景特征值为 $\pm\mathrm{i}\lambda$，解析域按 $\operatorname{Im}\lambda$ 定号，谱边界为实轴与圆周。GI 的背景指数涉及 $k\lambda$，其谱边界不用于本稿。
- 全微分形式以矩阵同时满足完整Lax对为条件，并固定谱点；时间部分避开两个共振原像。
- 空间Jost解由固定 $t$ 的积分方程独立构造，其谱域包含共振原像。第4节继续独立证明构造解的时间方程。
- 定理3.1的四列方向、积分核界、Neumann阶乘估计与唯一性证明均保留；分支点继续使用加权条件及退化核极限处理。
- $R_\infty$ 在共振共轭关系与第7节三场重构中仍有使用，因此在第4节引入。

## 核查

[精确矩阵核验脚本](../scripts/verify_background_jost.py)独立检查谱曲线、特征向量、对角化、逆矩阵、相位导数及全微分的两个系数。结果与源文件哈希见 [矩阵核验数据](background-jost-identities.json)。两册公式编号、源文件同步与PDF版面检查以当前验收记录为准。
