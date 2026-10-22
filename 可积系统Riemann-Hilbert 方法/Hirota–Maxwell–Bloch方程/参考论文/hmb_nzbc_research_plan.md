# Hirota–Maxwell–Bloch 方程非零边界 RHP：课题总结与学习计划

> 整理日期：2026 年 10 月 20 日  
> 定位：硕士毕业论文预选题；以 **非零光场背景 + 无反射 RHP + 孤子解** 为主体，以 **共振点相容条件** 为可选补充。  
> 文献状态：以下“尚未查到”的判断基于目前公开文献检索，不等于已证明不存在相关论文；正式开题前仍需查阅全文、参考文献及不同名称下的等价模型。

## 一、拟定题目与研究目标

**中文题目：** Hirota–Maxwell–Bloch 方程在非零边界条件下的 Riemann–Hilbert 问题及孤子解。

**英文题目：** *Riemann–Hilbert Problem and Soliton Solutions for the Hirota–Maxwell–Bloch System with Nonzero Boundary Conditions*.

考虑系统

$$
\begin{dcases}
q_z=\mathrm{i}\alpha\left(\dfrac{1}{2}q_{xx}+|q|^2q\right)+\beta(q_{xxx}+6|q|^2q_x)+2p,\\
p_x=2\eta q-2\mathrm{i}ap,\\
\eta_x=-(\overline q p+q\overline p),
\end{dcases}
$$

其中 $q,p$ 为复值函数，$\eta$ 为实值函数，$a,\alpha,\beta\in\mathbb R$。这里的**非零边界**专指光场 $q$ 在 $x\to\pm\infty$ 时具有非零振幅，不是仅有 $\eta\to\eta_0\ne0$。

可先固定最简单的相容平面波背景（$a\ne0$）：

$$
\begin{aligned}
q_{\mathrm{bg}}&=Qe^{-\mathrm{i}\Omega z},\\
p_{\mathrm{bg}}&=-\dfrac{\mathrm{i}\eta_0Q}{a}e^{-\mathrm{i}\Omega z},\\
\eta_{\mathrm{bg}}&=\eta_0,\qquad
\Omega=-\alpha Q^2+\dfrac{2\eta_0}{a},
\end{aligned}
$$

其中 $Q>0$、$\eta_0\ne0$，并研究趋于此背景的扰动。后续再视需要考虑两端具有不同常相位的边界。

**核心目标：** 借鉴范恩贵老师参考书第六章的非零边界聚焦 NLS 谱分析方法，构造 HMB 的单值化谱问题、Jost 解、散射矩阵与矩阵 RHP；在无反射条件下求简单离散谱对应的显式解，并同时恢复 $q,p,\eta$。

## 二、已有研究及其与选题的区别

### 已下载的四篇（按截图文件名）

| 论文 | 已有工作 | 重点阅读内容 |
| --- | --- | --- |
| **Li–Xia–Wei (2021)**, [*The N-soliton solutions to the Hirota and Maxwell–Bloch equation via the Riemann–Hilbert approach*](https://doi.org/10.1142/S0217979221501538) | HMB 的**零光场边界** RHP 与 $N$-孤子解 | Lax 对、Volterra 方程、RHP、无反射求解；用于熟悉原方程记号 |
| **Zhao–Zhaqilao (2026)**, [*Breathers and rogue waves on the periodic background in the Hirota and the Maxwell-Bloch equations*](https://doi.org/10.1016/j.wavemoti.2025.103696) | 同一 HMB 的平面波及周期背景；主要使用 Darboux 变换研究呼吸子、怪波等 | **平面波种子、完整 Lax 对、背景特征函数**；用于检查新解是否与已有 Darboux 解等价 |
| **Liu–Liu (2026)**, [*Inverse scattering transform and soliton solutions for the coupled Hirota-Maxwell-Bloch system with nonzero background*](https://doi.org/10.1007/s11071-026-12724-0) | 使用 IST/RHP，但其“非零背景”只施加在 $\eta$ 上，$q,p\to0$ | 比较两种“非零背景”的本质差别，学习介质变量的处理方式 |
| **Biondini–Gabitov–Kovačič–Li (2019)**, [*Inverse scattering transform for two-level systems with nonzero background*](https://doi.org/10.1063/1.5084720) | Maxwell–Bloch 系统的**非零光场** IST 和显式解；并非完整 HMB | 非零背景下介质状态的相容性、散射理论与共振问题；[预印本](https://arxiv.org/abs/1907.06231) |

**当前文献判断：** 零光场 HMB–RHP 已做；非零光场 HMB 的 Darboux 特殊解已做；普通 Hirota 的非零边界 RHP 已做。本次尚未查到与“**完整 HMB + 真正非零光场边界 + 标准无反射矩阵 RHP**”完全对应的论文。不可据此宣称绝对首创。

### 建议补充下载的两篇（优先）

1. **Li–He–Porsezian (2013)**, [*Rogue waves of the Hirota and the Maxwell-Bloch equations*](https://doi.org/10.1103/PhysRevE.87.012913)，*Physical Review E*, **87**, 012913；[arXiv 全文](https://arxiv.org/abs/1205.1191)。**用途：** 核对 HMB 的 Lax 对、Darboux/dressing 矩阵、共振点处 $p,\eta$ 的变换公式。
2. **Zhang–Chen–Yan (2020)**, [*Focusing and defocusing Hirota equations with non-zero boundary conditions: Inverse scattering transforms and soliton solutions*](https://doi.org/10.1016/j.cnsns.2019.104927)，*Communications in Nonlinear Science and Numerical Simulation*, **80**, 104927。**用途：** 直接参考非零边界普通 Hirota 方程的 Jost 解析性、单值化、散射对称性、无反射 RHP、trace formula 和 theta condition；区别在于该方程没有 Maxwell–Bloch 共振耦合。

**主要教材：** 范恩贵《可积系统、正交多项式和随机矩阵——Riemann–Hilbert 方法》，**第五章作零边界对照，第六章作为主要技术模板**。以上六篇加该教材，足以构成第一轮学习材料；暂时不必大量增加论文。

## 三、拟开展的工作与创新点

**核心工作 1：非零光场边界下的谱分析与 RHP。** 核对背景及 Lax 对，通过 $\rho(\lambda)=\sqrt{\lambda^2+Q^2}$ 引入双叶黎曼曲面，并使用单值化坐标 $w=\lambda+\rho(\lambda)$。建立 Jost 解、散射关系、解析区域、对称性、归一化条件，以及简单离散谱的留数条件。**潜在创新点在于将 HMB 的零光场 RHP 推进到真正非零光场边界。**

**核心工作 2：无反射 RHP 与孤子解。** 令反射系数恒为零，求解离散谱所对应的亚纯 RHP；至少完成一组简单离散谱的显式解，争取推广到两组或一般 $N$ 组。除光场 $q$ 外，还要根据共振矩阵恢复 $p,\eta$，并检验原方程、背景条件和 Bloch 约束。已有 Darboux/dressing 的有限极点构造**不等于**完整 NZBC 散射 RHP 已经建立，必须补充二者的严格对应。

**可选补充工作 3：共振点的散射相容条件。** HMB 的时间 Lax 矩阵在实谱参数 $\lambda=a$ 处存在极点。由介质方程可得

$$
R_x=[U(a),R],\qquad
R=\begin{pmatrix}
\eta&-p\\
-\overline p&-\eta
\end{pmatrix}.
$$

在**两端介质状态相容、Jost 基底在 $a$ 处正则等明确假设**下，可尝试证明散射矩阵在 $a$ 处的非对角元为零，进而得到必要条件 $r(a)=0$。这仍是**待严证和查重的候选命题**；因为主课题采用 $r\equiv0$，它并非无反射情形自动产生的独立新结论，不应在开题时写作“已证明的创新成果”。

## 四、边界与难点

- **必须完成：** 背景相容性、非零背景 Jost 解和单值化、正确的 RHP 及简单极点条件、至少一个无反射解、$q,p,\eta$ 的重构与验证。
- **重点核查：** 时间 Lax 矩阵的共振点 $\lambda=a$ 位于连续谱上；背景规范化、时间相位、介质背景以及离散谱对称性必须一致。
- **不作为毕业论文必做：** 任意非零反射数据的完整 IST 存在唯一性、Deift–Zhou 长时间渐近、双重极点、怪波极限。
- **创新性边界：** 不宣称首次得到 HMB 的非零背景孤子；强调新边界条件下的标准散射 RHP 构造及介质变量的严格重构。正式确定创新点前还需排查等价 AKNS 层级和其他命名下的已有文献。

## 五、简要学习与研究计划

| 阶段 | 学习材料 | 本阶段任务与成果 |
| --- | --- | --- |
| **1. 熟悉 HMB 的零边界 RHP**（约 1–2 周） | 范老师第五章；**Li–Xia–Wei (2021)** | 看懂方程对应的 Lax 对、Volterra/Jost 解、散射矩阵与无反射解法；整理一页符号对照表。 |
| **2. 确认非零背景与谱结构**（约 2–3 周） | 范老师第六章；**Zhao–Zhaqilao (2026)**；补充 **Li–He–Porsezian (2013)** | 验证平面波背景，求背景特征值与单值化变量，核对共振项及介质约束。先完成技术可行性检查。 |
| **3. 构造自己的非零边界 RHP**（约 3–5 周） | 补充 **Zhang–Chen–Yan (2020)**；参考 **Biondini 等 (2019)** | 推导 Jost 解析性、谱对称性、跳跃关系、留数条件；明确哪些公式来自 NLS/Hirota，哪些必须重新证明。 |
| **4. 无反射解与论文整理**（约 3–4 周） | 回看六篇论文；**Liu–Liu (2026)** 用于对照 | 求一组及两组简单离散谱解，重构三个变量并验证；视时间研究 $r(a)=0$，撰写毕业论文初稿。 |

以上时间是**初步估计**，各阶段可以重叠。最重要的阶段验收是：**在阶段 2 结束时确认背景时间演化与共振点可处理，再决定是否正式按完整 RHP 路线推进**。

## 六、当前最直接的待办事项

1. **先读已有的 2021 年 HMB 零边界 RHP**，同时对照范老师第五章；不用从头细读所有其他论文。
2. **补充下载 2013 年 HMB Darboux/Lax 对论文、2020 年 Hirota 非零边界 RHP 论文**。
3. 接着以范老师第六章为主线，验证自己选定背景的谱结构、单值化及共振项；完成这一步后再开始正式推导 RHP。

> **一句话定位：** 在已有 HMB 的零光场 RHP 与非零光场 Darboux 结果之间，补上“真正非零光场边界下的标准无反射 RHP 与完整介质重构”这一环节；共振点 $r(a)=0$ 的必要条件作为可选补充，而不是主任务。
