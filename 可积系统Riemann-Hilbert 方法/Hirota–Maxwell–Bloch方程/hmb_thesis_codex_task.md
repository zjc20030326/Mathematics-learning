# Codex 任务书：HMB 方程非零光场边界 RHP 硕士论文初稿

> 文档性质：供 Codex 在本地工作区直接执行的完整任务说明。目标是生成一套**可以编译、可以阅读、数学状态诚实标注**的复旦硕士学位论文 LaTeX 初稿，用来预览结构并检验研究路线；**不是**要求伪造已经完成的科研证明。
>
> 研究日期基准：2026-10-21。论文可以不投稿；希望与现有研究不完全重复。本文档中的“尚未检索到”绝不等于严格确认“世界上无人研究过”。

## 0. 总指令及交付标准

你是辅助数学研究与 LaTeX 写作的 Codex。请**实际创建/编辑文件、核验数学步骤并尽可能编译生成 PDF**，不要只回复建议或论文提纲。最终至少交付：

1. 一套可独立编译的中文硕士论文 `.tex` 工程，沿用指定模板；包含中英文摘要、目录、正文、参考文献。首版内容要尽可能充实，而非只有目录和空白章节。
2. 可编译时输出 `fduthesis.pdf`（或按模板实际主文件同名生成的 PDF），记录使用的编译命令和结果。若没有 TeX 环境，保留完整源码并明确说明未编译，**不要声称成功**。
3. 工程内新增 `RESEARCH_STATUS.md`：按章节/定理记录**已由原始文献核实、已独立证明、仅形式推导、尚待证明、需要文献确认**的状态；将所有当前不能严格完成的数学步骤列为待办清单。
4. `REFERENCES.md`：各论文准确书目信息、DOI、主要用途、核查来源、是否拿到全文及准确的式号/页码（拿不到全文不要虚构式号）。
5. `BUILD.md`：环境依赖、构建命令、错误排除和生成文件路径；若新增绘图脚本，说明复现方式。

**执行优先顺序：数学正确性与来源真实性 > 核心内容完整性 > 首版 PDF 可编译 > 排版精细化。**

若使用现有工作区，先检查 `git status` 和适用的 `AGENTS.md`，不要删除、覆盖用户的无关文件。优先新建独立目录 `hmb-nzbc-thesis/`；不主动推送或提交远程 GitHub。可以分阶段多次编译，保持可恢复进度。遇到难题先尝试独立研究、计算核验，再把无法完成的环节写进状态报告；不要为了让文章看起来完整而编造证明。

## 1. 模板及数学排版规范

- 指定模板仓库：<https://github.com/yzbrlan/fudan-thesis-latex-template.git>。
- **使用 `graduation-thesis/`，不要使用 `course-thesis/`。** 已核实其主文件为 `graduation-thesis/fduthesis.tex`，默认包含 `\documentclass[type=master,oneside]{fduthesis}`、`\fdusetup`、章节文件 `section/*.tex` 及 `citelist.bib`。
- 优先克隆或复制其完整毕业论文目录，再改写主文件的正文 `\chapter` 与 `\input` 部分。保留模板需要的 `.cls`、`.def`、`.sty`、图片与引文依赖，不要从零手写简陋的替代模板。
- 模板 `fduthesis.tex` 实际配置了 `bib-backend=bibtex`、`bib-resource={citelist.bib}`、`\bibliography{citelist}`。首版优先按 **XeLaTeX → BibTeX → XeLaTeX → XeLaTeX** 编译，不能在未检查类文件的情况下把 BibTeX 与 biber 混用。系统字体缺失时选用已安装的 CJK 字体（例如 Windows 环境可尝试 `cjk-font=windows`），不盲目指定不存在的字体。
- 若工作区有 `AGENTS.md`，优先执行其指令。数学 LaTeX 风格参考用户仓库 <https://github.com/zjc20030326/mathematics-competition> 的当前版本：`.agents/skills/writing-mathematical-latex/SKILL.md` 及 `references/math-latex-style.md`，先读取再生成复杂数学源码；若无法读取，说明原因并按照以下基本约定工作。
- `.tex` 中行内公式用 `$...$`，独立公式用 `\[...\]` 或 `equation/align` 环境；正文不使用 `$$...$$`。普通分式 `\dfrac{...}{...}`，上下标/指数中分式 `\frac{...}{...}`，保持括号、复共轭和上下标一致。数学 Markdown 说明文件中行内公式使用 `$...$`，行间公式使用独占行的 `$$...$$`。
- 中文论文以数学研究论文的表达为准：定理、引理、命题、证明清晰分层；关键推导逐步列出依据；对借用文献结论标注引用，对自推结论明确适用条件。定理引用、公式标签和参考文献交叉引用必须真实。
- 此版仅试排：姓名、学号、导师、院系、致谢等暂不确定的个人信息请填写“待填写”并在 `RESEARCH_STATUS.md` 注明，**不可猜测**；保留可编译的模板占位封面数据。

## 2. 拟定论文题目与范围

中文：**《Hirota–Maxwell–Bloch 方程在非零边界条件下的 Riemann–Hilbert 问题及孤子解》**。

英文：*Riemann–Hilbert Problem and Soliton Solutions for the Hirota–Maxwell–Bloch System with Nonzero Boundary Conditions*。

为了与已有文献准确区分，正文应反复强调研究的是**光场 $q$ 的真正非零背景**（nonzero optical-field boundary），而不是仅有介质变量 $\eta$ 非零。核心任务集中于非零背景空间谱、Jost 解与散射矩阵、RHP 的构造、无反射简单离散谱解以及完整重构 $q,p,\eta$。共振条件 $r(a)=0$ 为附加研究，允许未完成。

不将一般含反射 RHP 的全局可解性、Deift–Zhou 长时间渐近、双重极点、怪波极限、无限维适定性作为首版必做内容。

## 3. 研究方程与约定（必须原样保留数学含义）

取 $q=q(x,z),p=p(x,z)\in\mathbb C$，$\eta=\eta(x,z)\in\mathbb R$，参数 $\alpha,\beta,a\in\mathbb R$，其中 $x$ 为推迟时间，$z$ 为传播变量。方程组为：

$$
\begin{cases}
q_z=\mathrm{i}\alpha\left(\dfrac{1}{2}q_{xx}+|q|^2q\right)
+\beta\left(q_{xxx}+6|q|^2q_x\right)+2p,\\
p_x=2\eta q-2\mathrm{i}ap,\\
\eta_x=-(\overline q p+q\overline p).
\end{cases}
\tag{HMB}
$$

### 3.1 初步选定的非零平面波背景

为简化起步，先取相同左右背景相位的情形，并假设 $Q>0$、$a\ne0$、$\eta_0\ne0$：

$$
\begin{aligned}
q_{\mathrm{bg}}(x,z)&=Qe^{-\mathrm{i}\Omega z},\\
p_{\mathrm{bg}}(x,z)&=-\dfrac{\mathrm{i}\eta_0Q}{a}e^{-\mathrm{i}\Omega z},\\
\eta_{\mathrm{bg}}(x,z)&=\eta_0,\\
\Omega&=-\alpha Q^2+\dfrac{2\eta_0}{a}.
\end{aligned}
\tag{BG}
$$

应**逐项直接代回** HMB 验证，而不是只说“容易验证”。先研究扰动满足 $q-q_{\mathrm{bg}}\to0$、$p-p_{\mathrm{bg}}\to0$、$\eta-\eta_0\to0$（$x\to\pm\infty$）的情形；后续可以探讨两端相位 $\phi_\pm$ 不同但模相同的条件及 theta condition。时间演化的相位不能从 NLS 直接照抄。

### 3.2 Bloch 约束（应独立证明）

由后两个方程严格推得：

$$
\partial_x\left(|p|^2+\eta^2\right)=0,
\qquad
|p(x,z)|^2+\eta(x,z)^2=C(z).
\tag{B}
$$

如采用归一化 $C(z)=1$，必须额外说明与背景参数的关系，不得暗自将 $Q,a,\eta_0$ 视为完全独立。注意 $q\ne0$ 但 $p=\eta=0$ 时会退化为普通 Hirota 演化；本课题应保持 Maxwell–Bloch 耦合的非平凡性。

## 4. Lax 对、空间谱与单值化：已知线索和需验证的计算

可采用文献中的 AKNS 型空间 Lax 对（若原文谱参数的平移、缩放或符号不同，请先做精确转换，不能混用）：

$$
\Phi_x=U(\lambda)\Phi,\qquad
U(\lambda)=-\mathrm{i}\lambda\sigma_3+
\begin{pmatrix}0&q\\-\overline q&0\end{pmatrix},
\qquad \sigma_3=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\tag{L1}
$$

在适当的时间相位规范变换后，对相同相位背景，预期得到：

$$
U_\infty(\lambda)=\begin{pmatrix}-\mathrm{i}\lambda&Q\\-Q&\mathrm{i}\lambda\end{pmatrix},
\qquad U_\infty(\lambda)^2=-(\lambda^2+Q^2)I.
\tag{S1}
$$

引入 $\rho(\lambda)=\sqrt{\lambda^2+Q^2}$，背景特征值 $\pm\mathrm{i}\rho(\lambda)$，分支点 $\pm\mathrm{i}Q$。预期连续谱为 $\mathbb R\cup\mathrm{i}[-Q,Q]$（须在采用的谱约定下重新证明）。单值化参数：

$$
w=\lambda+\rho(\lambda),\qquad
\lambda=\dfrac{1}{2}\left(w-\dfrac{Q^2}{w}\right),\qquad
\rho=\dfrac{1}{2}\left(w+\dfrac{Q^2}{w}\right).
\tag{S2}
$$

预期跳跃轮廓变为 $\mathbb R\cup\{w:|w|=Q\}$；需明确两叶的映射、定向、谱对称、解析区域、边界值和 $w=0,\infty$ 两处归一化。不得直接把普通 NLS 的 RHP 结果照抄为 HMB 定理。

### 4.1 时间 Lax 对与共振点

现有 HMB 的时间 Lax 矩阵具有如下结构：

$$
\Phi_z=V(\lambda)\Phi,\qquad
V(\lambda)=V_{\mathrm{poly}}(\lambda)+\dfrac{\mathrm{i}}{\lambda-a}R,
\qquad
R=\begin{pmatrix}\eta&-p\\-\overline p&-\eta\end{pmatrix}.
\tag{L2}
$$

$V_{\mathrm{poly}}$ 最高三次，但**必须以原论文完整 Lax 对为准**逐项抄录、变换和检验。由 (HMB) 可直接核实 $R_x=[U(a),R]$；这一等式应详细证明，并用零曲率条件再次核对。时间共振极点为 $\lambda=a\in\mathbb R$，与非零背景空间连续谱相交；它不是普通离散特征值。

此前讨论中曾得到以下**待重新独立验算**的背景时间矩阵公式（只作为校验候选，不得不经检验认定正确）：

$$
\widehat V_\infty(\lambda)=-F(\lambda)U_\infty(\lambda),\qquad
F(\lambda)=4\beta\lambda^2-\alpha\lambda-2\beta Q^2
+\dfrac{\eta_0}{a(\lambda-a)}.
\tag{TC}
$$

请从所采用的完整 $V$ 和规范变换 $\Phi=G(z)\Psi$ 出发重算该式，核实所有正负号、系数、规范相位、谱参数和零曲率条件；如不一致，以自洽推导为准并在状态报告中说明。若成立，再定义背景振荡相位 $\theta(\lambda;x,z)=\rho(\lambda)[x-F(\lambda)z]$；若不成立，必须相应修改。

### 4.2 必须澄清的几个谱问题

- 在单值化谱平面确切确定 Jost 解哪些列在哪些区域解析及其 Volterra 积分方程，证明所用积分核和背景特征向量的性质。
- 时间演化中共振点 $\lambda=a$ 在 $w$ 平面的原像为 $w_\pm=a\pm\sqrt{a^2+Q^2}$，位于实谱轴上；不得因单值化而假设共振极点消失。
- 写出散射矩阵的对称性、行列式约束、零点结构；检查非零背景下典型的谱点轨道是否包含共轭与 $w\mapsto-Q^2/w$ 等变换。确切对称关系须推导，不能凭记忆假定。
- 需要在 $w\to\infty$ 与 $w\to0$ 两种极限中分别恢复正确的背景行为，防止只有一种归一化。

## 5. 核心研究任务（必须尽可能实际推导）

| 研究内容 | 优先级 | 首版目标 |
| --- | --- | --- |
| 非零背景谱问题、Jost 解、散射矩阵及 RHP 构造 | 核心 | 以定理、引理、证明的形式完整推导可验证部分 |
| 无反射 RHP、简单离散极点 | 核心 | 至少求一组离散谱对应的非平凡解，争取二孤子/一般 $N$-孤子 |
| $q,p,\eta$ 的重构与原 PDE 验证 | 核心 | 三个变量都必须重构，严格或通过独立代数验证其满足 HMB |
| 共振点相容条件 $r(a)=0$ | 补充 | 在清楚列出的背景、Jost 正则性及散射矩阵定义下尝试证明；不能证明则明确列为条件命题 |

### 5.1 无反射 RHP 的目标

建立具有正确谱对称性、归一化、离散谱留数条件的矩阵 RHP。将反射系数置零，求解相应的亚纯逆问题。需要处理的不是随意一个单点极点，而是满足物理约化及非零背景对称性的**完整谱点集合**。

可以参考先前讨论的一次 dressing 候选式：若 $\lambda_1\in\mathbb C^+$ 且离开连续谱、分支点和共振点，取背景 Lax 解 $y\ne0$，定义

$$
P=\dfrac{yy^\dagger}{y^\dagger y},\qquad
D(\lambda)=I-\dfrac{\lambda_1-\overline{\lambda_1}}{\lambda-\overline{\lambda_1}}P,
\qquad P^2=P=P^\dagger.
\tag{D1}
$$

它满足 $D\to I$、$D^{-1}(\lambda)=I+\dfrac{\lambda_1-\overline{\lambda_1}}{\lambda-\lambda_1}P$；由于 $a\in\mathbb R$ 且 $\operatorname{Im}\lambda_1>0$，$D(a)$ 存在且可逆。可尝试利用

$$
R^{[1]}=D(a)R^{[0]}D(a)^{-1}
\tag{D2}
$$

重构介质变量，并从空间 Lax 对的无穷远展开恢复 $q$。**(D1)–(D2) 只是已有 Darboux/dressing 构造的技术线索，不可把它们冒充“已经完成标准非零边界散射 RHP 的证明”。** 必须检查它与单值化谱平面上的极点轨道、留数条件、theta condition、对称性及边界相位是否一致；未完成此处则标记为“尚待证明”。

### 5.2 解的验证

- 明确 $q^{[1]}$, $p^{[1]}$, $\eta^{[1]}$ 的可计算表达式。若得到一般 $N$-孤子行列式或线性系统表达式，应先严格解释符号并检查 $N=1$。
- 直接代回 HMB 三个方程；如果代数巨大，可以提供 SymPy/Mathematica 脚本，分别计算三条 PDE 残差，示例参数必须满足背景与离散谱条件。
- 核对 $\partial_x(|p|^2+\eta^2)=0$、非零背景渐近、极点无碰撞及物理约化 $\eta\in\mathbb R$。
- 和 Zhao–Zhaqilao（2026）、Li–He–Porsezian（2013）的 Darboux 解作合理比较；若解在适当变换下等价，应如实说明，创新点是 RHP/逆散射表述，而非新的波形。
- 如有绘图，生成 $|q|$、$|p|$、$\eta$ 的示例，并说明参数与图对应的数学结果；绘图不是数学证明。

## 6. 补充课题：共振点的反射条件（必须注意逻辑状态）

假设 $\lambda=a\in\mathbb R$ 位于连续谱，且介质背景与光场相容。由 (L2) 定义的 $R$ 满足

$$
R_x=[U(a),R],\qquad
\partial_x\left(\Phi(x,a)^{-1}R(x)\Phi(x,a)\right)=0.
\tag{RC1}
$$

此前曾在**两端介质归一化一致、Jost 基底在 $a$ 处正则可逆、两端谱对角化匹配且散射矩阵在该点可定义**等假设下，形式上推得 $S_{12}(a)=S_{21}(a)=0$，继而 $r(a)=0$。这是一条**需要独立复核和严格列明假设的候选命题**，不是可以不经证明直接引用的文献定理。

Codex 应：

1. 从 $R_x=[U(a),R]$ 证明积分常数矩阵关系。
2. 分别讨论左右背景相位与介质末态，核实何时两端的共轭常数矩阵确实相同。
3. 在正确定义的 Jost 归一化下，推导 $S(a)$ 必须与何种非退化对角矩阵对易；检验该推导是否足以推出 $r(a)=0$。
4. 澄清即使 $r(a)=0$，含有 $\exp(c/(\lambda-a))$ 的散射数据演化仍可能破坏微分正则性；因此**不能**从 $r(a)=0$ 直接推出一般含反射 IST 完全可解。
5. 无法严格完成时，在正文写“条件下的讨论/有待验证的命题”，并将完整的推导障碍写入 `RESEARCH_STATUS.md`。不要虚构一般 RHP 存在唯一性证明。

## 7. 论文结构建议（按数学论文重组模板）

建议正文六章左右，按照实际成果灵活拆并：

1. **绪论**：研究背景、HMB 与 NLS/Hirota/Maxwell–Bloch 的联系、已有 IST/Darboux/RHP 文献、研究目标及创新点的准确表述。
2. **模型、Lax 对与非零背景谱分析**：方程、相容的平面波、Bloch 约束、背景规范变换、黎曼曲面、单值化变量、谱对称性。
3. **Jost 解和直接散射理论**：积分方程、解析性、散射系数和对称性、时间演化、共振点处的必要限制。
4. **非零边界矩阵 RHP**：未知矩阵、解析区域、跳跃矩阵、两端归一化、极点与留数、迹公式或 theta condition（依实际可证明情况）。
5. **无反射 RHP 与孤子解**：简单离散谱、有限维代数问题、一孤子及可能的 $N$-孤子、$q,p,\eta$ 重构、结果验证/图像。
6. **总结与展望**：严格区分已完成证明的结论、条件性的讨论以及一般反射数据/高阶极点等未来问题。

可在附录收集长矩阵计算、单极点构造与验证代码说明。参考范恩贵《可积系统、正交多项式和随机矩阵——Riemann–Hilbert 方法》第五章（零边界聚焦 NLS）与第六章（非零边界聚焦 NLS）作为方法模板；具体定理引用时须核对用户持有的 PDF 版本和编号。

**写作要求**：正文第一版应当在可证明的模型背景、Bloch 约束、谱特征值、单值化恒等式、基本 dressing 恒等式等方面展示完整推导。还未建立的 Jost 解析性、正规 RHP 与 $N$-孤子定理，不能“空口声称”完成；可以把真正尚待解决的技术点合理写入末尾的“进一步研究”及 `RESEARCH_STATUS.md`，同时保证首版论文仍然是可读的数学初稿，而不是堆砌占位符。

## 8. 已有论文与需要补充的论文

### 8.1 用户已经下载的四篇（尽量优先读取本地 PDF）

| 简称 | 准确标题 | DOI / 链接 | 主要用途 |
| --- | --- | --- | --- |
| Li–Xia–Wei (2021) | *The N-soliton solutions to the Hirota and Maxwell–Bloch equation via the Riemann–Hilbert approach* | <https://doi.org/10.1142/S0217979221501538> | **HMB、零光场边界、RHP、无反射 $N$-孤子**；参考空间 Lax 对和零边界流程 |
| Zhao–Zhaqilao (2026) | *Breathers and rogue waves on the periodic background in the Hirota and the Maxwell-Bloch equations* | <https://doi.org/10.1016/j.wavemoti.2025.103696> | **同一类 HMB、非零平面波及周期背景、Darboux**；参考背景与精确解，但不能说它已做完整 NZBC RHP |
| Liu–Liu (2026) | *Inverse scattering transform and soliton solutions for the coupled Hirota-Maxwell-Bloch system with nonzero background* | <https://doi.org/10.1007/s11071-026-12724-0> | 其所谓 nonzero background 施加在 $\eta$ 上，$q,p\to0$；用于区分研究边界与排除重复 |
| Biondini–Gabitov–Kovačič–Li (2019) | *Inverse scattering transform for two-level systems with nonzero background* | <https://doi.org/10.1063/1.5084720>；<https://arxiv.org/abs/1907.06231> | **Maxwell–Bloch、真正非零光场、IST**；处理介质背景相容性与共振谱问题 |

本地截图所见文件名（Codex 可据此寻找文件，但不能假定它们一定已挂载在工作区）：

- `Breathers and rogue waves on the periodic background in the Hirota and the Maxwell-Bloch equations.pdf`
- `Inverse scattering transform for two-level systems with nonzero background.pdf`
- `Inverse scattering transform and soliton solutions for the coupled Hirota-Maxwell-Bloch system with nonzero background.pdf`
- `the-n-soliton-solutions-to-the-hirota-and-maxwell-bloch-equation-via-the-riemann-hilbert-approach.pdf`

### 8.2 建议必须补充的两篇

- **Li–He–Porsezian (2013)**, *Rogue waves of the Hirota and the Maxwell-Bloch equations*, *Physical Review E* **87**, 012913. DOI <https://doi.org/10.1103/PhysRevE.87.012913>，arXiv <https://arxiv.org/abs/1205.1191>。核心：HMB Lax 对、$n$ 重 Darboux、通过共振点共轭重构介质变量。请同时注意该刊 2013 年的出版者更正。
- **Zhang–Chen–Yan (2020)**, *Focusing and defocusing Hirota equations with non-zero boundary conditions: Inverse scattering transforms and soliton solutions*, *Communications in Nonlinear Science and Numerical Simulation* **80**, 104927. DOI <https://doi.org/10.1016/j.cnsns.2019.104927>。核心：**普通 Hirota 真正 NZBC 下的 IST/RHP、简单/双重极点、trace formula、theta condition**，可与范老师第六章对照；其模型不含 $p,\eta$，不可直接充作 HMB 的论文。

### 8.3 可选背景原始文献

- **Porsezian–Nakkeeran (1995)**, *Optical Soliton Propagation in an Erbium Doped Nonlinear Light Guide with Higher Order Dispersion*, *Physical Review Letters* **74**, 2941–2944. DOI <https://doi.org/10.1103/PhysRevLett.74.2941>。用于模型起源及物理背景。
- 范恩贵：《可积系统、正交多项式和随机矩阵——Riemann–Hilbert 方法》。用户持有的 PDF 版本可能与其他版本符号、页码不同；**未经核对不得虚构版次、出版年份、页码或公式编号**。

请对所有参考文献实际建立 BibTeX 条目，作者、标题、年份、期刊、卷号、文章编号及 DOI 逐一核验；获取不到全文时只能引用已验证的元数据和公开摘要结论，不能捏造证明细节。

## 9. 已知与未证的“科研状态”边界

必须准确区分：

- **已知既有工作**：HMB 的零光场无反射 RHP；HMB 非零光场背景上的 Darboux 特殊解；普通 Hirota 的真正非零边界 IST/RHP；Maxwell–Bloch 非零光场 IST。
- **可直接独立检验的初步数学事实**：背景 (BG) 是原 HMB 的精确解、Bloch 约束 (B)、空间矩阵平方恒等式 (S1)、单值化代数恒等式 (S2)、由介质方程推得的 (RC1) 第一式。
- **需从原文重算并核对规范的计算**：时间 Lax 矩阵全部系数、(TC) 的 $F(\lambda)$、背景指数相位、各谱点对称性、所有 RHP 留数公式和介质重构系数。
- **可以尝试而不能预先宣布完成的核心成果**：在单值化平面严格建立标准 NZBC RHP；将已知 dressing 因子与真正的无反射散射 RHP 一一对应；构造符合全部谱对称性和背景条件的孤子解；完成 $q,p,\eta$ 的严格重构。
- **条件性补充命题**：在相容介质端态、共振点 Jost 正则性等假设下推导 $r(a)=0$；不应在无反射问题中将其称为额外独立的非平凡结果（因为 $r\equiv0$ 已包含此条件）。
- **未解决且本稿不承诺的内容**：所有一般含反射散射数据的时间演化正则性和完整 IST/RHP 存在唯一性；共振点本性振荡可能产生的复杂局部分析。

禁止把以前聊天中出现过的“高精度数值残差已验证”等陈述当成你已经实际执行的结果。若本轮没有脚本输出、参数和程序，必须重新计算，才能在正文中声称进行了数值验证。

## 10. 推荐实际执行流程

**阶段 A：模板与文献核对**。检查本地文件、克隆模板、读取目录与模板 README、读取数学 LaTeX Skill；索引本地四篇 PDF；下载/读取另外两篇（如果有网络权限）。写出文献对照与状态表，确保题目边界没有混淆。

**阶段 B：先出可编译论文骨架**。修改真实毕业论文模板，建立数学论文目录、中英文摘要、参考文献；真实的模型和文献综述必须写入，不得只生成标题。先尝试编译一次，解决 CJK、交叉引用和 BibTeX 的问题。

**阶段 C：核心数学内容逐节填充**。从最容易严格验证的 (BG)、(B)、(S1)–(S2) 开始；然后精确核对 Lax 对和相位；依次处理 Jost、RHP、简单谱点、无反射解、三场重构。每一步有原文/计算来源，按定理—证明组织。

**阶段 D：验收与可选数值验证**。检查原 PDE 三式、背景渐近、Hermitian 约束、谱对称性和图像结果。运行实际可用的 LaTeX 编译及 BibTeX 检查；输出 `RESEARCH_STATUS.md` 和 `BUILD.md`。

最终回复应包含：生成/修改的文件清单、LaTeX 编译状态、主要已证结论、尚未完成的关键数学证明、下一步最值得解决的三项具体技术问题。不要把待证命题写成“本文已经证明”。

---

**一句话验收标准：** 我要的是**可以先打开 PDF 预览的真实硕士论文数学初稿**，并且能够清楚看出哪些推导确实成立、哪些是尚待研究的核心课题；不需要一次性达到可送审的学位论文质量，更不能通过虚构证明来追求表面完整。
