# 第三章：Jost 解、解析性、散射矩阵与共振相容性

对照正文：../hmb-nzbc-thesis/section/03-scattering.tex  
覆盖 §3.1–§3.4，特别是 thm:jost、prop:scattering、prop:quartet、prop:evolution、prop:resonance、prop:highenergy。

**学习目标**：能从二阶 AKNS 系统独立写出 Volterra 方程，判断四个 Jost 列的解析域，使用 Cramer 法则求出散射系数，并说明共振点反射为零的具体前提。

## 3.1 第一段：为什么先固定 z？

HMB 的 x 方程用于定义散射数据，z 方程用于决定散射数据如何传播。最初构造 Jost 解时，我们可以把 z 固定为一个参数，研究

$$
\chi_x=\widehat U(x,z,\lambda)\chi,
\qquad
\widehat U=-\mathrm{i}\lambda\sigma_3+A(u).
$$

当 x→±∞ 时 u→Q，故 \widehat U→U∞。把两者之差记为

$$
\Delta A(x)=A(u(x))-A(Q)
=\begin{pmatrix}0&u-Q\\-(\overline u-Q)&0\end{pmatrix}.
$$

若 u−Q∈L¹∩L∞，这个扰动沿实 x 可积。L¹ 保证 Volterra 积分的总强度有限；L∞ 帮助控制局部乘积和矩阵核。

注意这时**不是假设 q 本身在无穷远趋零**，而是减掉常数 Q 后的扰动可积。这是非零背景与第五章零背景的根本差别。

## 3.2 Jost 归一化为什么要乘一个指数？

背景常系数系统的基本解为 E(w)e^{−iρxσ₃}。所以自然定义从左、右逼近背景的空间基本解：

$$
\chi_\pm(x,w)e^{\mathrm{i}\rho x\sigma_3}
\longrightarrow E(w),\qquad x\to\pm\infty.
$$

把指数移到矩阵内部定义

$$
\mu_\pm=\chi_\pm e^{\mathrm{i}\rho x\sigma_3}.
$$

这样 μ± 的端态是 E，而不是快速振荡的 E e^{−iρxσ₃}；后续把左右两套解析列拼成 RHP 更方便。

下标 ± 表示归一化端点，不是取复平面的上下边界。μ_{+,1} 是右 Jost 矩阵第一列，μ_{-,2} 是左 Jost 矩阵第二列。

## 3.3 thm:jost 和 eq:volterra：完整的变分常数推导

定义 W±=E⁻¹χ±。利用 U∞E=−iρEσ₃，原方程给出

$$
\begin{aligned}
(W_\pm)_x
&=E^{-1}\chi_{\pm,x}\\
&=E^{-1}(U_\infty+\Delta A)E W_\pm\\
&=-\mathrm{i}\rho\sigma_3 W_\pm
+E^{-1}\Delta A E W_\pm.
\end{aligned}
$$

先在有限起点 x₀ 应用变分常数：

$$
W_\pm(x)
=e^{-\mathrm{i}\rho(x-x_0)\sigma_3}W_\pm(x_0)
+\int_{x_0}^{x}e^{-\mathrm{i}\rho(x-y)\sigma_3}
E^{-1}\Delta A(y)E W_\pm(y)\,\mathrm{d}y.
$$

左右乘回 E 和 e^{iρxσ₃}，再分别令 x₀→−∞ 或 +∞。最终得原文 eq:volterra：

$$
\boxed{
\begin{aligned}
\mu_\pm(x,w)=E(w)
+\int_{\pm\infty}^{x}
&E(w)e^{-\mathrm{i}\rho(x-y)\sigma_3}E(w)^{-1}\\
&\cdot\Delta A(y)\mu_\pm(y,w)
e^{\mathrm{i}\rho(x-y)\sigma_3}\,\mathrm{d}y.
\end{aligned}}
$$

为什么右侧最后还有一个右乘指数？因为 μ=χ e^{iρxσ₃}；右侧乘指数不可以与左侧矩阵自由交换。e^{iρxσ₃} 是对角矩阵，左、右两次指数属于不同乘法位置。

“从 +∞ 到 x 的积分”按

$$
\int_{+\infty}^{x}f(y)\,\mathrm{d}y
=-\int_x^{+\infty}f(y)\,\mathrm{d}y
$$

解释。这决定后面哪一个指数在积分区域内衰减。

## 3.4 eq:columnkernel：逐列计算指数有界性

先设

$$
\sigma_3=\operatorname{diag}(1,-1),\qquad
s_1=1,\quad s_2=-1.
$$

观察 μ 的第 j 列：右侧 e^{iρ(x−y)σ₃} 作用在 e_j 上，只给它乘一个标量 e^{i s_j ρ(x−y)}。于是核简化为

$$
K_j(x-y,w)=e^{\mathrm{i}s_j\rho(x-y)}
E e^{-\mathrm{i}\rho(x-y)\sigma_3}E^{-1}.
$$

对 j=1：

$$
e^{\mathrm{i}\rho(x-y)}
e^{-\mathrm{i}\rho(x-y)\sigma_3}
=\operatorname{diag}(1,e^{2\mathrm{i}\rho(x-y)}).
$$

对 j=2：

$$
e^{-\mathrm{i}\rho(x-y)}
e^{-\mathrm{i}\rho(x-y)\sigma_3}
=\operatorname{diag}(e^{-2\mathrm{i}\rho(x-y)},1).
$$

假设 w∈D₊，即 Imρ>0。

- 若选择负端归一化 χ₋，积分变量 y≤x，故 x−y≥0；对第一列 |e^{2iρ(x−y)}|≤1。**μ_{-,1} 在 D₊ 解析**。
- 若选择正端归一化 χ₊，积分变量 y≥x，故 x−y≤0；对第二列 |e^{−2iρ(x−y)}|≤1。**μ_{+,2} 在 D₊ 解析**。

如果 Imρ<0，上述指数不再按同一选择有界，需要交换列与端点。结论为

$$
\begin{array}{c|c}
D_+:\operatorname{Im}\rho>0 & \mu_{-,1},\ \mu_{+,2}\\
D_-:\operatorname{Im}\rho<0 & \mu_{+,1},\ \mu_{-,2}
\end{array}
$$

这里是帮助理解的二维表，不表示完整矩阵在整个 D₊ 或 D₋ 同时解析。

## 3.5 为何 Neumann 级数收敛？1/n! 从何而来？

以负端第一列为例，设积分算子为 K，则

$$
\mu_{-,1}=E_1+KE_1+K^2E_1+\cdots.
$$

第 n 次迭代的积分区域不是整个 n 维立方体，而是有序单纯形

$$
-\infty<y_n<\cdots<y_2<y_1<x.
$$

若核有界于 C_K∥ΔA(y)∥，由积分对称性得到

$$
\begin{aligned}
\|(K^nE_1)(x)\|
&\le\|E_1\|C_K^n
\int_{y_n<\cdots<y_1<x}
\prod_{\ell=1}^{n}\|\Delta A(y_\ell)\|
\,\mathrm{d}y_1\cdots\mathrm{d}y_n\\
&=\dfrac{\|E_1\|C_K^n}{n!}
\left(\int_{-\infty}^{x}\|\Delta A(y)\|\,\mathrm{d}y\right)^n.
\end{aligned}
$$

最后一个等号因为对称被积函数在 n! 个变量排列区域上的积分相等。有阶乘分母，级数绝对收敛；同样估计作用于两解之差，说明唯一性。每次迭代项都随 w 解析，紧集一致收敛由 Weierstrass 定理保证解析性。

**关键限制**：E⁻¹ 含 1/γ(w)，在 w=±iQ 退化，不能直接把此证明延伸到分支点。分支点可能需要带权的空间衰减假设与单独的阈值展开。

## 3.6 eq:Sdef：散射矩阵究竟表示什么？

χ₊、χ₋ 是同一二阶空间系统的两套可逆基本矩阵，所以

$$
S_{\mathrm{sp}}(w,z)=\chi_+(x,w)^{-1}\chi_-(x,w)
$$

与 x 无关。原文约定

$$
\chi_-=\chi_+S_{\mathrm{sp}}.
$$

这是一个重要的方向约定。如果改成 χ₊=χ₋S，所有反射系数及 RHP 跳跃符号都要相应改变。

为了使用完整传播相位，定义同时归一化的基本矩阵

$$
\Psi_\pm=\chi_\pm e^{\mathrm{i}\rho Fz\sigma_3},
\qquad
\Psi_-=\Psi_+S.
$$

S 与 S_sp 的区别：S 是两部分 Lax 对同时演化的散射矩阵，在后面证明的正则性条件下与 z 无关；S_sp 则可随 z 通过对角相位共轭变化。

## 3.7 eq:detjost：为什么散射矩阵行列式为 1？

对基本矩阵 χₓ=Uχ，有 Liouville（Jacobi）行列式公式

$$
\dfrac{\mathrm{d}}{\mathrm{d}x}\det\chi
=\operatorname{tr}U\det\chi.
$$

这里 tr U=0，所以 detχ 与 x 无关；利用两端归一化

$$
\det\chi_\pm=\det E=\gamma(w).
$$

于是

$$
\det S_{\mathrm{sp}}
=\dfrac{\det\chi_-}{\det\chi_+}=1.
$$

同时归一化右乘的对角指数行列式也为 1，所以 detS=1。这个等式是以后计算跳跃矩阵 detJ=1 的原因之一。

## 3.8 eq:wronskian：Cramer 法则逐列推导

把 Ψ₋=Ψ₊S 写成列：

$$
\begin{aligned}
\Psi_{-,1}&=s_{11}\Psi_{+,1}+s_{21}\Psi_{+,2},\\
\Psi_{-,2}&=s_{12}\Psi_{+,1}+s_{22}\Psi_{+,2}.
\end{aligned}
$$

对第一式两边与 Ψ₊,₂ 构成行列式：

$$
\det(\Psi_{-,1},\Psi_{+,2})
=s_{11}\det(\Psi_{+,1},\Psi_{+,2}),
$$

因为包含两列 Ψ₊,₂ 的那项为零。分母 detΨ₊=γ，得到

$$
s_{11}=\dfrac{\det(\Psi_{-,1},\Psi_{+,2})}{\gamma}.
$$

同理，与 Ψ₊,₁ 构成适当的行列式，得到

$$
\begin{aligned}
s_{22}&=\dfrac{\det(\Psi_{+,1},\Psi_{-,2})}{\gamma},\\
s_{21}&=\dfrac{\det(\Psi_{+,1},\Psi_{-,1})}{\gamma},\\
s_{12}&=\dfrac{\det(\Psi_{-,2},\Psi_{+,2})}{\gamma}.
\end{aligned}
$$

最后一个分子列序不能交换，否则会多一个负号。

由于 s₁₁ 的分子是 D₊ 的两列解析函数的行列式，且 γ 在普通点非零，s₁₁ 可以延拓进 D₊；s₂₂ 类似在 D₋ 解析。反射相关的 s₁₂、s₂₁ 一般只在连续谱定义，不应无条件宣称它们在某个开域解析。

## 3.9 eq:conjugation、eq:sconj：为何存在复共轭对称？

定义

$$
K=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
$$

空间 Lax 矩阵的实参数约化意味着，如果 Ψ(x,z,\overline w) 是某个基本解，把它逐项共轭后再做 K 共轭，也产生具有相应归一化的基本解。Jost 解唯一，故

$$
\Psi_\pm(w)=K\overline{\Psi_\pm(\overline w)}K^{-1}.
$$

对散射矩阵也有相同形式。若 S=\begin{pmatrix}s₁₁&s₁₂\\s₂₁&s₂₂\end{pmatrix}，直接相乘可得到

$$
s_{22}(w)=\overline{s_{11}(\overline w)},
\qquad
s_{12}(w)=-\overline{s_{21}(\overline w)}.
$$

这一对称把 D₊ 的透射零点送到 D₋ 的共轭零点。它与换叶对称是两种不同操作，不可以混为一谈。

## 3.10 eq:Esheet、eq:musheet：换叶与共轭不是一回事

从第二章

$$
h(w)=-\dfrac{Q^2}{w},\quad
\lambda(h(w))=\lambda(w),\quad \rho(h(w))=-\rho(w)
$$

出发，对 E 直接计算可得

$$
E(h(w))=E(w)C(w),
\qquad
C(w)=\dfrac{\mathrm{i}w}{Q}\sigma_1.
$$

C σ₃=−σ₃ C，正好把两种背景指数模式交换。于是

$$
\Psi_\pm(h(w))=\Psi_\pm(w)C(w),
$$

并由 Ψ₋=Ψ₊S 导出

$$
S(h(w))=\sigma_1S(w)\sigma_1.
$$

后者相当于交换透射和反射矩阵的列、行位置。请动手算 σ₁Sσ₁，确认它确实使 s₁₁↔s₂₂、s₁₂↔s₂₁。

**与共轭区别**：h 保持 λ 不变但改变 ρ 符号；复共轭通常把 λ、ρ 变成共轭值，改变相应解析域。两者共同生成四点轨道。

## 3.11 eq:reflection：反射系数只是比例

原文采用

$$
r=\dfrac{s_{21}}{s_{11}},
\qquad
\widetilde r=\dfrac{s_{12}}{s_{22}}.
$$

这是散射矩阵两个非对角系数相对于对应透射系数的比值。若 s₁₁ 或 s₂₂ 在连续谱为零，比例就可能没有定义，因此论文额外假定普通连续谱无透射零点。

反射为零时 r=\widetilde r=0，意味着 S 在普通连续谱是对角矩阵，并不意味着 S=I 或没有离散零点。第五章的 S 恰为一个非平凡的对角矩阵。

## 3.12 prop:quartet 与 eq:quartet：四点离散谱的来历

若 ζ∈D₊ 且 s₁₁(ζ)=0，那么由共轭对称得到 s₂₂(\overline ζ)=0；由换叶对称得到 s₂₂(h(ζ))=0；再组合两种对称得到 s₁₁(h(\overline ζ))=0。四点为

$$
\left\{\zeta,\overline\zeta,-\dfrac{Q^2}{\zeta},-\dfrac{Q^2}{\overline\zeta}\right\}.
$$

其中 ζ、h(\overline ζ) 在 D₊，另外两个在 D₋。一般简单零点的四点互不相同。第五章外部纯虚谱点 iW 的 D₊ 配对点是 −iQ²/W。

为什么透射零点对应本征函数？在 s₁₁(ζ)=0 处，

$$
\det(\Psi_{-,1}(\zeta),\Psi_{+,2}(\zeta))=0.
$$

两列线性相关，可以形成同一向量解：它向左、右分别由所选的衰减模态控制，因而与局域本征态相联系。线性相关常数 b_j 是“范数常数”的原始来源。

## 3.13 b_j、c_j 与 C_j 各表示什么？

原文定义

$$
\Psi_{-,1}(\zeta_j)=b_j\Psi_{+,2}(\zeta_j),
\qquad
c_j=\dfrac{b_j}{s_{11}'(\zeta_j)}.
$$

由于 s₁₁ 在 ζ_j 有简单零点，1/s₁₁ 的留数为 1/s'₁₁(ζ_j)。当第一列被除以 s₁₁ 后，留数系数就出现 c_j。去振荡后的列还引入

$$
C_j(x,z)=c_j e^{2\mathrm{i}\theta(\zeta_j;x,z)}.
$$

请区分：c_j 是在正确的同时归一化约定下与 x、z 无关的常数；C_j 由于含 θ，仍然随 x、z 变化。在第五章有限维系统里真正出现的是 C_j，不是把 c_j 直接代进去。

## 3.14 prop:evolution：两种散射矩阵为什么演化方式不同？

背景时间矩阵为 V∞=−FU∞，对角模式传播含 e^{iρFzσ₃}。定义 Ψ±=χ±e^{iρFzσ₃} 之后，若 Ψ± 同时满足两部分 Lax 方程，比较同一组基本解得到 S 不依赖 x、z。

由 χ₋=χ₊S_sp 与 Ψ₋=Ψ₊S，代入 χ=Ψe^{−iρFzσ₃}，整理：

$$
S_{\mathrm{sp}}(w,z)
=e^{\mathrm{i}\rho Fz\sigma_3}S(w)
e^{-\mathrm{i}\rho Fz\sigma_3}.
$$

写 E_z=\operatorname{diag}(e^{iρFz},e^{-iρFz})，则

$$
E_z\begin{pmatrix}s_{11}&s_{12}\\s_{21}&s_{22}\end{pmatrix}E_z^{-1}
=\begin{pmatrix}
s_{11}&e^{2\mathrm{i}\rho Fz}s_{12}\\
e^{-2\mathrm{i}\rho Fz}s_{21}&s_{22}
\end{pmatrix}.
$$

因此对角透射系数不变，而空间归一化反射系数满足原文 eq:timer：

$$
r_{\mathrm{sp}}(w,z)
=e^{-2\mathrm{i}\rho Fz}r_{\mathrm{sp}}(w,0),
\quad
\widetilde r_{\mathrm{sp}}(w,z)
=e^{2\mathrm{i}\rho Fz}\widetilde r_{\mathrm{sp}}(w,0).
$$

对 S 不依赖 z 的证明，用零曲率对 χ_{±,z}−Vχ± 求 x 导数，识别其端态，再利用 Jost 唯一性。这一步要求对 z 微分合法，因此论文额外引入 W^{3,1}、C¹ 等充分正则性。不能把这些演化公式无条件用于所有粗糙数据。

## 3.15 prop:resonance：为什么在共振点 r_sp(a)=0？

这不是从“F 在 a 发散”直接得出的，而是从介质方程严格推导：

$$
R_x=[\widehat U(a),R].
$$

对共振处可逆的空间基本解 χ，利用 χₓ=\widehat U(a)χ：

$$
\begin{aligned}
(\chi^{-1}R\chi)_x
&=-\chi^{-1}\widehat U(a)R\chi
+\chi^{-1}R_x\chi
+\chi^{-1}R\widehat U(a)\chi\\
&=\chi^{-1}\bigl(R_x-[\widehat U(a),R]\bigr)\chi=0.
\end{aligned}
$$

所以 χ⁻¹Rχ 沿 x 为常数。左右端态及 E(a) 的对角化，确定该常矩阵都是

$$
d_a\sigma_3,\qquad
d_a=\dfrac{\eta_0\rho(w_a)}{a}.
$$

利用 χ₋=χ₊S_sp：

$$
S_{\mathrm{sp}}^{-1}d_a\sigma_3S_{\mathrm{sp}}
=d_a\sigma_3.
$$

因为 dₐ≠0，等价于 [σ₃,S_sp]=0。设 S_sp=\begin{pmatrix}a&b\\c&d\end{pmatrix}，则

$$
[\sigma_3,S_{\mathrm{sp}}]
=\begin{pmatrix}0&2b\\-2c&0\end{pmatrix}=0
\ \Longrightarrow\ b=c=0.
$$

所以反射系数在共振点为零。

**所有前提都必须记住**：固定 z，介质方程成立；两端有相容背景；Jost 基在 a 正则可逆；η₀≠0；透射分母非零。这不是一般无条件结论，也不是“r(a)=0 足以保证 RHP 一般可解”。

对于不同介质端态，左右常矩阵可能不同，得到的是相似关系而非与 σ₃ 对易。

## 3.16 prop:highenergy：为什么需要做高能展开？

第四章需要 M(w)=I+M⁽¹⁾/w+O(w⁻²)，第五章要从 M⁽¹⁾ 读取光场 u。光说“由渐近可得”不足以计算系数，因此本章在额外条件 u−Q∈W^{2,1} 下推导高能展开。

原文引入与空间势对应的 Hermitian 矩阵

$$
B(u)=\begin{pmatrix}0&u\\\overline u&0\end{pmatrix},
\qquad A(u)=\sigma_3B(u).
$$

及累计积分

$$
J_\pm(x)=\int_{\pm\infty}^{x}(Q^2-|u(y)|^2)\,\mathrm{d}y.
$$

为什么这个积分是有限的？因为 u−Q 可积而 u 有界，且

$$
Q^2-|u|^2
=(Q-u)\overline u+Q(\overline Q-\overline u)
$$

在 Q 实值时可由 |u−Q| 控制，因此可积。注意 J_± 端点不同，后面差值变成全实轴积分。

论文高能结论 eq:mu-highenergy：

$$
\mu_\pm=I+\dfrac{-\mathrm{i}B(u)+\mathrm{i}J_\pm\sigma_3}{w}
+O(w^{-2}),
\qquad w\to\infty.
$$

这里公式矩阵式需按具有相应解析性的列理解。

## 3.17 为什么要设 T=I−iB/w？

直接从原积分方程展开，背景 E 与非对角指数同时出现，第一阶系数不易提取。做“局部规范”：

$$
T(x,w)=I-\dfrac{\mathrm{i}B(u(x))}{w},
\qquad \mu_\pm=Tn_\pm.
$$

由于 B²=|u|²I，

$$
\begin{aligned}
T\left(I+\dfrac{\mathrm{i}B}{w}\right)
&=I+\dfrac{B^2}{w^2}
=\left(1+\dfrac{|u|^2}{w^2}\right)I,\\
T^{-1}
&=\dfrac{I+\mathrm{i}B/w}{1+|u|^2/w^2}.
\end{aligned}
$$

原文 n-volterra-system 最关键的一步是将 UT+iρTσ₃−Tₓ 化简。使用 A=σ₃B、Bσ₃=−σ₃B、λ+ρ=w 和 ρ−λ=Q²/w，得到

$$
UT+\mathrm{i}\rho T\sigma_3-T_x
=\dfrac{\mathrm{i}}{w}\bigl((Q^2-|u|^2)\sigma_3+B_x\bigr).
$$

所以 n 满足

$$
(n_\pm)_x=-\mathrm{i}\rho[\sigma_3,n_\pm]
+\mathcal Cn_\pm,
$$

其中

$$
\mathcal C=\dfrac{1}{w}T^{-1}
\,\mathrm{i}\bigl((Q^2-|u|^2)\sigma_3+B_x\bigr).
$$

由于可积的扰动乘了 1/w，积分迭代自然得到 n−I=O(1/w)。

## 3.18 第一迭代如何抽取 J_±，非对角项为什么更小？

第一次迭代中对角项没有快速振荡，主项就是

$$
\dfrac{\mathrm{i}}{w}\int_{\pm\infty}^{x}
(Q^2-|u(y)|^2)\,\mathrm{d}y
=\dfrac{\mathrm{i}J_\pm(x)}{w}.
$$

非对角项会出现

$$
\dfrac{1}{w}\int_{-\infty}^{x}
e^{2\mathrm{i}\rho(x-y)}\overline{u_y(y)}\,\mathrm{d}y.
$$

以积分内核做分部积分（不引入机械 u、v 记号）：

$$
\begin{aligned}
&\int_{-\infty}^{x}
e^{2\mathrm{i}\rho(x-y)}\overline{u_y(y)}\,\mathrm{d}y\\
&=-\dfrac{\overline{u_x(x)}}{2\mathrm{i}\rho}
+\dfrac{1}{2\mathrm{i}\rho}
\int_{-\infty}^{x}
e^{2\mathrm{i}\rho(x-y)}\overline{u_{yy}(y)}\,\mathrm{d}y.
\end{aligned}
$$

因为在允许的解析区域指数模≤1、uₓ 有界、uₓₓ 可积，右侧为 O(1/ρ)=O(1/w)。外面再乘 1/w，得到 O(w⁻²)；这就解释第一阶非对角修正实际上已由 T 中的 −iB/w 捕获。

## 3.19 透射高能系数和第二个无穷远

由 s₁₁=det(μ₋,₁,μ₊,₂)/γ，把两列展开，得到

$$
s_{11}(w)=1+\dfrac{\mathrm{i}J_{\mathrm{tot}}}{w}
+O(w^{-2}),
$$

其中

$$
J_{\mathrm{tot}}=
\int_{-\infty}^{+\infty}(Q^2-|u|^2)\,\mathrm{d}x.
$$

理由是 J₋(x)−J₊(x) 等于整个实轴积分，与 x 无关。由于 γ=1+O(w⁻²)，不会改变 1/w 的系数。

换叶 h(w)=−Q²/w 将大 w 送到小 w。把大 w 的 Jost 展开和换叶对称结合，可得 eq:smallmu：

$$
\mu_\pm(w)=-\dfrac{\mathrm{i}Q}{w}\sigma_1+O(1),
\qquad
S(w)=I+O(w),\qquad w\to0.
$$

这说明非零背景 RHP 有**两种归一化**：w→∞ 时趋于 I，w→0 时有指定 1/w 主项。第四章的部分分式式因此必须包含固定的 −iQσ₁/w 项。

## 3.20 自测与阅读关卡

1. 不看正文，自己从 W=E⁻¹χ 推导 Volterra 方程。
2. 判断 e^{2iρ(x−y)} 在 x−y≥0、Imρ>0 下的模。
3. 用 Cramer 法则独立求出 s₁₁ 和 s₁₂，并检查列顺序的负号。
4. 把 S_sp 的对角相位共轭做完整 2×2 乘法，得出两个反射系数的指数。
5. 用两次对称操作找出 ζ 的四点轨道。
6. 写出 r_sp(a)=0 定理的全部条件，特别是左右介质端态和 Jost 正则性。
7. 亲自完成 T⁻¹ 的运算与一次分部积分，说明 O(w⁻²) 从哪里来。

读懂这一章后，第四章的 RHP 矩阵不再是某种“凭空定义”的未知函数，而是我们已经知道解析域的 Jost 列的自然拼接。
