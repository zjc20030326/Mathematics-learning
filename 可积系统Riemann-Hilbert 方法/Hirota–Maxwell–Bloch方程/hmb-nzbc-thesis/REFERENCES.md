# 文献核查记录

核查基准：2026-10-21。六篇主体论文均已在上级“参考论文”目录取得本地全文。本文页码默认是实际PDF从1开始的页码，不将预印本式号伪装为期刊终版式号。本地文件页数、SHA-256见[local-papers.json](evidence/local-papers.json)，文本提取保存在evidence/paper01.txt至paper06.txt，关键Lax页另有PNG截图。

书目信息由本地首页、DOI出版社页面与Crossref登记交叉核对。Crossref原始响应保存在evidence/doi-*.json。网上检索不构成“此前无人研究”的证明。

## 1. Li–Xia–Wei（2021）

Jian Li, Tiecheng Xia, Hanyu Wei. **The N-soliton solutions to the Hirota and Maxwell–Bloch equation via the Riemann–Hilbert approach**. International Journal of Modern Physics B, 35(11), 2150153 (2021).

- DOI：[10.1142/S0217979221501538](https://doi.org/10.1142/S0217979221501538)。
- 本地全文：paper06，9页，期刊版，页脚2150153-1至2150153-9。
- 已核对：第2–3页式(2.1)完整特化Lax对；第3页式(2.6)、(2.8)–(2.9)的零光场归一化。
- 用途：HMB零光场RHP工作的边界与流程比较。
- 限制：其特化符号与一般参数版本存在需要注意的组合差别，本稿不直接照抄为自己的完整Lax对。非零光场谱曲线不能由此直接得出。
- 日期差异：本地首页写Published 19 June 2021，Crossref记录2021-04-30；本稿只使用一致的年份，不虚构统一日级出版日期。

## 2. Zhao–Zhaqilao（2026）

Yijie Zhao, Zhaqilao. **Breathers and rogue waves on the periodic background in the Hirota and the Maxwell-Bloch equations**. Wave Motion, 142, 103696 (2026).

- DOI：[10.1016/j.wavemoti.2025.103696](https://doi.org/10.1016/j.wavemoti.2025.103696)。
- 本地全文：paper01，16页，期刊版；Crossref卷142、文章号103696、出版年2026一致。DOI中的2025不是本稿所用期刊年份。
- 已核对：第2页式(1a)–(1c)、(2)–(3)；第3页式(12)一重/二重光场Darboux表达；第3–4页平面波与周期背景讨论。
- 变量转换：原文传播变量t改为本稿z，原文omega=-a，u改为q；不改变alpha、beta。
- 用途：非零光场Darboux解确属既有工作，比较而不宣称新波形。
- 限制：不将该文视作一般NZBC散射RHP证明。

## 3. Liu–Liu（2026）

Shuai Liu, Yaqing Liu. **Inverse scattering transform and soliton solutions for the coupled Hirota-Maxwell-Bloch system with nonzero background**. Nonlinear Dynamics, 114(12), article 821 (2026).

- DOI：[10.1007/s11071-026-12724-0](https://doi.org/10.1007/s11071-026-12724-0)。
- 本地全文：paper03，16页；出版社页面与首页均写2026-06-14在线发表。
- 已核对：第1页摘要明确说非零背景仅加在eta分量，u、p在无穷远趋零；不是本稿Q>0的边界。
- 用途：准确排除“题名相同即已解决同一光场边界问题”的混淆。
- 没有据此宣称本稿研究世界首次，也没有移用其零光场Jost归一化。

## 4. Biondini–Gabitov–Kovačič–Li（2019）

Gino Biondini, Ildar Gabitov, Gregor Kovačič, Sitai Li. **Inverse scattering transform for two-level systems with nonzero background**. Journal of Mathematical Physics, 60(7), 073510 (2019).

- DOI：[10.1063/1.5084720](https://doi.org/10.1063/1.5084720)；[arXiv:1907.06231](https://arxiv.org/abs/1907.06231)。
- 本地全文：paper04，63页，预印本；首页出现December 7, 2021排版日期及arXiv标识。正式书目仍按2019期刊条目，不将本地排版日期当作期刊年份。
- 已核对：第1页摘要、引言中真正非零光场及介质端态相容性的研究范围。
- 用途：说明非零光场MB逆散射已有重要基础。
- 限制：不直接借用该文的展宽/频率相关模型作为单一实极点HMB定理。本稿未引用其未逐式核实的定理编号。

## 5. Li–He–Porsezian（2013）

Chuanzhong Li, Jingsong He, K. Porsezian. **Rogue waves of the Hirota and the Maxwell-Bloch equations**. Physical Review E, 87(1), 012913 (2013).

- DOI：[10.1103/PhysRevE.87.012913](https://doi.org/10.1103/PhysRevE.87.012913)；[arXiv:1205.1191v3](https://arxiv.org/abs/1205.1191v3)。
- 本地全文：paper05，17页，arXiv v3，不是期刊分页；arXiv页面确认v3于2013-01-25更新。
- 已核对：第1页(1.1)–(1.3)模型；第2页(2.2)–(2.4)完整一般参数Lax矩阵；第3页(3.5)–(3.20)交织关系与变换；第4页(4.2)–(4.4)三场一重变换。
- 精确转换：源文t=x、E=q、omega=-a；传播z与谱lambda不变。
- 用途：本文完整Lax对和已知Darboux定理的主要来源。重构q的符号和介质共振点共轭另行独立验算，零曲率对一般参数已由SymPy检验。
- 原预印本首页拼作Porseizan；书目采用当前出版社与Crossref核对后的Porsezian。

### 2013年出版者更正

Chuanzhong Li, Jingsong He, K. Porsezian. **Publisher's Note: Rogue waves of the Hirota and the Maxwell-Bloch equations [Phys. Rev. E 87, 012913 (2013)]**. Physical Review E, 87(5), 059903 (2013).

- DOI：[10.1103/PhysRevE.87.059903](https://doi.org/10.1103/PhysRevE.87.059903)。
- 书目已由Crossref与APS原论文Corrections栏目核实；独立BibTeX条目LHPnote2013已建立。
- 全文状态：未取得。APS页面访问失败，harvest全文接口返回401。不得据此声称已核对更正全文或其全部数学影响。
- 辅助索引PubMed注明作者姓名更正，但本稿没有把辅助索引当作已阅读出版者更正正文的替代品。

## 6. Zhang–Chen–Yan（2020）

Guoqiang Zhang, Shuyan Chen, Zhenya Yan. **Focusing and defocusing Hirota equations with non-zero boundary conditions: Inverse scattering transforms and soliton solutions**. Communications in Nonlinear Science and Numerical Simulation, 80, 104927 (2020).

- DOI：[10.1016/j.cnsns.2019.104927](https://doi.org/10.1016/j.cnsns.2019.104927)。
- 本地全文：paper02，22页，期刊版；卷80与文章号104927核实。
- 已核对：第3页式(10)–(13)单值化域、Jost归一化与积分式；第4页命题1与式(14)–(20)；第5页式(21)–(31)的四点轨道、留数和矩阵拼接；第6页命题2的RHP表述。迹公式所在文字已定位，但本稿无反射迹公式重新在明确拼接假设下推导。
- 注意：该文空间Lax及散射方向与本稿约定不同，不能机械抄写解析列与跳跃正负号。
- 用途：普通Hirota NZBC结构参照，不能充当HMB共振点与介质重构定理。

## 7. 可选物理背景

K. Porsezian, K. Nakkeeran. **Optical Soliton Propagation in an Erbium Doped Nonlinear Light Guide with Higher Order Dispersion**. Physical Review Letters, 74(15), 2941–2944 (1995).

- DOI：[10.1103/PhysRevLett.74.2941](https://doi.org/10.1103/PhysRevLett.74.2941)。
- 本轮核实：作者、题名、期刊、卷期页码及年份的Crossref元数据。
- 全文未取得；仅作高阶色散掺铒光纤模型的历史背景引用，不引用具体式号、证明或定理。

## 8. 范恩贵本地二校稿

范恩贵：《可积系统、正交多项式和随机矩阵——Riemann-Hilbert 方法》。

- 本轮读取的文件：上两级目录中的“可积系统正交多项式和随机矩阵Riemann-Hilbert 方法二校20220303.pdf”，463页。
- 版权页含占位ISBN与CIP号码，未将2021.10等试排信息当作正式版本证据；BibTeX使用misc条目并说明本地二校稿性质。
- 目录核对：第5章从书内69页开始，第6章从88页开始。
- 实际读取：PDF第84页（书内70页）式(5.1.5)–(5.1.10)；PDF第103页（书内89页）的单值化解释。摘取文本保存在evidence/fan-body-p84.txt及fan-body-p103.txt。
- 用途：组织“直接散射—RHP—解重构”的讲述路径。没有给未核对版本编造出版社年份、ISBN或定理编号。

## 9. 模板与检索边界

指定模板：[yzbrlan/fudan-thesis-latex-template](https://github.com/yzbrlan/fudan-thesis-latex-template)，实际复制graduation-thesis目录，提交d4a4a1e468c91dd3ef15dd0d021b56263f936988。保留原.cls、.def、.sty、logo、bst及其他依赖。项目内兼容修改见BUILD.md及原类文件备份。

本轮检索是围绕任务指定文献的来源核验与交叉比较，未完成穷尽式系统文献综述。因此正文不作“此前无人研究”的排他性创新声明。
