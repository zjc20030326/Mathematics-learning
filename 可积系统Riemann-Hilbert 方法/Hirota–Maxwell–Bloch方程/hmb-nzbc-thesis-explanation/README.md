# HMB 论文逐式推导详细版

本册与相邻 hmb-nzbc-thesis 的当前正文同步重写。保留正文的定义、命题、引理、定理、证明及编号公式，增加不编号的矩阵计算、Volterra 估计、留数、消元、特殊点分析与导数核验。两册采用相同八节结构。

## 阅读与源码

- hmb-explanation.tex：入口。
- chapters/00-reading-guide.tex：记号与正向、逆向条件的说明。
- chapters/01-introduction.tex 至 chapters/08-numerics.tex：逐节扩展。
- chapters/10-formula-index.tex：当前正文公式索引。
- hmb-explanation.pdf：成品。

主要散射对象为空间 Jost 解和随时间演化的散射矩阵。初始常数与演化范数数据区分，共振取值使用空间归一化。一般含反射表示、多轨道条件性公式与指定单轨道的全局结果分别说明。

原先独立的导数补充并入第8节。摘要和最终学位格式暂缓处理。

## 编译

在本目录运行：

    .\build.ps1

脚本通过 ASCII 临时目录完成完整 XeLaTeX 与 BibTeX 编译，再复制 PDF 与 SyncTeX 回项目。

## 验收

checks/check_source.py 比较当前正文和详细版的标签、编号环境、公式、正文保留顺序、定界符与公式索引，并核对两册引用编号。同步基准是本次重写的正文，不要求保留旧九节或旧编号。

checks/verify_calculations.py 与 checks/verify_derivatives.py 检查补充计算及导数。checks/inspect_pdf.py 检查页面几何并生成全页缩略图；链接核对与实际目视检查另记在 checks/acceptance.json。

同步哈希和引用映射在两册完成编译后更新。有限样点检查不证明一般含反射 RHP 或一般多轨道的全局可解性。
