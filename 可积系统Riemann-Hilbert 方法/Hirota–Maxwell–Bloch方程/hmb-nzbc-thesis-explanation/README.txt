HMB 非零背景 Riemann–Hilbert 方法论文全文与逐式精读（LaTeX）

主文件：hmb-explanation.tex
解析章节：chapters/00-preliminaries.tex 至 chapters/07-formula-index.tex
原文对应：../hmb-nzbc-thesis/section/01-introduction.tex 至 06-conclusion.tex

这份解析书籍在同一 PDF 中先收录对应原论文六章正文，再按同一章序讲解背景、推导、矩阵代数、留数条件、逆重构及共振正则性。第五章采用最新版直接 RHP 求解，Darboux 只用于结果对照。

编译（在本文件夹中）：
  latexmk -xelatex -interaction=nonstopmode -halt-on-error hmb-explanation.tex

XeLaTeX + ctexbook，采用 TeX Live Fandol 字体。
原论文所有源文件保持不变。若论文修改，请重新编译并核对公式标签与讲解。
