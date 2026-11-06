# HMB 论文逐式推导详细版

本工程以相邻 `hmb-nzbc-thesis` 目录中的当前论文源码为底稿重新编写. 保留原文九部分结构、原有符号、原公式编号、定理和证明, 在相应位置插入不编号的计算和说明. 原论文已有的规范变换、谱单值化及辅助量按原样解释, 没有另外设置坐标或谱参数代换来缩短推导.

入口文件为 [hmb-explanation.tex](hmb-explanation.tex). 生成的阅读文件为 [hmb-explanation.pdf](hmb-explanation.pdf). `chapters/00-reading-guide.tex` 是阅读约定, `chapters/01` 至 `09` 分别对应原论文九部分. `09-derivative-details.tex` 在第9部分后继续展开三阶导数、介质求导和端态误差. `10-formula-index.tex` 提供127条原公式的跳转和本册页码. `figures/` 保留原论文五幅图, `citelist.bib` 保留原参考文献. 编译模板沿用原项目的 `fduthesis.cls`、`fduthesis.def` 和 `gbt7714-unsrt.bst`, 没有新增 LaTeX 宏包依赖.

详细版的计算重点包括:

- 平面波逐条代入、Bloch 求导抵消, 零曲率每次幂及各独立矩阵元, 规范和背景矩阵的全部相乘.
- 谱曲面、单值化反解、连续谱与解析域、背景特征向量与两共振根.
- 四类 Jost 核的指数方向、Neumann 阶乘估计、唯一性和分支点极限.
- 散射 Wronskian、时间演化、共振条件、两种约化、离散谱和留数的链式导数.
- 两处渐近的系数比较、RHP 跳跃消元、Cauchy 表示、迹公式及同相约束.
- 一轨道有限系统、Schur 正性、显式矩阵的留数与约化、三场矩阵乘法、端态和完整 Lax 对恢复.
- 数值参数、传播速度、内部周期、原点三场、振幅和衰减, 以及残差表的形成方法.

在含反射部分, 表示式仍以满足原文正则性要求的 RHP 解存在为前提. 对指定纯虚一轨道的全局光滑解, 则保留其独立存在唯一性与逆验证证明. 数值残差不替代解析证明.

## 编译

在本目录运行:

```powershell
.\build.ps1
```

需要原模板使用的 XeLaTeX、BibTeX 和 Windows 字体. 默认辅助文件写入系统临时目录 `hmb-explanation-20261104-build`, 成品 PDF 和 SyncTeX 写回本工程. 也可通过 `-OutputDirectory` 指定独立的辅助文件目录. 编译不会覆盖相邻原论文的源文件或 PDF.

## 核验

`checks/source-hashes.json` 与 `checks/source-main-and-pdf-hashes.json` 记录开工时原章节、入口及 PDF 的 SHA-256. `checks/source-reference.json` 从开工时原论文的辅助文件提取编号, 作为原编号对照. 它们用于核对本次工作没有改动原论文, 并在源文档后续更新时提示需要重新同步.

```powershell
python .\checks\check_source.py "$env:TEMP\hmb-explanation-20261104-build\hmb-explanation.aux"
python .\checks\verify_calculations.py
python .\checks\verify_derivatives.py
```

第一个程序检查原文件哈希、原公式显示块、标签和编号、原文顺序、环境与定界符, 报告写入 `checks/source-coverage.json`. 第二个程序从相邻原论文读取其符号及高精度验证程序, 将输出重定向至本工程 `calculations/`, 不向原论文写报告. 它需要现有 Python 环境中的 SymPy 和 mpmath. 独立计算覆盖零曲率、背景谱恒等式、一轨道消元、矩阵解、三条原方程、两条 Lax 方程、两共振重构和边界, 并核对本详细版新增的原点数值. 第三个程序用 SymPy 精确验证本册新增的导数和端态分量恒等式, 逐项结果在 `calculations/derivative-checks.json`.

应在已安装相应依赖的 Python 环境中执行数学核验. 本次实际使用的解释器是 `D:\Software\Python\anaconda3\python.exe`; 默认 `python` 若指向另一个环境, 可以在 PowerShell 中显式调用该解释器. `checks/inspect_pdf.py` 还需要 PyMuPDF 和 Pillow, 用于渲染逐页联系表并检查 PDF 边界, 不参与文档编译.

最终编译、编号和 PDF 版面验收摘要见 [checks/acceptance.md](checks/acceptance.md). 所有检查都有具体范围: 公式覆盖与编译不等于自动证明每个新增解释正确, 数值结果也只对应报告列出的参数与样点.
