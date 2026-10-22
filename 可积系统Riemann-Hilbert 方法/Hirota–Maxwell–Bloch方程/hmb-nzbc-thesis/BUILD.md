# 构建与计算复现

研究日期基准为2026-10-21，工程验收更新于2026-10-22。编译、数学计算及排版检查是三项不同的验证；PDF能够生成不表示一般逆散射定理已经证明。

## 1. 工程位置与入口

工程绝对路径：

    D:\Desktop\数学学习\可积系统Riemann-Hilbert 方法\Hirota–Maxwell–Bloch方程\hmb-nzbc-thesis

- 主源码：[fduthesis.tex](fduthesis.tex)。
- 编译结果：[fduthesis.pdf](fduthesis.pdf)。
- 自动构建：[build.ps1](build.ps1)。
- 正文目录：section/00-frontmatter.tex至section/07-appendix.tex；主文件仅载入这八个新正文文件。
- 文献数据库：[citelist.bib](citelist.bib)。
- 数学限制：[RESEARCH_STATUS.md](RESEARCH_STATUS.md)。
- 文献来源：[REFERENCES.md](REFERENCES.md)。

模板原有的其他section文件和图片作为原始模板资源保留，不参与当前论文正文编译。作者、导师等信息均为待填写。

## 2. 已验证环境与依赖

本机Windows环境实际使用：

| 组件 | 已验证版本或配置 | 用途 |
| --- | --- | --- |
| XeTeX | 0.999998，TeX Live 2026 | 中文、系统字体和LaTeX编译 |
| BibTeX | 0.99e，TeX Live 2026 | 模板指定文献后端 |
| 文档类 | 本地fduthesis，type=master，oneside | 复旦硕士论文版式 |
| CJK字体配置 | cjk-font=windows | 使用已安装的Windows中文字体 |
| Python | 3.13.7，Anaconda | 数学验证和绘图 |
| SymPy | 1.14.0 | 零曲率与精确PDE代回 |
| mpmath | 1.3.0 | 60位精度微分与残差 |
| NumPy | 2.4.6 | 绘图网格计算 |
| Matplotlib | 3.10.8 | 导出三场PDF/PNG图像 |
| PyMuPDF | 1.28.2 | PDF文本和页面检查 |

当前Python可执行文件为 D:\Software\Python\anaconda3\python.exe。图像已经随工程提供，因此仅编译论文不需要Python。

主文件显式载入booktabs、longtable、graphicx、array，用于三线表、长符号表、图像及表格列格式；数学环境和字体依赖由原模板提供。新增宏ii、dd、conj、Res、diag、tr等分别统一虚数单位、微分、共轭及矩阵运算记号；status仅用于明确标识数学状态。

## 3. 一键编译

在PowerShell中运行：

    Set-Location -LiteralPath 'D:\Desktop\数学学习\可积系统Riemann-Hilbert 方法\Hirota–Maxwell–Bloch方程\hmb-nzbc-thesis'
    .\build.ps1

也可给出可写的绝对ASCII临时目录：

    .\build.ps1 -OutputDirectory 'C:\Temp\hmb-thesis-build'

脚本执行XeLaTeX → BibTeX → XeLaTeX → XeLaTeX，逐步检查退出码，最后复制PDF与日志回工程。默认中间文件目录是当前用户TEMP下的hmb-nzbc-thesis-build。每次运行会更新当前工程PDF和构建记录；不会修改其他论文。

主要命令相当于：

    xelatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=<absolute-build-directory> fduthesis.tex

然后将citelist.bib及本地bst复制至构建目录，在该目录执行：

    bibtex fduthesis

最后回到源码目录执行两次同样的XeLaTeX命令。不要将BibTeX替换为biber，也不要在未检查当前目录的情况下对绝对aux路径运行BibTeX。

## 4. 实际排除的两个构建问题

1. 原模板加载ctexbook类，却使用包版本检查查询ctex，导致当前TeX环境误报版本过旧。本工程仅把本地fduthesis.cls内这一次检查改为ctexbook类版本检查，保留原最低日期条件。原类文件完整保存在[evidence/template-fduthesis.cls](evidence/template-fduthesis.cls)，修改处有注释；没有修改TeX安装目录或全局配置。
2. BibTeX对绝对路径aux写出bbl时被openout_any=p限制。解决方式是进入构建目录运行BibTeX，并把该次构建必需的bib与bst放在那里；没有放宽TeX文件访问策略。

指定模板来源为[yzbrlan/fudan-thesis-latex-template](https://github.com/yzbrlan/fudan-thesis-latex-template)，采用graduation-thesis子目录，提交d4a4a1e468c91dd3ef15dd0d021b56263f936988。原入口、文献表、README及类文件均保留备份。

## 5. 复现数学计算与图像

在工程目录中依次执行以下命令；若python不在PATH中，将其替换为上述本机Python绝对路径：

    python scripts/verify_symbolic.py
    python scripts/verify_exact_example.py
    python scripts/verify_solitons.py
    python scripts/verify_rhp.py
    python scripts/plot_fields.py
    python scripts/audit_sources.py

对应输出：

| 程序 | 证据 | 判定范围 |
| --- | --- | --- |
| verify_symbolic | [symbolic-results.json](evidence/symbolic-results.json) | 一般参数的矩阵代数恒等式 |
| verify_exact_example | [exact-example-results.json](evidence/exact-example-results.json) | 一个有理参数例的三条PDE与Bloch偏差精确为零 |
| verify_solitons | [soliton-results.json](evidence/soliton-results.json) | 五个时空点的一重、二重高精度残差及两个端点检查 |
| verify_rhp | [rhp-results.json](evidence/rhp-results.json) | 一个完整谱轨道的留数、谱对称、行列式、空间Lax和介质重构 |
| plot_fields | figures内的PDF/PNG | 251×401网格上的三场图像 |
| audit_sources | [source-audit.json](evidence/source-audit.json) | 作者新增TeX的定界符、环境、标签、引用及报告链接检查 |

高精度脚本用十进制字符串构造参数，避免先经过双精度舍入。数值误差是绝对误差，不声称给出一般定理。绘图用双精度；其Bloch误差约为1e-15，与高精度证明辅助计算分开。

论文原始参数、采样点、各方程最大残差均记录在JSON中。复现同一依赖版本时最后几位浮点误差仍可能受平台影响；检查阈值与正文中的舍入上界用途不同。

## 6. 构建证据与检查范围

最终PDF为40页。XeLaTeX与BibTeX均成功退出；最终日志中没有Overfull、缺字、未定义引用或BibTeX警告。保留1条Underfull排版提示、1条xeCJK初始化提示及2条封面换行的hyperref元数据提示。

最终构建记录位于[evidence/build/fduthesis.log](evidence/build/fduthesis.log)及[evidence/build/fduthesis.blg](evidence/build/fduthesis.blg)。源码检查脚本是辅助检查，不能替代TeX解析或数学审稿。

最终验收数字、PDF页数、字体嵌入和原始参考文献哈希核对见[evidence/delivery-audit.json](evidence/delivery-audit.json)。页面检查图保存在evidence/preview目录，覆盖封面、摘要、目录、长矩阵、RHP、残差表、图像及参考文献。

保留的两类非致命提示分别是初始化时的xeCJK字体族提示，以及封面手动换行被hyperref从PDF元数据中移除。最终页面的中文字体与正文数学已另行渲染核验；此类提示不等于缺字。排版检查还记录Underfull提示，不将其隐瞒为“完全零警告”。正文嵌入字体文件，Matplotlib图像中的Type3字体以字形程序嵌入。

若在其他操作系统上编译，需要修改cjk-font为该平台实际可用的字体配置；不要通过指定不存在的字体解决错误。本文的已验证编译结果限于上述Windows环境。

## 7. 源码归档

父目录中的hmb-nzbc-thesis-source.zip是可独立编译的源码归档，含模板依赖、章节、文献库、图像、计算脚本、状态报告及精选核验记录。编译PDF单独交付；原始参考论文及其全文抽取、页面预览缓存不打包。

运行python scripts/package_project.py可重新生成归档。压缩包内SOURCE_MANIFEST.json列出各文件SHA-256；父目录hmb-nzbc-thesis-package.json记录归档哈希及完整性检查。打包之前应先完成本节前述构建和核验。

页面检查记录见[evidence/VISUAL_REVIEW.md](evidence/VISUAL_REVIEW.md)。辅助脚本audit_delivery.py还会比对本机相邻参考论文目录的原始哈希；迁移至只有源码的机器时，缺少该参考目录不影响LaTeX编译，但无法复现这一项原文件保全检查。
