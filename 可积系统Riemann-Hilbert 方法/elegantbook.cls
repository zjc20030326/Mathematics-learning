% -----------------------------------------------------------------------------
% ElegantBook 文档类（中文注释增强版）
% Copyright 2013-2022 ElegantLaTeX (elegantlatex2e@gmail.com)
% Copyright 2026-     ElegantLaTeX (maintained by ranwang.osbert@outlook.com)
% 原项目许可证：LaTeX Project Public License (LPPL) 1.3c 或其后续版本。
% 原项目地址：https://github.com/ElegantLaTeX/ElegantBook
%
% 本文件是在 ElegantBook v4.7 基础上为“数学学习笔记”做的本地定制。
% 代码逻辑尽量保持 ElegantBook 原样；主要新增内容用“【本地修改】”标出。
%
% ============================ 阅读导航 ============================
% 1. 文档类选项：颜色、语言、定理样式、设备尺寸、数学字体、水印等
% 2. 基础类与超链接：加载 book 类，配置 hyperref
% 3. 页面尺寸：geometry，支持 normal / pad
% 4. 字体：西文字体、中文字体、数学字体
% 5. 颜色主题：blue / green / cyan / gray / black
% 6. 标题页信息：title、subtitle、author、institute、version、logo、cover
% 7. 列表与图表标题：enumerate/itemize、caption
% 8. 数学宏包与水印：【本地修改】增加数学工具及 watermark 开关
% 9. 多语言名称：目录、图表、定理、证明、附录等名称
% 10. 章节标题：chapter / section / subsection / subsubsection 样式
% 11. 定理环境：fancy / simple 两套模式
% 12. 例题、练习、问题、证明、remark 等辅助环境
% 13. 页眉页脚：fancyhdr 默认样式（可在 preamble.tex 中覆盖）
% 14. 封面：\maketitle 的实现
% 15. 代码高亮、版本历史、目录、附录目录修正
% ================================================================
%
% ====================== 当前项目只需重点看这些 ======================
% 你的第五章使用：
%   lang=cn, scheme=chinese, mode=simple, device=normal,
%   color=blue, math=cm, thmcnt=section, 12pt, watermark
%
% 因此平时真正可能需要查看/修改的主要区域是：
%   ① 文档类选项与水印开关；
%   ② normal 页面尺寸；
%   ③ 中文字体与 math=cm；
%   ④ blue 颜色主题；
%   ⑤ simple 定理环境；
%   ⑥ 页眉页脚默认值（你的个人样式在 preamble.tex 覆盖）；
%   ⑦ 目录与附录。
%
% 其他语言（it/fr/nl/hu/de/mn/pt/jp）、fancy 定理框、mtpro2/newtx
% 等代码是 ElegantBook 的通用能力。当前数学组会笔记基本不用改，可以略读。
% ====================================================================
% -----------------------------------------------------------------------------
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{elegantbook}[2026/05/01 v4.7 ElegantBook document class]
% 【本地修改】面向数学学习笔记的通用扩展：
% - 增加数学宏包：mathtools、extarrows、arydshln、tikz-cd、physics2/ab
% - 增加 needspace，用于避免标题或列表项被分页拆开
% - 目录深度与章节编号深度默认设为 3（到 subsubsection）
% - 增加 watermark / nowatermark 文档类选项

%%%
% ================= 1. 文档类选项系统 =================
% kvoptions 与 etoolbox 用于解析 \documentclass[...] 中的选项。
% 常用选项包括：color、lang、mode、device、math、thmcnt、watermark 等。
\RequirePackage{kvoptions}
\RequirePackage{etoolbox}
\SetupKeyvalOptions{family=ELEGANT, prefix=ELEGANT@, setkeys=\kvsetkeys}
\newcommand{\ekv}[1]{\kvsetkeys{ELEGANT}{#1}}
% ================= 文档类选项 =================
% 主要选项及默认值：
%   color=blue       主题色，可选 green/cyan/blue/gray/black
%   lang=en          语言；当前项目用 lang=cn
%   result=answer    是否显示 proof/solution；noanswer 可隐藏答案
%   mode=fancy       定理样式；当前项目用 mode=simple
%   device=normal    页面尺寸；normal=A4，pad=平板尺寸
%   math=cm          数学字体；也支持 newtx、mtpro2
%   marginpar=false  是否预留宽旁注区域
%   toc=onecol       目录单栏；twocol 为双栏目录
%   scheme           中文章节编号方案；当前项目用 scheme=chinese
%   thmcnt=chapter   定理编号层级；当前项目用 thmcnt=section
%   bibend=biber     参考文献后端，可切换 bibtex
%   watermark        【本地新增】开启水印；默认 nowatermark
\DeclareStringOption[blue]{color}
\DeclareStringOption[en]{lang}
\DeclareStringOption[answer]{result}
\DeclareStringOption[fancy]{mode}
\DeclareStringOption[normal]{device}
\DeclareStringOption[cm]{math}
\DeclareStringOption[marginfalse]{marginpar}
\DeclareStringOption[onecol]{toc}
\DeclareStringOption{scheme}

% ----- 【本地修改】水印开关 -----
% 默认不显示水印；文档类选项 `watermark` 用于开启水印，
% `nowatermark` 用于显式关闭水印。
\newif\ifELEGANT@watermark
\ELEGANT@watermarkfalse
\DeclareVoidOption{watermark}{\ELEGANT@watermarktrue}
\DeclareVoidOption{nowatermark}{\ELEGANT@watermarkfalse}

% ----- 兼容旧版写法 -----
\DeclareVoidOption{green}{\ekv{color=green}}
\DeclareVoidOption{cyan}{\ekv{color=cyan}}
\DeclareVoidOption{blue}{\ekv{color=blue}}
\DeclareVoidOption{gray}{\ekv{color=gray}}
\DeclareVoidOption{black}{\ekv{color=black}}
\DeclareVoidOption{nocolor}{\ekv{color=none}}

\DeclareVoidOption{chinese}{\ekv{scheme=chinese}}

\DeclareStringOption[ctexfont]{chinesefont}
\DeclareVoidOption{ctexfont}{\ekv{chinesefont=ctexfont}}
\DeclareVoidOption{founder}{\ekv{chinesefont=founder}}
\DeclareVoidOption{nofont}{\ekv{chinesefont=nofont}}

\DeclareVoidOption{en}{\ekv{lang=en}}
\DeclareVoidOption{cn}{\ekv{lang=cn}}
\DeclareVoidOption{it}{\ekv{lang=it}}
\DeclareVoidOption{fr}{\ekv{lang=fr}}
\DeclareVoidOption{nl}{\ekv{lang=nl}}
\DeclareVoidOption{hu}{\ekv{lang=hu}}
\DeclareVoidOption{de}{\ekv{lang=de}}
\DeclareVoidOption{mn}{\ekv{lang=mn}}
\DeclareVoidOption{pt}{\ekv{lang=pt}}
\DeclareVoidOption{jp}{\ekv{lang=jp}}

\DeclareVoidOption{fancy}{\ekv{mode=fancy}}
\DeclareVoidOption{simple}{\ekv{mode=simple}}
\DeclareBoolOption[false]{usesamecnt}

\DeclareVoidOption{answer}{\ekv{result=answer}}
\DeclareVoidOption{noanswer}{\ekv{result=noanswer}}

\DeclareVoidOption{normal}{\ekv{device=normal}}
\DeclareVoidOption{pad}{\ekv{device=pad}}

\DeclareStringOption[numeric-comp]{citestyle}
\DeclareStringOption[numeric]{bibstyle}

\DeclareVoidOption{newtx}{\ekv{math=newtx}}
\DeclareVoidOption{mtpro2}{\ekv{math=mtpro2}}
\DeclareVoidOption{cm}{\ekv{math=cm}}

\DeclareVoidOption{margintrue}{\ekv{marginpar=margintrue}}
\DeclareVoidOption{marginfalse}{\ekv{marginpar=marginfalse}}


\DeclareVoidOption{onecol}{\ekv{toc=onecol}}
\DeclareVoidOption{twocol}{\ekv{toc=twocol}}

\DeclareStringOption[chapter]{thmcnt}
\DeclareVoidOption{chapter}{\ekv{thmcnt=chapter}}
\DeclareVoidOption{section}{\ekv{thmcnt=section}}

\DeclareStringOption[biber]{bibend}
\DeclareVoidOption{biber}{\ekv{bibend=biber}}
\DeclareVoidOption{bibtex}{\ekv{bibend=bibtex}}


% ----- 其他基础命令 -----
\newcommand\mailto[1]{\href{mailto:#1}{\nolinkurl{#1}}}

% ================= 标题样式选项 =================
\DeclareStringOption[hang]{titlestyle}[hang]
% ----- 兼容旧版写法 -----
\DeclareVoidOption{hang}{\ekv{titlestyle=hang}}
\DeclareVoidOption{display}{\ekv{titlestyle=display}}
% ----- 未识别选项的默认处理 -----
\DeclareDefaultOption{\PassOptionsToClass{\CurrentOption}{book}}

\ProcessKeyvalOptions*\relax
% \ProcessOptions*\relax
\LoadClass[a4paper,oneside]{book}

% ================= 2. 基础排版辅助 =================
% setspace：控制行距；csquotes：配合 biblatex 处理引号与引文格式。
\RequirePackage{setspace}

\RequirePackage{csquotes}


% ================= 3. PDF 超链接与书签 =================
% hyperref 负责目录跳转、公式/章节引用链接、URL、PDF 书签等。
% 这里保留 ElegantBook 的通用默认设置；个人链接颜色不写入 class。
\RequirePackage{hyperref}
\hypersetup{
  breaklinks,
  unicode,
  linktoc=all,
  bookmarksnumbered=true,
  bookmarksopen=true,
  pdfkeywords={ElegantBook},
  colorlinks,
  linkcolor=winered,
  citecolor=winered,
  urlcolor=winered,
  plainpages=false,
  pdfstartview=FitH,
  pdfborder={0 0 0},
  linktocpage
}
\newif\ifELEGANT@Hy@later
\@ifpackagelater{hyperref}{2022/05/16}
{\ELEGANT@Hy@latertrue}
{\ELEGANT@Hy@laterfalse}

%% ================= 页面尺寸与版心设置 =================
% ================= 4. 页面尺寸与版心 =================
% device=normal：A4 页面；device=pad：平板阅读尺寸。
% 这里保持 ElegantBook 原生 geometry 参数，不做个人化覆盖。
\RequirePackage{geometry}
\ifdefstring{\ELEGANT@device}{normal}{
  \geometry{
    a4paper,
    top=25.4mm, bottom=25.4mm,
    left=20mm, right=20mm,
    headheight=2.17cm,
    headsep=4mm,
    footskip=12mm
  }
  \ifdefstring{\ELEGANT@marginpar}{margintrue}{
    \geometry{
      marginparwidth=5cm, marginparsep=5mm,
      left=2cm,right=7cm}}{\relax}}{
  \relax}

\ifdefstring{\ELEGANT@device}{pad}{
  \geometry{
    paperwidth=7.5in,
    paperheight=10in,
    margin=16mm,
    headheight=2.17cm,
    footskip=4mm
  }}{\relax}


\RequirePackage{indentfirst,comment}
% ================= 字体设置 =================
\ifdefstring{\ELEGANT@math}{mtpro2}{
  \let\Bbbk\relax
  \RequirePackage[lite]{mtpro2}
}{\relax}

% 【本地修改】目录显示到 subsubsection；章节编号也统一到 subsubsection。
\setcounter{tocdepth}{3}
\renewcommand{\baselinestretch}{1.2}

% ================= 5. 西文字体与中文字体 =================
% pdfLaTeX 与 XeLaTeX/LuaLaTeX 使用不同的字体方案。
% lang=cn 时由 ctex/xeCJK 接管中文字体；当前建议用 XeLaTeX 编译中文。
\PassOptionsToPackage{no-math}{fontspec}
\PassOptionsToPackage{quiet}{fontspec}
\RequirePackage{iftex}

\ifpdftex
  \RequirePackage{newtxtext}
  \RequirePackage[scaled=.90]{helvet}
\else
  \RequirePackage[no-math]{fontspec}
  \setmainfont{TeXGyreTermesX}[
    UprightFont = *-Regular ,
    BoldFont = *-Bold ,
    ItalicFont = *-Italic ,
    BoldItalicFont = *-BoldItalic ,
    Extension = .otf ,
    Scale = 1.0]

  \setsansfont{texgyreheros}[
    UprightFont = *-regular ,
    BoldFont = *-bold ,
    ItalicFont = *-italic ,
    BoldItalicFont = *-bolditalic ,
    Extension = .otf ,
    Scale = 0.9]
\fi

\ifdefstring{\ELEGANT@lang}{cn}{
  \ifdefstring{\ELEGANT@chinesefont}{founder}{
    \RequirePackage[UTF8, scheme=plain, fontset=none]{ctex}
    \setCJKmainfont[BoldFont={FZHei-B01},ItalicFont={FZKai-Z03}]{FZShuSong-Z01}
    \setCJKsansfont[BoldFont={FZHei-B01}]{FZKai-Z03}
    \setCJKmonofont[BoldFont={FZHei-B01}]{FZFangSong-Z02}
    \setCJKfamilyfont{zhsong}{FZShuSong-Z01}
    \setCJKfamilyfont{zhhei}{FZHei-B01}
    \setCJKfamilyfont{zhkai}[BoldFont={FZHei-B01}]{FZKai-Z03}
    \setCJKfamilyfont{zhfs}[BoldFont={FZHei-B01}]{FZFangSong-Z02}
    \newcommand*{\songti}{\CJKfamily{zhsong}}
    \newcommand*{\heiti}{\CJKfamily{zhhei}}
    \newcommand*{\kaishu}{\CJKfamily{zhkai}}
    \newcommand*{\fangsong}{\CJKfamily{zhfs}}}{\relax}

  \ifdefstring{\ELEGANT@chinesefont}{nofont}{
    \RequirePackage[UTF8,scheme=plain,fontset=none]{ctex}
    \xeCJKsetup{AutoFakeBold=true}
  }{\relax}

  \ifdefstring{\ELEGANT@chinesefont}{ctexfont}{
    \RequirePackage[UTF8,scheme=plain]{ctex}
    \xeCJKsetup{AutoFakeBold=true}
  }{\relax}

  \AfterEndPreamble{
    \setlength\parindent{2\ccwd}}
}{\relax}

\ifcsname heiti\endcsname
  \newcommand{\cbfseries}{\heiti}
\else
  \newcommand{\cbfseries}{\bfseries}
\fi


\ifcsname kaishu\endcsname
  \newcommand{\citshape}{\kaishu}
\else
  \newcommand{\citshape}{\itshape}
\fi
\ifcsname kaishu\endcsname
  \newcommand{\cnormal}{\kaishu}
\else
  \newcommand{\cnormal}{\normalfont}
\fi

\ifcsname fangsong\endcsname
  \newcommand{\cfs}{\fangsong}
\else
  \newcommand{\cfs}{\normalfont}
\fi

% ================= 6. 数学字体选项 =================
% math=cm / newtx / mtpro2 用于选择数学字体方案，默认 math=cm。
\RequirePackage{anyfontsize}
\ifdefstring{\ELEGANT@math}{newtx}{
  \let\oldencodingdefault\encodingdefault
  \let\oldrmdefault\rmdefault
  \let\oldsfdefault\sfdefault
  \let\oldttdefault\ttdefault
  \def\encodingdefault{T1}
  \renewcommand{\rmdefault}{ntxtlf}
  \renewcommand{\sfdefault}{qhv}
  \renewcommand{\ttdefault}{ntxtt}
  \RequirePackage{newtxmath}
  \let\encodingdefault\oldencodingdefault
  \let\rmdefault\oldrmdefault
  \let\sfdefault\oldsfdefault
  \let\ttdefault\oldttdefault
  \let\Bbbk\relax
  \RequirePackage{esint}
  %%% 如需 yhmath 数学字体，可取消下面几行注释
  % \let\oldwidering\widering
  % \let\widering\undefined
  % \RequirePackage{yhmath}
  % \let\widering\oldwidering

  %%% 如需 esvect 向量箭头，可取消下一行注释
  % \RequirePackage{esvect}

  \DeclareSymbolFont{CMlargesymbols}{OMX}{cmex}{m}{n}
  \let\sumop\relax\let\prodop\relax
  \DeclareMathSymbol{\sumop}{\mathop}{CMlargesymbols}{"50}
  \DeclareMathSymbol{\prodop}{\mathop}{CMlargesymbols}{"51}
}{\relax}


% ================= 颜色主题 =================
%% 章节以及页脚图形

\RequirePackage[table]{xcolor}
\ifdefstring{\ELEGANT@color}{green}{
  \definecolor{structurecolor}{RGB}{0,120,2}%
  \definecolor{main}{RGB}{0,120,2}%
  \definecolor{second}{RGB}{230,90,7}%
  \definecolor{third}{RGB}{0,160,152}%
}{\relax}
\ifdefstring{\ELEGANT@color}{cyan}{
  \definecolor{structurecolor}{RGB}{31,186,190}%
  \definecolor{main}{RGB}{59,180,5}%
  \definecolor{second}{RGB}{175,153,8}%
  \definecolor{third}{RGB}{244,105,102}%
}{\relax}
\ifdefstring{\ELEGANT@color}{blue}{
  \definecolor{structurecolor}{RGB}{60,113,183}
  \definecolor{main}{RGB}{0,166,82}%
  \definecolor{second}{RGB}{255,134,24}%
  \definecolor{third}{RGB}{0,174,247}%
}{\relax}
\ifdefstring{\ELEGANT@color}{gray}{
  \definecolor{structurecolor}{RGB}{150,150,150}
  \definecolor{main}{RGB}{150,150,150}%
  \definecolor{second}{RGB}{150,150,150}%
  \definecolor{third}{RGB}{150,150,150}%
}{\relax}
\ifdefstring{\ELEGANT@color}{black}{
  \definecolor{structurecolor}{RGB}{0,0,0}
  \definecolor{main}{RGB}{0,0,0}%
  \definecolor{second}{RGB}{0,0,0}%
  \definecolor{third}{RGB}{0,0,0}%
}{\relax}

% 绿色主题备用颜色
\definecolor{structure1}{RGB}{0,120,2}%
\definecolor{main1}{RGB}{0,120,2}%
\definecolor{second1}{RGB}{230,90,7}%
\definecolor{third1}{RGB}{0,160,152}%
% 青色主题备用颜色
\definecolor{structure2}{RGB}{31,186,190}%
\definecolor{main2}{RGB}{59,180,5}%
\definecolor{second2}{RGB}{175,153,8}%
\definecolor{third2}{RGB}{244,105,102}%
% 蓝色主题备用颜色
\definecolor{structure3}{RGB}{60,113,183}
\definecolor{main3}{RGB}{0,166,82}%
\definecolor{second3}{RGB}{255,134,24}%
\definecolor{third3}{RGB}{0,174,247}%
% 灰色主题备用颜色
\definecolor{structure4}{RGB}{150,150,150}
\definecolor{main4}{RGB}{150,150,150}%
\definecolor{second4}{RGB}{150,150,150}%
\definecolor{third4}{RGB}{150,150,150}%
% 黑色主题备用颜色
\definecolor{structure5}{RGB}{0,0,0}
\definecolor{main5}{RGB}{0,0,0}%
\definecolor{second5}{RGB}{0,0,0}%
\definecolor{third5}{RGB}{0,0,0}%

% 其他通用颜色定义
\definecolor{winered}{rgb}{0.5,0,0}
\definecolor{bule}{RGB}{18,29,57}
\colorlet{coverlinecolor}{second}


% ================= 标题样式选项 =================
\ifdefstring{\ELEGANT@titlestyle}{hang}{\def\style{hang}}{\relax}
\ifdefstring{\ELEGANT@titlestyle}{display}{\def\style{display}}{\relax}

% 参考实现：
% https://tex.stackexchange.com/questions/58506/how-to-make-a-new-command-similar-to-author
% https://pastebin.com/C8W4axzV
% ================= 7. 标题页信息接口 =================
% 定义 email、subtitle、institute、version、logo、cover 等命令，
% 最终由后面的 \maketitle 统一生成封面。
\newcommand\email[1]{\href{mailto:#1}{\nolinkurl{#1}}}

\global\let\@title\@empty
\global\let\@author\@empty
\global\let\@date\@empty
\newcommand{\subtitle}[1]{\gdef\@subtitle{#1}}
\newcommand{\institute}[1]{\gdef\@institute{#1}}
\newcommand{\version}[1]{\gdef\@version{#1}}
\newcommand{\extrainfo}[1]{\gdef\@extrainfo{#1}}

\RequirePackage{mwe}
\newcommand{\logo}[1]{\gdef\@logo{#1}}
\newcommand{\cover}[1]{\gdef\@cover{#1}}

\newcommand{\question}[1]{{\par\citshape #1}\\[0.2ex]}

% ================= 8. 列表、图表标题与显示公式间距 =================
% enumitem 控制 enumerate/itemize；caption 控制图表标题；
% ElegantBook 还在 AtBeginDocument 中设置默认的显示公式与脚注间距。
% 这些均保持模板原样，若个人项目需要调整，应在 preamble.tex 中覆盖。
\RequirePackage{enumerate}
% ================= 列表环境 =================
\RequirePackage[shortlabels,inline]{enumitem}
\setlist{nolistsep}

\RequirePackage[labelfont={bf,color=structurecolor}]{caption}
\captionsetup[table]{skip=3pt}
\captionsetup[figure]{skip=3pt}

\AtBeginDocument{
  \setlength{\abovedisplayskip}{3pt}
  \setlength{\belowdisplayskip}{3pt}
  \RequirePackage[flushmargin,stable]{footmisc}
  \setlength{\footnotesep}{12pt}
}



\RequirePackage{graphicx}

% ================= 【本地修改】水印支持 =================
% 文档类只负责“开/关”水印；具体水印文字、透明度、缩放比例等
% 放在具体项目的 preamble.tex 中通过 \SetWatermarkText 等命令设置。
% ================= 9. 【本地修改】水印机制 =================
% class 只提供 watermark / nowatermark 开关。
% 水印文字、缩放比例、明暗度等个人设置放在 preamble.tex。
\RequirePackage{draftwatermark}
\AtBeginDocument{%
  \ifELEGANT@watermark
    \relax
  \else
    \SetWatermarkText{}%
  \fi
}
%\RequirePackage{enumerate}
\RequirePackage{amsmath,mathrsfs,amsfonts,amssymb}

% ================= 【本地修改】数学学习常用扩展 =================
% ================= 10. 【本地修改】数学学习扩展宏包 =================
% mathtools：amsmath 增强；extarrows：长箭头；arydshln：矩阵虚线；
% tikz-cd：交换图；physics2 + ab：绝对值、范数等数学定界符；
% needspace：避免列表项/标题被分页拆开。
\RequirePackage{mathtools}
\RequirePackage{extarrows}
\RequirePackage{arydshln}
\RequirePackage{tikz-cd}
\RequirePackage{physics2}
\usephysicsmodule{ab}
\RequirePackage{needspace}
% \providecommand\qed{}
% \renewcommand{\qed}{\hfill\ensuremath{\square}}
% ================= 11. 表格与多栏辅助 =================
% booktabs：三线表；multicol：多栏；multirow：跨行单元格；
% fancyvrb / makecell 等用于代码、表格和示例文档排版。
\RequirePackage{booktabs}
\RequirePackage{multicol,multirow}

\RequirePackage{fancyvrb}
\RequirePackage{makecell,lipsum,hologo}
%% ================= 各语言环境中的结构名称 =================


%% 中文模式：目录、图表、定理、证明、附录等名称
% ================= 12. 章节标题、附录与参考文献 =================
% titlesec：chapter/section 等标题格式；appendix：附录；
% biblatex：参考文献，后端由 bibend=biber/bibtex 控制。
\RequirePackage[center]{titlesec}
\RequirePackage[title,titletoc,header]{appendix}


\RequirePackage[
  backend=\ELEGANT@bibend,
  citestyle=\ELEGANT@citestyle,
  bibstyle=\ELEGANT@bibstyle]{biblatex}

% ----- 中文模式：目录、图表、定理等中文名称 -----
\ifdefstring{\ELEGANT@lang}{cn}{
  \renewcommand{\baselinestretch}{1.3}
  \renewcommand{\contentsname}{目录}
  \renewcommand{\figurename}{图}
  \renewcommand{\tablename}{表}
  \renewcommand{\partname}{\color{structurecolor}}
  \renewcommand{\thepart}{第\zhnumber{\arabic{part}}部分}
  \renewcommand{\listfigurename}{插图目录}
  \renewcommand{\listtablename}{表格目录}
  \renewcommand{\bibname}{参考文献}
  \newcommand{\ebibname}{参考文献}
  \renewcommand{\appendixname}{附录}
  \renewcommand{\appendixtocname}{附录}
  \renewcommand{\indexname}{索\hspace{2em}引}
  \newcommand\figref[1]{\textbf{图}~\ref{#1}}
  \newcommand\tabref[1]{\textbf{表}~\ref{#1}}
  \newcommand{\authorname}{\citshape 作者：}
  \newcommand{\institutename}{\citshape 组织：}
  \newcommand{\datename}{\citshape 时间：}
  \newcommand{\versionname}{\citshape 版本：}
  \newcommand{\notename}{笔记}
  \renewcommand*{\proofname}{证明}
  \newcommand{\definitionname}{定义}
  \newcommand{\theoremname}{定理}
  \newcommand{\axiomname}{公理}
  \newcommand{\postulatename}{公设}
  \newcommand{\lemmaname}{引理}
  \newcommand{\propositionname}{命题}
  \newcommand{\corollaryname}{推论}
  \newcommand{\examplename}{例题} %
  \newcommand{\instancename}{示例} %
  \newcommand{\problemname}{问题} % 问题
  \newcommand{\exercisename}{练习} % 练习=习题
  \newcommand{\remarkname}{注}
  \newcommand{\assumptionname}{假设}
  \newcommand{\conclusionname}{结论}
  \newcommand{\solutionname}{解}
  \newcommand{\propertyname}{性质}
  \newcommand{\introductionname}{内容提要}
  \newcommand\bioinfo[2]{\gdef\@bioinfo{{\citshape #1}：#2}}
  \newcommand{\updatename}{更新：}
  \newcommand{\historyname}{版本更新历史}
  \newcommand{\beforechap}{第}
  \newcommand{\afterchap}{章}
}{\relax}


% ----- 英文模式：界面与结构名称 -----
\ifdefstring{\ELEGANT@lang}{en}{
  \setlength\parindent{2em}
  \newcommand\figref[1]{\textbf{Figure}~\ref{#1}}
  \newcommand\tabref[1]{\textbf{Table}~\ref{#1}}
  \renewcommand{\chaptername}{Chapter}
  \renewcommand{\partname}{\color{structurecolor} Part}
  \newcommand{\authorname}{\textbf{Author: }}
  \newcommand{\institutename}{\textbf{Institute: }}
  \newcommand{\datename}{\textbf{Date: }}
  \newcommand{\versionname}{\textbf{Version: }}
  \newcommand{\notename}{Note}
  \newcommand{\proofname}{Proof}
  \newcommand{\problemname}{Problem}
  \newcommand{\definitionname}{Definition}
  \newcommand{\theoremname}{Theorem}
  \newcommand{\axiomname}{Axiom}
  \newcommand{\postulatename}{Postulate}
  \newcommand{\lemmaname}{Lemma}
  \newcommand{\propositionname}{Proposition}
  \newcommand{\corollaryname}{Corollary}
  \newcommand{\examplename}{Example}
  \newcommand{\exercisename}{Exercise}
  \newcommand{\remarkname}{Remark}
  \newcommand{\assumptionname}{Assumption}
  \newcommand{\conclusionname}{Conclusion}
  \newcommand{\solutionname}{Solution}
  \newcommand{\propertyname}{Property}
  \newcommand{\introductionname}{Introduction}
  \renewcommand{\appendixname}{Appendix}
  \newcommand{\ebibname}{Bibliography}
  % \newcommand{\problemsetname}{Exercise}
  \newcommand\bioinfo[2]{\gdef\@bioinfo{\textbf{#1}: #2}}
  \newcommand{\updatename}{Updates:}
  \newcommand{\historyname}{Version History}
}{\relax}

% ----- 意大利语模式 -----
\ifdefstring{\ELEGANT@lang}{it}{
  % \RequirePackage[utf8]{inputenc}
  \RequirePackage[italian]{babel}
  % \RequirePackage[T1]{fontenc}
  \setlength\parindent{2em}
  \newcommand\figref[1]{\textbf{Figura}~\ref{#1}}
  \newcommand\tabref[1]{\textbf{Tabella}~\ref{#1}}
  \renewcommand{\chaptername}{Chapter}
  \newcommand{\authorname}{\textbf{Autore: }}
  \newcommand{\institutename}{\textbf{Istituto: }}
  \newcommand{\datename}{\textbf{Data: }}
  \newcommand{\versionname}{\textbf{Versione: }}
  \newcommand{\notename}{Nota}
  \newcommand{\proofname}{Dimostrazione}
  \newcommand{\problemname}{Problema}
  \newcommand{\definitionname}{Definizione}
  \newcommand{\theoremname}{Teorema}
  \newcommand{\axiomname}{Assioma}
  \newcommand{\postulatename}{Postulato}
  \newcommand{\lemmaname}{Lemma}
  \newcommand{\propositionname}{Proposizione}
  \newcommand{\corollaryname}{Corollario}
  \newcommand{\examplename}{Esempio}
  \newcommand{\exercisename}{Esercizio}
  \newcommand{\remarkname}{Commento}
  \newcommand{\assumptionname}{Assunto}
  \newcommand{\conclusionname}{Conclusione}
  \newcommand{\solutionname}{Soluzione}
  \newcommand{\propertyname}{Proprietà}
  \newcommand{\introductionname}{Introduzione}
  \newcommand\bioinfo[2]{\gdef\@bioinfo{\textbf{#1}: #2}}
  \newcommand{\updatename}{Aggiornamenti:}
  \newcommand{\ebibname}{Bibliografia}
  \newcommand{\historyname}{Cronologia delle versioni}
}{\relax}

% ----- 法语模式 -----
\ifdefstring{\ELEGANT@lang}{fr}{
  \RequirePackage[utf8]{inputenc}
  \RequirePackage[T1]{fontenc}
  \RequirePackage[french]{babel}
  \setlength\parindent{2em}
  \renewcommand\figref[1]{\textbf{Figure}\ref{#1}}
  \renewcommand\tabref[1]{\textbf{Tableau}\ref{#1}}
  \renewcommand{\chaptername}{Chapitre}
  \renewcommand{\authorname}{\textbf{Auteur: }}
  \renewcommand{\institutename}{\textbf{Institut: }}
  \renewcommand{\datename}{\textbf{Date: }}
  \renewcommand{\versionname}{\textbf{Version: }}
  \renewcommand{\notename}{Note}
  \renewcommand{\proofname}{Démonstration}
  \renewcommand{\problemname}{Problème}
  \renewcommand{\definitionname}{Définition}
  \renewcommand{\theoremname}{Théorème}
  \renewcommand{\axiomname}{Axiôme}
  \renewcommand{\postulatename}{Postulat}
  \renewcommand{\lemmaname}{Lemme}
  \renewcommand{\propositionname}{Proposition}
  \renewcommand{\corollaryname}{Corollaire}
  \renewcommand{\examplename}{Exemple}
  \renewcommand{\exercisename}{Exercice}
  \renewcommand{\remarkname}{Remarque}
  \renewcommand{\assumptionname}{Supposition}
  \renewcommand{\conclusionname}{Conclusion}
  \renewcommand{\solutionname}{Solution}
  \renewcommand{\propertyname}{Propriété}
  \renewcommand{\introductionname}{Introduction}
  \renewcommand{\problemsetname}{Exercice}
  \renewcommand{\updatename}{Mises à jour:}
  \newcommand{\ebibname}{Bibliographie}
  \renewcommand{\historyname}{Historique des versions}
}{\relax}


% ----- 荷兰语模式 -----
\ifdefstring{\ELEGANT@lang}{nl}{
  \RequirePackage[dutch]{babel}
  \setlength\parindent{2em}
  \newcommand\figref[1]{\textbf{Figuur}~\ref{#1}}
  \newcommand\tabref[1]{\textbf{Tabel}~\ref{#1}}
  \renewcommand{\chaptername}{Hoofdstuk}
  \newcommand{\authorname}{\textbf{Auteur: }}
  \newcommand{\institutename}{\textbf{Instituut: }}
  \newcommand{\datename}{\textbf{Datum: }}
  \newcommand{\versionname}{\textbf{Versie: }}
  \newcommand{\notename}{Opmerking}
  \newcommand{\proofname}{Bewijs}
  \newcommand{\problemname}{Probleem}
  \newcommand{\definitionname}{Definitie}
  \newcommand{\theoremname}{Stelling}
  \newcommand{\axiomname}{Uitgangspunt}
  \newcommand{\postulatename}{Hypothese}
  \newcommand{\lemmaname}{Lemma}
  \newcommand{\propositionname}{Voorstel}
  \newcommand{\corollaryname}{Gevolgtrekking}
  \newcommand{\examplename}{Voorbeeld}
  \newcommand{\exercisename}{Oefening}
  \newcommand{\remarkname}{Commentaar}
  \newcommand{\assumptionname}{Veronderstelling}
  \newcommand{\conclusionname}{Besluit}
  \newcommand{\solutionname}{Oplossing}
  \newcommand{\propertyname}{Eigenschap}
  \newcommand{\introductionname}{Introductie}
  \newcommand{\problemsetname}{Probleemcomplex}
  \newcommand\bioinfo[2]{\gdef\@bioinfo{\textbf{#1}: #2}}
  \newcommand{\updatename}{Update:}
  \newcommand{\ebibname}{Bibliographie}
  \newcommand{\historyname}{Versie Geschiedenis}
}{\relax}


% ----- 匈牙利语模式 -----
\ifdefstring{\ELEGANT@lang}{hu}{
  \RequirePackage[magyar]{babel}
  \setlength\parindent{2em}
  \newcommand\figref[1]{\ref{#1}~\textbf{Ábra}}
  \newcommand\tabref[1]{\ref{#1}~\textbf{Táblázat}}
  \renewcommand{\chaptername}{Fejezet}
  \newcommand{\authorname}{\textbf{Szerző: }}
  \newcommand{\institutename}{\textbf{Intézmény: }}
  \newcommand{\datename}{\textbf{Dátum: }}
  \newcommand{\versionname}{\textbf{Verziószám: }}
  \newcommand{\notename}{Jegyzet}
  \newcommand{\proofname}{Bizonyítás}
  \newcommand{\problemname}{Probléma}
  \newcommand{\definitionname}{Definíció}
  \newcommand{\theoremname}{Tétel}
  \newcommand{\axiomname}{Axióma}
  \newcommand{\postulatename}{Követelmény}
  \newcommand{\lemmaname}{Lemma}
  \newcommand{\propositionname}{Előzmény}
  \newcommand{\corollaryname}{Következmény}
  \newcommand{\examplename}{Példa}
  \newcommand{\exercisename}{Feladat}
  \newcommand{\remarkname}{Megjegyzés}
  \newcommand{\assumptionname}{Sejtés}
  \newcommand{\conclusionname}{Összefoglalás}
  \newcommand{\solutionname}{Megoldás}
  \newcommand{\propertyname}{Tulajdonság}
  \newcommand{\introductionname}{Bevezetés}
  \newcommand{\problemsetname}{Feladatok}
  \newcommand\bioinfo[2]{\gdef\@bioinfo{\textbf{#1}: #2}}
  \newcommand{\updatename}{Frissítve:}
  \newcommand{\ebibname}{Bibliográfia}
  \newcommand{\historyname}{Korábbi verziók}
}{\relax}

% ----- 德语模式 -----
\ifdefstring{\ELEGANT@lang}{de}{
\setlength\parindent{2em}
\renewcommand{\contentsname}{Inhaltsverzeichnis}
\newcommand\figref[1]{\textbf{Figur}~\ref{#1}}
\newcommand\tabref[1]{\textbf{Tabelle}~\ref{#1}}
\renewcommand{\partname}{\color{structurecolor} Teil}
\renewcommand{\listfigurename}{Abbildungsverzeichnis}
\renewcommand{\listtablename}{Tabellenverzeichnis}
\renewcommand{\bibname}{Bibliographie}
\renewcommand{\appendixname}{Anhang}
\renewcommand{\chaptername}{Kapitel}
\newcommand{\authorname}{\textbf{Autor: }}
\newcommand{\institutename}{\textbf{Institut: }}
\newcommand{\datename}{\textbf{Datum: }}
\newcommand{\versionname}{\textbf{Version: }}
\newcommand{\notename}{Notiz}
\newcommand{\proofname}{Beweis}
\newcommand{\problemname}{Problem}
\newcommand{\definitionname}{Definition}
\newcommand{\theoremname}{Theorem}
\newcommand{\axiomname}{Axiom}
\newcommand{\postulatename}{Postulat}
\newcommand{\lemmaname}{Lemma}
\newcommand{\propositionname}{Behauptung}
\newcommand{\corollaryname}{Corollary}
\newcommand{\examplename}{Beispiel}
\newcommand{\instancename}{Beispiel} %
\newcommand{\exercisename}{\"{U}bung}
\newcommand{\remarkname}{Anmerkung}
\newcommand{\assumptionname}{Annahme}
\newcommand{\conclusionname}{Abschluss}
\newcommand{\solutionname}{Aufl\"{o}sung}
\newcommand{\propertyname}{Property}
\newcommand{\introductionname}{Einleitung}
\newcommand\bioinfo[2]{\gdef\@bioinfo{\textbf{#1}: #2}}
\newcommand{\updatename}{Updates:}
\renewcommand{\ebibname}{Bibliographie}
\newcommand{\historyname}{Version der Geschichte}
}{\relax}

\ifdefstring{\ELEGANT@lang}{es}{
  \RequirePackage[spanish]{babel}
  \setlength\parindent{2em}
  \newcommand\figref[1]{\textbf{Figura}~\ref{#1}}
  \newcommand\tabref[1]{\textbf{Tabla}~\ref{#1}}
  \renewcommand{\spanishchaptername}{Capítulo}
  \newcommand{\authorname}{\textbf{Autor: }}
  \newcommand{\institutename}{\textbf{Instituto: }}
  \newcommand{\datename}{\textbf{Fecha: }}
  \newcommand{\versionname}{\textbf{Versión: }}
  \newcommand{\notename}{Nota}
  \newcommand{\proofname}{Demostración}
  \newcommand{\problemname}{Problema}
  \newcommand{\definitionname}{Definición}
  \newcommand{\theoremname}{Teorema}
  \newcommand{\axiomname}{Axioma}
  \newcommand{\postulatename}{Postulado}
  \newcommand{\lemmaname}{Lema}
  \newcommand{\propositionname}{Proposición}
  \newcommand{\corollaryname}{Corolario}
  \newcommand{\examplename}{Ejemplo}
  \newcommand{\exercisename}{Ejercicio}
  \newcommand{\remarkname}{Comentario}
  \newcommand{\assumptionname}{Asunto}
  \newcommand{\conclusionname}{Conclusión}
  \newcommand{\solutionname}{Solución}
  \newcommand{\propertyname}{Propiedad}
  \newcommand{\introductionname}{Introducción}
  \newcommand{\problemsetname}{Ejercicio}
  \newcommand\bioinfo[2]{\gdef\@bioinfo{\textbf{#1}: #2}}
  \newcommand{\updatename}{Actualización:}
  \newcommand{\ebibname}{Bibliografía}
  \newcommand{\historyname}{Historial de versiones}
}{\relax}

% ----- 蒙古语模式 -----
\ifdefstring{\ELEGANT@lang}{mn}{
  \RequirePackage[mongolian]{babel}
  \RequirePackage[utf8]{inputenc}
  \RequirePackage[T2A]{fontenc}
  \setlength\parindent{2em}
  \newcommand\figref[1]{\textbf{Зураг}\ref{#1}}
  \newcommand\tabref[1]{\textbf{Хүснэгт}\ref{#1}}
  \renewcommand{\chaptername}{Бүлэг}
  \newcommand{\authorname}{\textbf{Зохиогч: }}
  \newcommand{\institutename}{\textbf{Сургууль: }}
  \newcommand{\datename}{\textbf{Огноо: }}
  \newcommand{\versionname}{\textbf{Хувилбар: }}
  \newcommand{\notename}{Тэмдэглэл}
  \newcommand{\proofname}{Баталгаа}
  \newcommand{\problemname}{Бодлого}
  \newcommand{\definitionname}{Тодорхойлолт}
  \newcommand{\theoremname}{Томьёо}
  \newcommand{\axiomname}{Аксиом}
  \newcommand{\postulatename}{Постулат}
  \newcommand{\lemmaname}{Лемма}
  \newcommand{\propositionname}{Таамаглал}
  \newcommand{\corollaryname}{Үр дүн}
  \newcommand{\examplename}{Жишээ}
  \newcommand{\exercisename}{Дасгал}
  \newcommand{\remarkname}{Тайлбар}
  \newcommand{\assumptionname}{Урьдчилсан нөхцөл}
  \newcommand{\conclusionname}{Дүгнэлт}
  \newcommand{\solutionname}{Хариулт}
  \newcommand{\propertyname}{Шинж чанар}
  \newcommand{\introductionname}{Удиртгал}
  \newcommand{\problemsetname}{Бодлогууд}
  \newcommand\bioinfo[2]{\gdef@bioinfo{\textbf{#1}: #2}}
  \newcommand{\updatename}{Шинэчлэлтүүд:}
  \newcommand{\ebibname}{Ном зүй}
  \newcommand{\historyname}{Хувилбарын түүх}
}{\relax}


% ----- 葡萄牙语模式 -----
\ifdefstring{\ELEGANT@lang}{pt}{
  \RequirePackage[portuguese]{babel}
  \setlength\parindent{2em}
  \newcommand\figref[1]{\textbf{Figura}~\ref{#1}}
  \newcommand\tabref[1]{\textbf{Tabela}~\ref{#1}}
  \renewcommand{\chaptername}{Capítulo}
  \newcommand{\authorname}{\textbf{Autor: }}
  \newcommand{\institutename}{\textbf{Instituição: }}
  \newcommand{\datename}{\textbf{Data: }}
  \newcommand{\versionname}{\textbf{Versão: }}
  \newcommand{\notename}{Observação}
  \newcommand{\proofname}{Demostração}
  \newcommand{\problemname}{Problema}
  \newcommand{\definitionname}{Definição}
  \newcommand{\theoremname}{Teorema}
  \newcommand{\axiomname}{Axioma}
  \newcommand{\postulatename}{Postulado}
  \newcommand{\lemmaname}{Lema}
  \newcommand{\propositionname}{Proposição}
  \newcommand{\corollaryname}{Corolário}
  \newcommand{\examplename}{Exemplo}
  \newcommand{\exercisename}{Exercício}
  \newcommand{\remarkname}{Comentário}
  \newcommand{\assumptionname}{Assunto}
  \newcommand{\conclusionname}{Conclusão}
  \newcommand{\solutionname}{Soluções}
  \newcommand{\propertyname}{Propiedade}
  \newcommand{\introductionname}{Introdução}
  \newcommand{\problemsetname}{Questão}
  \newcommand\bioinfo[2]{\gdef\@bioinfo{\textbf{#1}: #2}}
  \newcommand{\updatename}{Atualização: }
  \newcommand{\ebibname}{Bibliografia}
  \newcommand{\historyname}{Histórico de versões}
}{\relax}


% ----- 日语模式 -----
\ifdefstring{\ELEGANT@lang}{jp}{
  \RequirePackage{luatexja}
  \setlength\parindent{2em}
  \renewcommand{\contentsname}{目次}
  \newcommand\figref[1]{\textbf{図}~\ref{#1}}
  \newcommand\tabref[1]{\textbf{表}~\ref{#1}}
  \renewcommand{\partname}{\color{structurecolor}}
  \renewcommand{\listfigurename}{イラストカタログ}
  \renewcommand{\listtablename}{表カタログ}
  \renewcommand{\bibname}{参考文献}
  \renewcommand{\appendixname}{付録}
  \renewcommand{\chaptername}{第 \thechapter 章}
  \newcommand{\authorname}{\textbf{著者：}}
  \newcommand{\institutename}{\textbf{組織団体：}}
  \newcommand{\datename}{\textbf{日付：}}
  \newcommand{\versionname}{\textbf{\citshape バージョン：}}
  \newcommand{\notename}{ノート}
  \newcommand{\proofname}{証明}
  \newcommand{\problemname}{問題}
  \newcommand{\definitionname}{定義}
  \newcommand{\theoremname}{定理}
  \newcommand{\axiomname}{公理}
  \newcommand{\postulatename}{公準}
  \newcommand{\lemmaname}{補題}
  \newcommand{\propositionname}{命題}
  \newcommand{\corollaryname}{系}
  \newcommand{\examplename}{例題}
  \newcommand{\instancename}{例}
  \newcommand{\exercisename}{練習}
  \newcommand{\remarkname}{注}
  \newcommand{\assumptionname}{仮設}
  \newcommand{\conclusionname}{結論}
  \newcommand{\solutionname}{解答}
  \newcommand{\propertyname}{性質}
  \newcommand{\introductionname}{内容概要}
  \newcommand\bioinfo[2]{\gdef\@bioinfo{\textbf{#1}: #2}}
  \newcommand{\updatename}{更新：}
  \newcommand{\ebibname}{参考文献}
  \newcommand{\historyname}{バージョン更新履歴}
}{\relax}


\graphicspath{{./figure/}{./figures/}{./image/}{./images/}{./graphics/}{./graphic/}{./pictures/}{./picture/}}

% ================= 13. TikZ 图形与章节装饰 =================
% TikZ 用于章节标题装饰、圈号、封面元素等。
\RequirePackage{tikz} %% load tikz without tikz
\usetikzlibrary{backgrounds,calc,shadows,positioning,fit}
\newcommand*\circled[1]{\tikz[baseline=(char.base)]{
    \node[shape=circle,draw,inner sep=1pt] (char) {#1};}}
% 圈号数字实现参考：https://tex.stackexchange.com/questions/7032/good-way-to-make-textcircled-numbers

\newcommand*{\eitemi}{\tikz \draw [baseline, ball color=structurecolor,draw=none] circle (2pt);}
\newcommand*{\eitemii}{\tikz \draw [baseline, fill=structurecolor,draw=none,circular drop shadow] circle (2pt);}
\newcommand*{\eitemiii}{\tikz \draw [baseline, fill=structurecolor,draw=none] circle (2pt);}
\setlist[enumerate,1]{label=\color{structurecolor}\arabic*.}
\setlist[enumerate,2]{label=\color{structurecolor}(\alph*).}
\setlist[enumerate,3]{label=\color{structurecolor}\Roman*.}
\setlist[enumerate,4]{label=\color{structurecolor}\Alph*.}
\setlist[itemize,1]{label={\eitemi}}
\setlist[itemize,2]{label={\eitemii}}
\setlist[itemize,3]{label={\eitemiii}}

\RequirePackage{apptools}

% ================= 附录章节标题处理 =================

\ifdefstring{\ELEGANT@lang}{cn}{
  \ifdefstring{\ELEGANT@scheme}{chinese}{
    \newcommand{\xchaptertitle}{第\zhnumber{\arabic{chapter}}章} }{
    \newcommand{\xchaptertitle}{第 \thechapter{} 章}}
  % hyperref 2022/05/16 以后 nameref 的加载顺序发生变化，下面做版本兼容处理。
  % 相关变更：nameref 会更早、且始终被加载（issue #235）。
  \AtBeginDocument{
    \ifELEGANT@Hy@later
      \patchcmd{\NR@chapter}%
      {\numberline{\thechapter}}{\numberline{\xchaptertitle}}%
      {}{\fail}
    \else
      \patchcmd{\Hy@org@chapter}%
      {\numberline{\thechapter}}{\numberline{\xchaptertitle}}%
      {}{\fail}
    \fi
  }
}
{\newcommand{\xchaptertitle}{\chaptername~\thechapter~}}

% 【本地修改】章节编号深度设为 3，即编号到 subsubsection。
\setcounter{secnumdepth}{3}
\titleformat{\chapter}[\style]{\bfseries}{
  \filcenter\LARGE\enspace\bfseries{\color{structurecolor} \IfAppendix{\appendixname\;\thechapter\;}{\xchaptertitle\;}}}{1pt}{
  \LARGE\bfseries\color{structurecolor}\filcenter}[]
\titleformat{\section}[hang]{\bfseries}{
  \Large\bfseries{\color{structurecolor}\thesection}\enspace}{1pt}{%
  \color{structurecolor}\Large\bfseries\filright}
\titleformat{\subsection}[hang]{\bfseries}{
  \large\bfseries\color{structurecolor}\thesubsection\enspace}{1pt}{%
  \color{structurecolor}\large\bfseries\filright}
\titleformat{\subsubsection}[hang]{\bfseries}{
  \large\bfseries\color{structurecolor}\thesubsubsection\enspace}{1pt}{%
  \color{structurecolor}\large\bfseries\filright}

\titlespacing{\chapter}{0pt}{-20pt}{1.3\baselineskip}
%\titlespacing{\subsection}{0pt}{0.5\baselineskip}{-\baselineskip}



%% ================= 定理环境 =================
%% 定义 definition、theorem、postulate、axiom、corollary、lemma、proposition
%% 以及 proof 等数学环境。
\RequirePackage{pifont,manfnt,bbding}
\RequirePackage[many]{tcolorbox}
% \newlength{\normalparindent}
% \setlength{\normalparindent}{\parindent}
\newcommand{\ELEGANT@tmp@tl}{}
\ifELEGANT@usesamecnt
  \newcounter{ELEGANT@samecnt}
  \counterwithin{ELEGANT@samecnt}{\ELEGANT@thmcnt}
  \tcbset{
    new/usesamecnt/.style = {use counter = ELEGANT@samecnt}
  }
\else
  \tcbset{
    new/usesamecnt/.style = {}
  }
\fi

% ================= 14. 定理环境：fancy 模式 =================
% fancy 模式使用 tcolorbox 生成彩色定理框。
% 当前第五章采用 mode=simple，因此不会使用这一分支。
\ifdefstring{\ELEGANT@mode}{fancy}{
  \tcbset{
    common/.style={
        fontupper=\citshape,
        lower separated=false,
        % before upper={\setlength{\parindent}{\normalparindent}},
        coltitle=white,
        colback=gray!5,
        boxrule=0.5pt,
        fonttitle=\bfseries,
        enhanced,
        breakable,
        top=8pt,
        before skip=8pt,
        attach boxed title to top left={
            yshift=-0.11in,
            xshift=0.15in},
        boxed title style={
            boxrule=0pt,
            colframe=white,
            arc=0pt,
            outer arc=0pt},
        separator sign={.},},
    defstyle/.style={
        common,
        colframe=main,
        colback=main!5,
        colbacktitle=main,
        overlay unbroken and last={
            \node[anchor=south east, outer sep=0pt] at (\linewidth-width,0) {
              \textcolor{main}{$\clubsuit$}};}},
    thmstyle/.style={
        common,
        colframe=second,
        colback=second!5,
        colbacktitle=second,
        overlay unbroken and last={
            \node[anchor=south east, outer sep=0pt] at (\linewidth-width,0) {
              \textcolor{second}{$\heartsuit$}};}},
    prostyle/.style={
        common,
        colframe=third,
        colback=third!5,
        colbacktitle=third,
        overlay unbroken and last={
            \node[anchor=south east, outer sep=0pt] at (\linewidth-width,0) {
              \textcolor{third}{$\spadesuit$}};}},
    ELEGANT@title/.code n args={2}
      {

        \tcbset
        {
          title=
            {
              \csname #1name\endcsname~%
              \ifdef{\thetcbcounter}{\thetcbcounter}{}%
              \ifblank{#2}{}{\ (#2)}
            }
        }
      },
    ELEGANT@label/.code n args={2}
      {
        \ifblank{#2}
        {}{\tcbset{label={#1:#2}}}
  }
}

% fancy 模式内部命令：根据环境名、标签前缀和样式创建定理环境。
% #1：环境名；#2：label 前缀；#3：样式名。
% 可用样式：thmstyle、defstyle、prostyle。
% 例如：\ELEGANT@newtheorem{theorem}{thm}{thmstyle}
% 会同时生成有编号 theorem 与无编号 theorem* 两个环境。
% 多语言注意：该命令会自动读取 \theoremname 的定义，
% 因此各语言设置中必须先定义相应的“定理名称”。
\NewDocumentCommand \ELEGANT@newtheorem { m m m O{} }{
  \ifcsundef{#1name}{%
    \ClassError{elegantbook}{%
      \ #1name undefined, \MessageBreak
      Please check in language setting
    }{}
  }{\relax}
  \expandafter\ifblank\expandafter{#4}{
      \tcbset{
        new/usecnt/.style={}
      }
    }{
      \tcbset{
        new/usecnt/.style= {use counter from = #4}
      }
    }
    \DeclareTColorBox[auto counter,number within=\ELEGANT@thmcnt,usesamecnt,usecnt]{#1}{ g o t\label g }{
      common,#3,
      IfValueTF={##1}
        {ELEGANT@title={#1}{##1}}
        {
          IfValueTF={##2}
            {ELEGANT@title={#1}{##2}}
            {ELEGANT@title={#1}{}}
        },
      IfValueT={##4}
        {
          IfBooleanTF={##3}
            {label={##4}}
            {ELEGANT@label={#2}{##4}}
        }
    }
    \DeclareTColorBox{#1*}{ g o }{
      common,#3,
      IfValueTF={##1}
        {ELEGANT@title={#1}{##1}}
        {
          IfValueTF={##2}
            {ELEGANT@title={#1}{##2}}
            {ELEGANT@title={#1}{}}
        },
    }
    }
    % 批量建立 ElegantBook 预设的数学环境。
    % 这些环境显示的名称（如 \definitionname）已在前面的语言区定义。
    \ELEGANT@newtheorem{theorem}{thm}{thmstyle}
    \ELEGANT@newtheorem{definition}{def}{defstyle}
    \ELEGANT@newtheorem{postulate}{pos}{thmstyle}
    \ELEGANT@newtheorem{axiom}{axi}{thmstyle}
    \ELEGANT@newtheorem{corollary}{cor}{thmstyle}
    \ELEGANT@newtheorem{lemma}{lem}{thmstyle}
    \ELEGANT@newtheorem{proposition}{pro}{prostyle}
    }{\relax}


    % ================= 15. 定理环境：simple 模式 =================
    % simple 模式基于 amsthm，更接近传统数学书排版。
    % 当前第五章就是使用这一模式。
    \ifdefstring{\ELEGANT@mode}{simple}{
      \let\openbox\relax
      \RequirePackage{amsthm}
      \let\proof\relax
      % \let\proofname\relax
      \let\endproof\relax

      % simple 模式下声明三类 amsthm 样式：defstyle、thmstyle、prostyle。
      \newtheoremstyle{defstyle}{3pt}{3pt}{\citshape}{-3pt}{
        \bfseries\color{main}}{}{0.5em}{\thmname{#1} \thmnumber{#2} \thmnote{(#3)}}
      \newtheoremstyle{thmstyle}{3pt}{3pt}{\citshape}{-3pt}{
        \bfseries\color{second}}{}{0.5em}{\thmname{#1} \thmnumber{#2} \thmnote{(#3)}}
      \newtheoremstyle{prostyle}{3pt}{3pt}{\citshape}{-3pt}{
        \bfseries\color{third}}{}{0.5em}{\thmname{#1} \thmnumber{#2} \thmnote{(#3)}}

      % simple 模式内部命令：按环境名和样式创建定理环境。
      % #1：环境名；#2：样式名。
      % 可用样式：thmstyle、defstyle、prostyle。
      % 例如：\ELEGANT@newtheorem{theorem}{thmstyle}
      % 会同时生成有编号 theorem 与无编号 theorem* 两个环境。
      % 多语言注意：该命令会自动读取 \theoremname 的定义，
      % 因此各语言设置中必须先定义相应的“定理名称”。
      \NewDocumentCommand \ELEGANT@newtheorem { m m O{} } {
        \theoremstyle{#2}
        \ifELEGANT@usesamecnt
          \newtheorem{#1}[ELEGANT@samecnt]{\csname #1name\endcsname}
        \else
          \ifblank{#3}{
            \newtheorem{#1}{\csname #1name\endcsname}[\ELEGANT@thmcnt]
          }{
            \newtheorem{#1}[#3]{\csname #1name\endcsname}
          }
        \fi
        \newtheorem*{#1*}{\csname #1name\endcsname}
      }

      % 批量建立 ElegantBook 预设的数学环境。
      % 这些环境显示的名称（如 \definitionname）已在前面的语言区定义。
      \ELEGANT@newtheorem{theorem}{thmstyle}
      \ELEGANT@newtheorem{definition}{defstyle}
      \ELEGANT@newtheorem{postulate}{thmstyle}
      \ELEGANT@newtheorem{axiom}{thmstyle}
      \ELEGANT@newtheorem{corollary}{thmstyle}
      \ELEGANT@newtheorem{lemma}{thmstyle}
      \ELEGANT@newtheorem{proposition}{prostyle}
    }{\relax}

    % 用户接口 \elegantnewtheorem：供正文或 preamble 自定义新的定理环境。
    % #1：环境名；#2：显示标题。
    % #3：样式；#4：fancy 模式下的 label 前缀。
    % 可用样式：thmstyle、defstyle、prostyle。
    % simple 模式下如果误传 #4，会报错。
    % fancy 模式下如果未给 #4，会自动令前缀等于 #1，
    % 同时给出警告。
    \newcommand{\ELEGANT@newtheorem@warning}{}
    \newcommand{\ELEGANT@newtheorem@error}{}
    \newcommand{\ELEGANT@usesamecnt@warning}{}
    \NewDocumentCommand\elegantnewtheorem{ m m m g O{} }
    {
      \renewcommand{\ELEGANT@usesamecnt@warning}{
        \ClassWarning{elegantbook}{%
          [#5] won't make sence with option \MessageBreak
          `usesamecnt'.
        }
      }
      \renewcommand{\ELEGANT@newtheorem@warning}{
        \ClassWarning{elegantbook}{%
          Because you didn't provide a prefix, \MessageBreak
          we use #1 as the default prefix. \MessageBreak
          You have to use \MessageBreak
          \ref{#1:label} to refer a \MessageBreak
          \begin{#1}{name}{label} environment. \MessageBreak
        }
      }
      \renewcommand{\ELEGANT@newtheorem@error}{
        \ClassError{elegantbook}{%
          You can't set a prefix in mode ``simple''.\MessageBreak
          Just use \MessageBreak
          \elegantnewtheorem{#1}{#2}{#3}
        }{%
          We are using ``amsthm'' package to \MessageBreak
          generate theorem-like theorem. \MessageBreak
          The 4th parameter ``prefix'' isn't allowed.
        }
      }
      \ifELEGANT@usesamecnt
        \ifblank{#5}{}{\ELEGANT@usesamecnt@warning}
      \else
        \renewcommand{\ELEGANT@tmp@tl}{#5}
      \fi
      \expandafter\def\csname #1name\endcsname{#2}
      \ifdefstring{\ELEGANT@mode}{simple}{%
        \IfValueTF{#4}
        {\ELEGANT@newtheorem@error}
        {\ELEGANT@newtheorem{#1}{#3}[\ELEGANT@tmp@tl]}
      }{\relax}
      \ifdefstring{\ELEGANT@mode}{fancy}{%
        \IfValueTF{#4}
        {
          \def\ELEGANT@temp@prefix{#4}
          \ifdefempty{\ELEGANT@temp@prefix}
          {
            \ELEGANT@newtheorem{#1}{#1}{#3}[\ELEGANT@tmp@tl]
            \ELEGANT@newtheorem@warning
          }
          {\ELEGANT@newtheorem{#1}{#4}{#3}[\ELEGANT@tmp@tl]}
        }
        {
          \ELEGANT@newtheorem{#1}{#1}{#3}[\ELEGANT@tmp@tl]
          \ELEGANT@newtheorem@warning
        }
      }{\relax}
    }

    % main（绿色/定义类）：example、exercise、problem、solution
    % second（橙色/定理类）：proof、note、remark
    % third（蓝色/命题类）：assumption、property、conclusion、custom

    %% 带独立计数器的例题环境
    \newcounter{exam}[chapter]
    \setcounter{exam}{0}
    \renewcommand{\theexam}{\thechapter.\arabic{exam}}
    \newenvironment{example}[1][]{
      \refstepcounter{exam}
      \par\noindent\textbf{\color{main}{\examplename} \theexam #1 }\rmfamily}{
      \par\ignorespacesafterend}

    %% 带独立计数器的练习环境
    \newcounter{exer}[chapter]
    \setcounter{exer}{0}
    \renewcommand{\theexer}{\thechapter.\arabic{exer}}
    \newenvironment{exercise}[1][]{
      \refstepcounter{exer}
      \par\noindent\makebox[-3pt][r]{
        \scriptsize\color{red!90}\HandPencilLeft\quad}
      \textbf{\color{main}{\exercisename} \theexer #1 }\rmfamily}{
      \par\ignorespacesafterend}

    %% 带独立计数器的问题环境
    \newcounter{prob}[chapter]
    \setcounter{prob}{0}
    \renewcommand{\theprob}{\thechapter.\arabic{prob}}
    \newenvironment{problem}[1][]{
      \refstepcounter{prob}
      \par\noindent\textbf{\color{main}{\problemname} \theprob #1 }\rmfamily}{
      \par\ignorespacesafterend}

    \newenvironment{note}{
      \par\noindent\makebox[-3pt][r]{
        \scriptsize\color{red!90}\textdbend\quad}
      \textbf{\color{second}\notename} \citshape}{\par}

    \newenvironment{proof}{
      \par\noindent\textbf{\color{second}\proofname\;}
      \color{black!90}\cfs}{
      % \hfill$\Box$\quad
      \par}

    \newenvironment{solution}{\par\noindent\textbf{\color{main}\solutionname} \citshape}{\par}
    \newenvironment{remark}{\noindent\textbf{\color{second}\remarkname}}{\par}
    \newenvironment{assumption}{\par\noindent\textbf{\color{third}\assumptionname} \citshape}{\par}
    \newenvironment{conclusion}{\par\noindent\textbf{\color{third}\conclusionname} \citshape}{\par}
    \newenvironment{property}{\par\noindent\textbf{\color{third}\propertyname} \citshape}{\par}
    \newenvironment{custom}[1]{\par\noindent\textbf{\color{third} #1} \citshape}{\par}

    % ================= 16. 内容提要与习题集 =================
    % introduction：tcolorbox + 双栏列表；problemset：章末习题集合。
    \RequirePackage{multicol}
    \tcbset{
      introductionsty/.style={
          enhanced,
          breakable,
          colback=structurecolor!10,
          colframe=structurecolor,
          fonttitle=\bfseries,
          colbacktitle=structurecolor,
          fontupper=\citshape,
          attach boxed title to top center={yshift=-3mm,yshifttext=-1mm},
          boxrule=0pt,
          toprule=0.5pt,
          bottomrule=0.5pt,
          top=8pt,
          before skip=8pt,
          sharp corners
        },
    }


    \newenvironment{introduction}[1][\introductionname]{
      \begin{tcolorbox}[introductionsty,title={#1}]
        \begin{multicols}{2}
          \begin{itemize}[label=\textcolor{structurecolor}{\upshape\scriptsize\SquareShadowBottomRight}]}{
          \end{itemize}
        \end{multicols}
      \end{tcolorbox}}

    \RequirePackage{adforn}

    \newenvironment{problemset}[1][\xchaptertitle~\exercisename]{
      \vspace*{10pt}
      \begin{center}
        \phantomsection\addcontentsline{toc}{section}{\texorpdfstring{\xchaptertitle~\exercisename}{\exercisename}}
        % \markboth{#1}{\rightmark}
        \markright{#1}
        \textcolor{structurecolor}{\Large\bfseries\adftripleflourishleft~#1~\adftripleflourishright}
      \end{center}
      \begin{enumerate}}{
      \end{enumerate}}

    \def\relsec{\endgroup start}
    \def\endrelsec{end\begingroup\def \@currenvir {relsec}}

    \ifdefstring{\ELEGANT@result}{noanswer}{
      \AtBeginDocument{
        \excludecomment{solution}
        \excludecomment{proof}
        \excludecomment{inline}
      }
    }{\relax}


    % ================= 页眉页脚默认样式 =================
    % ================= 17. 页眉页脚默认样式 =================
    % 这里是 ElegantBook 原生 fancyhdr 设置。
    % 个人昵称、左右页眉、章首页样式等均在 preamble.tex 中覆盖。
    \RequirePackage{fancyhdr}
    \fancyhf{}


    \fancyfoot[c]{\color{structurecolor}\small\thepage}
    \if@twoside
      \fancyhead[EL]{\color{structurecolor}\cnormal\leftmark}
      \fancyhead[OR]{\color{structurecolor}\cnormal\rightmark}
    \else
      \fancyhead[R]{\color{structurecolor}\cnormal\rightmark}
    \fi

    \renewcommand{\headrule}{\color{structurecolor}\hrule width\textwidth}
    \pagestyle{fancy}
    \renewcommand{\headrulewidth}{1pt}
    % \renewcommand{\headrule}{}
    \fancypagestyle{plain}{\renewcommand{\headrulewidth}{0pt}\fancyhf{}\renewcommand{\headrule}{}}
    \renewcommand{\sectionmark}[1]{\markright{\thesection\, #1}{}}
    \renewcommand{\chaptermark}[1]{\markboth{\IfAppendix{\appendixname\;\thechapter\, #1}{\xchaptertitle\, #1}}{}}


    % ================= 18. 封面 / 标题页 =================
    % 重定义 \maketitle，处理封面图、标题、副标题、作者、机构、版本与 logo。
    \renewcommand*{\maketitle}{%
      \hypersetup{pageanchor=false}
      \pagenumbering{Alph}
      \begin{titlepage}
        \newgeometry{margin = 0in}
        \parindent=0pt
        \ifdefstring{\ELEGANT@device}{normal}{
          \ifcsname @cover\endcsname
            \includegraphics[width=\linewidth]{\@cover}
          \else
            \includegraphics[width=\linewidth]{example-image}
          \fi
        }{\relax}
        \ifdefstring{\ELEGANT@device}{pad}{
          \ifcsname @cover\endcsname
            \includegraphics[trim=0 26bp 0 26bp,clip=true, width=\linewidth]{\@cover}
          \else
            \includegraphics[trim=0 26bp 0 26bp,clip=true, width=\linewidth]{example-image}
          \fi
        }{\relax}
        \vskip-\baselineskip
        {\color{coverlinecolor}\hrule \@width \linewidth \@height .5in}
        \vfill
        \vskip-2ex
        \hspace{2em}
        \parbox{0.8\textwidth}{
          \bfseries\Huge
          \ifcsname @title\endcsname \@title \fi
          \par}
        \vfill
        \vspace{-1.0cm}
        \setstretch{2.5}
        \hspace{2.5em}
        \begin{minipage}[c]{0.67\linewidth}
          {\color{darkgray}\bfseries\Large
            \ifcsname @subtitle\endcsname\@subtitle\\[2ex]\fi}
          \color{gray}\normalsize
          {\renewcommand{\arraystretch}{0.618}
            \begin{tabular}{l}
              % \ifcsname @author\endcsname \authorname \@author\\\fi
              \ifx\@author\empty\else\authorname\cnormal\@author                \\ \fi
              \ifcsname @institute\endcsname \institutename \cnormal\@institute \\ \fi
              % \ifcsname @date\endcsname  \@date\\\fi
              \ifx\@date\empty\else\datename\cnormal\@date                      \\ \fi
              \ifcsname @version\endcsname \cnormal\versionname\@version        \\ \fi
              \ifcsname @bioinfo\endcsname \cnormal\@bioinfo                    \\ \fi
            \end{tabular}}
        \end{minipage}
        \begin{minipage}[c]{0.27\linewidth}
          \begin{tikzpicture}[remember picture,overlay]
            \begin{pgfonlayer}{background}
              \node[opacity=0.8,
                anchor=south east,
                outer sep=0pt,
                inner sep=0pt] at ($(current page.south east) +(-0.8in,1.5in)$) {
                \ifcsname @logo\endcsname\includegraphics[width=4.2cm]{\@logo}\fi};
            \end{pgfonlayer}
          \end{tikzpicture}
        \end{minipage}
        \vfill
        \begin{center}
          \setstretch{1.3}
          \parbox[t]{0.7\textwidth}{\centering \citshape
            \ifcsname @extrainfo\endcsname\@extrainfo\fi}
        \end{center}
        \vfill
      \end{titlepage}
      \restoregeometry
      \thispagestyle{empty}}


    \newcommand{\dollar}{\mbox{\textdollar}}

    % ================= 19. 代码高亮 =================
    % listings 用于模板示例中的 LaTeX/代码排版，不影响普通数学正文。
    \RequirePackage{listings}
    \renewcommand{\ttdefault}{cmtt}
    \lstdefinestyle{mystyle}{
    basicstyle=%
    \ttfamily
    \lst@ifdisplaystyle\small\fi
}

\lstset{basicstyle=\ttfamily,style=mystyle,breaklines=true}

\definecolor{lightgrey}{rgb}{0.9,0.9,0.9}
\definecolor{frenchplum}{RGB}{190,20,83}
\lstset{language=[LaTeX]TeX,
  texcsstyle=*\color{winered},
  numbers=none,
  mathescape=false,
  breaklines=true,
  keywordstyle=\color{winered},
  commentstyle=\color{gray},
  emph={elegantpaper,fontenc,fontspec,xeCJK,citestyle,FiraMono,xunicode,figure,fig,image,img,table,itemize,enumerate,ctex,microtype,description,times,booktabs,tabular,PDFLaTeX,XeLaTeX,type1cm,BibTeX,device,color,mode,lang,amsthm,tcolorbox,titlestyle,cite,ctex,listings,base,math,scheme,toc,esint,chinesefont,amsmath,bibstyle,natbib,pgfornament,addbibresource,printbibliography},
  emphstyle={\color{frenchplum}},
  morekeywords={DeclareSymbolFont,SetSymbolFont,toprule,midrule,bottomrule,institute,version,includegraphics,setmainfont,setsansfont,setmonofont ,setCJKmainfont,setCJKsansfont,setCJKmonofont,RequirePackage,figref,tabref,email,maketitle,keywords,definecolor,extrainfo,logo,cover,subtitle,appendix,chapter,section,hypersetup,mainmatter,frontmatter,tableofcontents,elegantpar,heiti,kaishu,lstset,pagecolor,zhnumber,marginpar,part,equote,marginnote,bioinfo,datechange,listofchange,lvert,lastpage,songti,heiti,fangsong,setCJKfamilyfont,textbf,elegantnewtheorem,thmcnt,colorlet,usesamecnt},
  frame=single,
  tabsize=2,
  rulecolor=\color{structurecolor},
  framerule=0.2pt,
  columns=flexible,
  % backgroundcolor=\color{lightgrey}
}

\newcommand\bmmax{0}
\RequirePackage{bm}

% ================= 版本更新历史 =================
\newcommand{\dateinfoline}[2]{
  \noindent\makebox[0pt][r]{%
    \makebox[-3em][r]{%
      \small
      \textbf{\textcolor{structurecolor}{#1}}}\;\;}%
  {\citshape\updatename\ignorespaces#2}}

\newcommand{\datechange}[2]{%
  \noindent{\makebox[\textwidth][r]{\color{structurecolor}\rule{1.15\textwidth}{.4pt}}}
  \dateinfoline{#1}{#2}
  \addcontentsline{dates}{section}{#1 -- #2}}

\newcommand{\listofchanges}{%
  \begingroup
  \renewcommand{\contentsname}{\historyname}
  \let\old@starttoc\@starttoc
  \def\@starttoc##1{
    \old@starttoc{dates}}
  \tableofcontents%
  \endgroup}
% 版本历史列表实现参考：https://tex.stackexchange.com/questions/472931/list-of-all-dates-in-a-document

\newenvironment{change}{
  \begin{enumerate}[label=\small\protect\circled{\arabic*}]}{
  \end{enumerate}}



% ================= 20. 目录排版 =================
% 重定义 \tableofcontents，支持单栏/双栏目录，并临时调整目录链接颜色。
\renewcommand\tableofcontents{%
  \hypersetup{linktoc=all, linkcolor=black}
  \if@twocolumn
    \@restonecoltrue\onecolumn
  \else
    \@restonecolfalse
  \fi
  \chapter*{\contentsname
    \@mkboth{%
      \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
  \ifdefstring{\ELEGANT@toc}{twocol}{
    \setlength{\columnsep}{2em}
    \begin{multicols}{2}%
      \@starttoc{toc}
    \end{multicols}}{
    \@starttoc{toc}}
  \if@restonecol\twocolumn\fi
  \hypersetup{linkcolor=winered}}


\renewcommand*{\cleardoublepage}{\clearpage\if@twoside \ifodd\c@page\else
      \hbox{}%
      \thispagestyle{empty}%
      \newpage%
      \if@twocolumn\hbox{}\newpage\fi\fi\fi}


% 附录目录修正参考：https://tex.stackexchange.com/questions/56839/chaptername-is-used-even-for-appendix-chapters-in-toc
% ================= 21. 附录目录编号修正 =================
% calc + tocloft 用于修正普通章节与附录在目录中的编号显示。
\RequirePackage{calc}
\RequirePackage[titles]{tocloft}
\ifdefstring{\ELEGANT@lang}{cn}{
  % \renewcommand{\cftchappresnum}{\beforechap\space}
  % \renewcommand{\cftchapaftersnum}{\space\afterchap}
  % \setlength{\cftchapnumwidth}{\widthof{\textbf{附录~十一}}}
  \renewcommand{\numberline}[1]{%
    \@cftbsnum #1\@cftasnum~\@cftasnumb ~%
  }
  \g@addto@macro\appendix{%
    % 在附录中把目录编号恢复为 \thechapter，避免仍使用普通章标题格式。
    \ifELEGANT@Hy@later
      \patchcmd{\NR@chapter}%
      {\numberline{\xchaptertitle}}{\numberline{\thechapter}}%
      {}{\fail}
    \else
      \patchcmd{\Hy@org@chapter}%
      {\numberline{\xchaptertitle}}{\numberline{\thechapter}}%
      {}{\fail}
    \fi
    \addtocontents{toc}{%
      \protect\renewcommand{\protect\cftchappresnum}{\appendixname\space}%
      \protect\renewcommand{\protect\cftchapaftersnum}{}%
    }%
  }
}{
  \renewcommand{\cftchappresnum}{\chaptername\space}
  \renewcommand{\cftchapaftersnum}{\space}
  \setlength{\cftchapnumwidth}{\widthof{\textbf{Appendix~9}}}
  \g@addto@macro\appendix{%
    \addtocontents{toc}{%
      \protect\renewcommand{\protect\cftchappresnum}{\appendixname\space}%
      \protect\renewcommand{\protect\cftchapaftersnum}{}%
      \setlength{\cftchapnumwidth}{\widthof{\textbf{Appendix~999}}}
    }%
  }
}

% 恢复等宽字体为 Latin Modern Typewriter，避免 listings 等环境字体异常。
\renewcommand\ttdefault{lmtt}
