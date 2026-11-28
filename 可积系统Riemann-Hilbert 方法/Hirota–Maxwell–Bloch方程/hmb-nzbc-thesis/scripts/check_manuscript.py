"""Check the final manuscript sources, build log, PDF outline, bibliography, and links.

Run after build.ps1, passing its output directory with --build-dir.
Visual inspection is recorded separately; this script does not replace it.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import pymupdf

ROOT = Path(__file__).resolve().parents[1]


def first_citation_order(text):
    """Read biblatex/natbib citation keys in source order, ignoring TeX comments."""
    text = re.sub(r'(?<!\\)%[^\n]*', '', text)
    pattern = re.compile(
        r'\\(?:cite(?:p|t|alp|alt|author|yearpar|year|num)?'
        r'|Cite(?:p|t|alp|alt|author)?)(?![A-Za-z])\*?\s*'
        r'(?:\[[^\]]*\]\s*){0,2}\{([^{}]+)\}')
    order = []
    seen = set()
    for match in pattern.finditer(text):
        for key in match[1].split(','):
            key = key.strip()
            assert key and key != '*', ('invalid citation key', key)
            if key not in seen:
                order.append(key)
                seen.add(key)
    return order


def bibliography_checks(first_order, setup_text, bcf, bbl, blg):
    """Match the current biber/biblatex build's GB/T 7714-2025 bibliography to source citations."""
    assert first_order and len(first_order) == len(set(first_order)), \
        'source citation order must contain distinct keys'
    assert re.search(r'bib-backend\s*=\s*biblatex', setup_text), 'expected bib-backend = biblatex in setup.tex'
    assert re.search(r'bib-style\s*=\s*gb7714-2025', setup_text), 'expected bib-style = gb7714-2025 in setup.tex'

    bcf_keys = []
    for key in re.findall(r'<bcf:citekey[^>]*>([^<]+)</bcf:citekey>', bcf):
        key = key.strip()
        if key not in bcf_keys:
            bcf_keys.append(key)
    assert bcf_keys == first_order, ('source/bcf citation order mismatch', first_order, bcf_keys)

    bbl_entries = re.findall(r'\\entry\{([^}]+)\}\{([^}]+)\}', bbl)
    bbl_keys = [k for k, _ in bbl_entries]
    assert bbl_keys == first_order, ('bbl entry order differs from first source citations', first_order, bbl_keys)
    assert 'WARN' not in blg and 'ERROR' not in blg, ('biber warnings/errors present', blg)

    entry_types = {k: t for k, t in bbl_entries}
    return {
        'style': 'gb7714-2025',
        'backend': 'biber',
        'cited_entries': len(first_order),
        'build_bcf_entries': len(bcf_keys),
        'build_bbl_entries': len(bbl_entries),
        'first_citation_order': first_order,
        'build_bbl_order': bbl_keys,
        'entry_types': entry_types,
        'source_bcf_and_bbl_order_match': True,
    }


def source_checks():
    files = [ROOT/'fduthesis.tex', ROOT/'setup.tex', *sorted((ROOT/'section').glob('*.tex'))]
    labels = []
    refs = []
    for path in files:
        text = path.read_text(encoding='utf-8-sig')
        text = re.sub(r'(?<!\\)%[^\n]*', '', text)
        forbidden = [r'\tfrac', r'\begin{smallmatrix}', r'\resizebox',
                     r'\scriptsize', '$$', r'\(', r'\)', r'\phys']
        assert not any(v in text for v in forbidden), (path.name, 'forbidden source form')
        assert len(re.findall(r'(?<!\\)\$', text)) % 2 == 0, path.name
        assert len(re.findall(r'(?<!\\)\\\[', text)) == len(re.findall(r'(?<!\\)\\\]', text)), path.name
        assert len(re.findall(r'\\left\b', text)) == len(re.findall(r'\\right\b', text)), path.name
        stack = []
        for m in re.finditer(r'\\(begin|end)\{([^}]+)\}', text):
            if m[1] == 'begin':
                stack.append(m[2])
            else:
                assert stack and stack.pop() == m[2], (path.name, m[0])
        assert not stack, (path.name, stack)
        braces = 0
        for c in re.sub(r'\\[{}]', '', text):
            braces += (c == '{')-(c == '}')
            assert braces >= 0, path.name
        assert braces == 0, path.name
        labels += re.findall(r'\\label\{([^}]+)\}', text)
        refs += re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', text)
    assert len(labels) == len(set(labels)), 'duplicate labels'
    assert not set(refs)-set(labels), sorted(set(refs)-set(labels))

    main = files[0].read_text(encoding='utf-8-sig')
    setup = files[1].read_text(encoding='utf-8-sig')
    assert len(re.findall(r'\\chapter\{', main)) == 7
    assert r'\listoffigures' in main
    assert r'\printbibliography' in main
    assert r'\PassOptionsToPackage{fleqn}{amsmath}' not in setup
    assert r'\setlength{\mathindent}{0pt}' not in setup

    section_files = files[2:]
    body = '\n'.join(p.read_text(encoding='utf-8-sig') for p in section_files)
    assert not re.search(r'\\(?:textbf|bfseries|paragraph)\b', body)
    assert not re.search(r'\\begin\{(?:theorem|proposition|lemma|definition|remark)\}\[', body)
    assert '故' not in body, 'found forbidden character 故 in thesis text'
    assert len(re.findall(r'\\subsection\{', body)) == 0, 'subsections should be section under chapter hierarchy'
    assert len(re.findall(r'\\begin\{definition\}', body)) == 2
    assert r'S_{\rm sp}' not in body and r'\rho_{\rm sp}' not in body
    assert not re.search(r'(?:[sABCKT]|\\mathcal\s*C)_\{?j,?k\}?|\\sum\\limits_\{k=|j,\s*k\s*=|[1-9]\\le\s*j,\s*k', body), \
        'found index subscript k clashing with spectral parameter k'
    assert r"R_\infty'" not in body and r"-\dfrac{z^2}{q_0^2}\det M(z)" not in body, \
        'found unresolved mathematical typo in manuscript body'

    # Regression checks: no mid-item 其中, no inline pmatrix in body, no non-script \frac{
    env_re = re.compile(r'\\begin\{(theorem|proposition|lemma|definition)\}(.*?)\\end\{\1\}', re.S)
    disp_re = re.compile(r'\\begin\{(equation\*?|align\*?|gather\*?|multline\*?)\}.*?\\end\{\1\}', re.S)
    for path in section_files:
        sec_text = re.sub(r'(?<!\\)%[^\n]*', '', path.read_text(encoding='utf-8-sig'))
        for m in env_re.finditer(sec_text):
            for qm in re.finditer(r'其中', m.group(2)):
                assert not re.search(r'\\item\b', m.group(2)[qm.end():]), (path.name, m.group(1), 'mid-item 其中')
        if path.name != '00-abstract.tex':
            no_disp = disp_re.sub(lambda mm: '\n' * mm.group(0).count('\n'), sec_text)
            assert r'\begin{pmatrix}' not in no_disp, (path.name, 'inline pmatrix outside display math')
            # Ensure \frac{ only appears inside exponents/subscripts (e.g. ^{\frac...} or _{\frac...})
            for fm in re.finditer(r'\\frac\{', sec_text):
                pre = sec_text[max(0, fm.start() - 20):fm.start()]
                assert re.search(r'[\^_]\{[^{}]*$|[\^_]$', pre), (path.name, 'non-script \\frac', sec_text[fm.start()-10:fm.start()+30])
            # Ensure no prose parentheses outside allowed acronyms / enums / subfigure labels
            for lineno, line in enumerate(no_disp.splitlines(), 1):
                s_clean = re.sub(r'\$[^$]*\$', 'MATH', line)
                for pm in re.finditer(r'（[^）]*）|\([^)]*\)', s_clean):
                    tok = pm.group(0)
                    assert re.fullmatch(r'\([0-9a-zA-Z]+\)|（[0-9a-zA-Z\-–—]+）', tok), (path.name, lineno, tok)

    # Supervisor item 3: no prose colons or pre-display colons in section files
    for path in section_files:
        in_math = False
        for lineno, line in enumerate(path.read_text(encoding='utf-8-sig').splitlines(), 1):
            line_nc = re.sub(r'(?<!\\)%[^\n]*', '', line)
            if re.search(r'\\begin\{(equation\*?|align\*?)\}', line_nc):
                in_math = True
            if not in_math:
                stripped = re.sub(r'\\(?:label|ref|eqref|cite|url)\{[^}]*\}', '', line_nc)
                stripped = re.sub(r'\$[^$]*\$', '', stripped)
                assert ':' not in stripped and '：' not in stripped, (path.name, lineno, line)
            if re.search(r'\\end\{(equation\*?|align\*?)\}', line_nc):
                in_math = False

    # Ordered parallel \eqref checks across each section file
    label_pos = {lbl: idx for idx, lbl in enumerate(labels)}
    for path in section_files:
        for lineno, line in enumerate(path.read_text(encoding='utf-8-sig').splitlines(), 1):
            for m in re.finditer(r'\\eqref\{([^}]+)\}(?:[、与及和]|,\s*)\\eqref\{([^}]+)\}', line):
                a, b = m.group(1), m.group(2)
                if a in label_pos and b in label_pos:
                    assert label_pos[a] <= label_pos[b], (path.name, lineno, a, b)

    direct = (ROOT/'section/03-direct-scattering.tex').read_text(encoding='utf-8')
    evolution = (ROOT/'section/04-evolution-resonance.tex').read_text(encoding='utf-8')
    assert r'\label{eq:jost-total-differential}' not in direct
    assert r'\label{eq:theta}' not in direct
    assert not re.search(r'\\(?:theta|phi)\b|V_\\infty|\\Delta V|\\dd t', direct)
    assert r'\label{eq:jost-conjugate-columns}' not in direct
    assert r'\label{eq:jost-sheet-columns}' not in direct
    assert r'\label{eq:jost-total-differential}' in evolution
    assert r'\label{eq:theta}' in evolution
    assert r'\label{eq:volterra}' in direct
    assert r'\label{eq:jost-normalization}' in direct
    assert r'\label{eq:modified-jost}' in direct
    assert r'\chi_+' in direct and r'\chi_-' in direct
    assert r'\label{ass:fields}' in body
    assert r'\label{prop:bloch}' not in body
    assert r'\label{def:general-rhp}' not in body
    assert r'\label{lem:RHP-determinant}' not in body
    assert r'\label{lem:medium-reduction}' not in body
    assert r'\label{eq:norming-time}' in body
    assert r'\label{eq:bloch}' not in body
    assert r'\label{eq:C}' not in body
    refl = (ROOT/'section/06-reflectionless.tex').read_text(encoding='utf-8')
    assert r'\label{eq:detM}' in refl and r'\label{eq:reconstructed-bloch}' not in refl

    citation_order = first_citation_order(body)
    citations = set(citation_order)
    bibliography = (ROOT/'citelist.bib').read_text(encoding='utf-8')
    keys = set(re.findall(r'@\w+\{([^,\s]+)', bibliography))
    assert not citations-keys, sorted(citations-keys)
    return {'tex_files': len(files), 'labels': len(labels), 'references': len(refs),
            'balanced_delimiters_environments_braces': True,
            'seven_chapters': True, 'forbidden_forms_absent': True,
            'zero_gu_characters': True, 'zero_prose_colons': True, 'zero_index_k_conflict': True,
            'sections_under_chapters': len(re.findall(r'\\section\{', body)),
            'spatial_Jost_then_time_evolution_with_single_scattering_matrix': True,
            'cited_bibliography_entries': len(citations),
            'citation_source_files': [p.name for p in section_files],
            'first_citation_order': citation_order}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--visual-review-note', default='Visual review is not inferred by this automated check.')
    args = parser.parse_args()
    sources = source_checks()
    log = (args.build_dir/'fduthesis.log').read_text(encoding='utf-8', errors='replace')
    log_lines = log.splitlines()
    body_line = next((i for i, l in enumerate(log_lines) if '第1章' in l), len(log_lines))
    raw_bad = []
    for pos, line in enumerate(log_lines):
        if not re.search(r'Overfull|Underfull|Warning|undefined|multiply defined|Missing character|^!|Error', line):
            continue
        # Before Chapter 1 only: ignore known fixed template Underfull \hbox from %更新区 Cover I (\makebox[7em][s] at line 6)
        # and fduthesis.cls / %更新区 __fdu_keywords:nNn parbox at end of section/00-abstract.tex (line 17).
        if pos < body_line and 'Underfull \\hbox' in line and (
            'detected at line 6' in line or 'in paragraph at lines 17--' in line
        ):
            continue
        raw_bad.append(line)
    bad = raw_bad
    assert not bad, bad
    pdf = ROOT/'fduthesis.pdf'
    assert pdf.stat().st_mtime >= max(p.stat().st_mtime for p in [ROOT/'fduthesis.tex', ROOT/'setup.tex', *(ROOT/'section').glob('*.tex')])
    document = pymupdf.open(pdf)
    outline = document.get_toc()
    chapters = [row for row in outline if re.match(r'^第 [1-9] 章 ', row[1])]
    if not chapters:
        chapters = [row for row in outline if row[0] == 1 and re.match(r'^(?:[1-9]\s|第\s*[1-9]\s*章)', row[1])]
    assert len(chapters) == 7, outline
    blanks = [
        i+1 for i, page in enumerate(document)
        if len(page.get_text().strip()) < 10
    ]
    assert not blanks, blanks

    eqnum_re = re.compile(r'^\(\d+\.\d+[a-z]?\)$')
    dropped_eqnums = []
    for pno, page in enumerate(document, 1):
        lines_info = []
        for b in page.get_text('dict')['blocks']:
            if b.get('type') != 0:
                continue
            for l in b['lines']:
                txt = ''.join(sp['text'] for sp in l['spans']).strip()
                x0, y0, x1, y1 = l['bbox']
                lines_info.append((y0, y1, txt))
        for i, (y0, y1, txt) in enumerate(lines_info):
            if eqnum_re.match(txt):
                ym = (y0 + y1) / 2
                if not any(j != i and abs(((oy0 + oy1) / 2) - ym) < 12 for j, (oy0, oy1, _) in enumerate(lines_info)):
                    dropped_eqnums.append((pno, txt))
    assert not dropped_eqnums, ('dropped equation numbers found', dropped_eqnums)
    invalid_links = []
    link_count = 0
    for i, page in enumerate(document):
        for link in page.get_links():
            link_count += 1
            if link['kind'] == pymupdf.LINK_GOTO and not 0 <= link['page'] < len(document):
                invalid_links.append((i+1, link))
    assert not invalid_links, invalid_links
    assert link_count > 30, link_count
    bcf = (args.build_dir/'fduthesis.bcf').read_text(encoding='utf-8', errors='replace')
    bbl = (args.build_dir/'fduthesis.bbl').read_text(encoding='utf-8', errors='replace')
    blg = (args.build_dir/'fduthesis.blg').read_text(encoding='utf-8', errors='replace')
    setup_text = (ROOT/'setup.tex').read_text(encoding='utf-8-sig')
    bibliography = bibliography_checks(sources['first_citation_order'], setup_text, bcf, bbl, blg)
    aux = (args.build_dir/'fduthesis.aux').read_text(encoding='utf-8', errors='replace')
    eq224 = re.search(r'\\newlabel\{eq:diagE\}\{\{([^}]+)\}\{([^}]+)\}', aux)
    assert eq224
    report = {
        'sources': sources, 'build': 'XeLaTeX, Biber, XeLaTeX, XeLaTeX',
        'bibliography': bibliography,
        'log_warnings_errors_overfull_underfull': bad, 'pdf_pages': len(document),
        'top_level_chapters': chapters, 'blank_pages': blanks,
        'pdf_links_checked': link_count, 'invalid_internal_links': invalid_links,
        'background_diagonalization_label': 'eq:diagE',
        'current_eq_number_and_printed_page': list(eq224.groups()),
        'pdf_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
        'visual_review': args.visual_review_note,
    }
    (ROOT/'calculations'/'pdf-acceptance.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
