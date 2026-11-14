"""Check the final manuscript sources, build log, PDF outline and links.

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


def source_checks():
    files = [ROOT/'fduthesis.tex', *sorted((ROOT/'section').glob('*.tex'))]
    labels = []
    refs = []
    for path in files:
        text = path.read_text(encoding='utf-8-sig')
        text = re.sub(r'(?<!\\)%[^\n]*', '', text)
        forbidden = [r'\tfrac', r'\begin{smallmatrix}', r'\resizebox',
                     r'\scriptsize', '$$', r'\(', r'\)', r'\phys']
        assert not any(v in text for v in forbidden), (path.name, 'forbidden source form')
        assert len(re.findall(r'(?<!\\)\$', text)) % 2 == 0, path.name
        # Do not treat the row-spacing command \\[6pt] as a display opener.
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
    assert len(re.findall(r'\\section\{', main)) == 8
    assert r'\chapter{' not in main
    assert r'\begin{abstract' not in main
    body = '\n'.join(p.read_text(encoding='utf-8-sig') for p in files[1:])
    assert r'S_{\rm sp}' not in body and r'\rho_{\rm sp}' not in body
    direct = (ROOT/'section/03-direct-scattering.tex').read_text(encoding='utf-8')
    assert r'\phi' not in direct and r'\theta' not in direct
    assert r'\label{ass:fields}' in body
    assert r'\label{eq:norming-time}' in body
    assert r'\label{lem:RHP-determinant}' in body
    assert r'\label{lem:medium-reduction}' in body
    citations = set()
    for match in re.finditer(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}', body):
        citations.update(match[1].split(','))
    bibliography = (ROOT/'citelist.bib').read_text(encoding='utf-8')
    keys = set(re.findall(r'@\w+\{([^,\s]+)', bibliography))
    assert not citations-keys, sorted(citations-keys)
    return {'tex_files': len(files), 'labels': len(labels), 'references': len(refs),
            'balanced_delimiters_environments_braces': True,
            'eight_continuous_sections': True, 'forbidden_forms_absent': True,
            'spatial_Jost_and_single_scattering_notation': True,
            'cited_bibliography_entries': len(citations)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--visual-review-note', default='Visual review is not inferred by this automated check.')
    args = parser.parse_args()
    sources = source_checks()
    log = (args.build_dir/'fduthesis.log').read_text(encoding='utf-8', errors='replace')
    # Keep matching explicit rather than relying on compiler exit status alone.
    bad = [line for line in log.splitlines() if re.search(
        r'Overfull|Underfull|Warning|undefined|multiply defined|Missing character|^!|Error', line)]
    assert not bad, bad
    pdf = ROOT/'fduthesis.pdf'
    assert pdf.stat().st_mtime >= max(p.stat().st_mtime for p in [ROOT/'fduthesis.tex', *(ROOT/'section').glob('*.tex')])
    document = pymupdf.open(pdf)
    outline = document.get_toc()
    sections = [row for row in outline if re.match(r'^[1-9]\. ', row[1])]
    assert len(sections) == 8 and all(row[0] == 1 for row in sections), sections
    assert not any(re.match(r'^0\.', row[1]) for row in outline)
    blanks = [i+1 for i, page in enumerate(document) if len(page.get_text().strip()) < 20]
    assert not blanks, blanks
    invalid_links = []
    link_count = 0
    for i, page in enumerate(document):
        for link in page.get_links():
            link_count += 1
            if link['kind'] == pymupdf.LINK_GOTO and not 0 <= link['page'] < len(document):
                invalid_links.append((i+1, link))
    assert not invalid_links, invalid_links
    assert link_count > 30, link_count
    aux = (args.build_dir/'fduthesis.aux').read_text(encoding='utf-8', errors='replace')
    eq224 = re.search(r'\\newlabel\{eq:diagE\}\{\{([^}]+)\}\{([^}]+)\}', aux)
    assert eq224
    report = {
        'sources': sources, 'build': 'XeLaTeX, BibTeX, XeLaTeX, XeLaTeX',
        'log_warnings_errors_overfull_underfull': bad, 'pdf_pages': len(document),
        'top_level_sections': sections, 'blank_pages': blanks,
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
