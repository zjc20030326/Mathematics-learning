"""Check the current rewritten source baseline, coverage and equation numbering.

This is a document acceptance tool, not a replacement for mathematical review.
The owning thesis is read only during this check. Reports are written beside
this script. The baseline and numbering are refreshed for the eight-section
rewrite; the companion must retain the current owning text and all its labels.
"""
from pathlib import Path
import collections
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / 'hmb-nzbc-thesis'


def clean(text):
    return re.sub(r'\s+', '', re.sub(r'(?<!\\)%[^\n]*', '', text))


def displays(text):
    pattern = r'\\\[.*?\\\]|\\begin\{(equation\*?|align\*?|gather\*?|multline\*?)\}.*?\\end\{\1\}'
    return [m.group() for m in re.finditer(pattern, text, re.S)]


def numbers(path):
    text = path.read_text(encoding='utf-8-sig')
    return dict(re.findall(r'\\newlabel\{([^{}]+)\}\{\{([^{}]*)\}', text))


def ordered_subsequence(a, b):
    """Allow an explanation to be inserted inside an original prose line."""
    original, expanded = clean(a), clean(b)
    pos = 0
    for ch in original:
        pos = expanded.find(ch, pos)
        if pos < 0:
            return False
        pos += 1
    return True


def main():
    failures = []
    hashes = json.loads((ROOT / 'checks/source-hashes.json').read_text(encoding='utf-8-sig'))
    hashes += json.loads((ROOT / 'checks/source-main-and-pdf-hashes.json').read_text(encoding='utf-8-sig'))
    for entry in hashes:
        path = Path(entry['Path'])
        if hashlib.sha256(path.read_bytes()).hexdigest().upper() != entry['Hash']:
            failures.append(f'Synchronized source baseline changed: {path}')
    sections = []
    all_labels = []
    for original in sorted((SOURCE / 'section').glob('*.tex')):
        target = ROOT / 'chapters' / original.name
        if not target.exists():
            failures.append(f'Missing chapter: {target.name}')
            continue
        a, b = original.read_text(encoding='utf-8'), target.read_text(encoding='utf-8')
        old_labels = re.findall(r'\\label\{([^{}]+)\}', a)
        new_labels = re.findall(r'\\label\{([^{}]+)\}', b)
        all_labels += new_labels
        if old_labels != new_labels:
            failures.append(f'{original.name}: label sequence changed')
        original_displays = [clean(x) for x in displays(a)]
        expanded_displays = collections.Counter(clean(x) for x in displays(b))
        missing = []
        for index, value in enumerate(original_displays, 1):
            if not expanded_displays[value]:
                missing.append(index)
            else:
                expanded_displays[value] -= 1
        if missing:
            failures.append(f'{original.name}: original display blocks missing/changed: {missing}')
        old_env = re.findall(r'\\begin\{(equation|align|gather|multline|theorem|lemma|proposition|remark|definition|corollary)\}', a)
        new_env = re.findall(r'\\begin\{(equation|align|gather|multline|theorem|lemma|proposition|remark|definition|corollary)\}', b)
        if old_env != new_env:
            failures.append(f'{original.name}: numbered environment sequence changed')
        old_headings = re.findall(r'^\\subsection\{[^\n]*', a, re.M)
        new_headings = re.findall(r'^\\subsection\{[^\n]*', b, re.M)
        if [clean(x) for x in old_headings] != [clean(x) for x in new_headings]:
            failures.append(f'{original.name}: subsection sequence differs from current source')
        if original.name.startswith(('07-', '08-')) and re.search(r'\\(?:mathcal|mathscr)\b', b):
            failures.append(f'{original.name}: calligraphic notation remains')
        old_lines = [x.strip() for x in a.splitlines() if x.strip() and not x.lstrip().startswith('%')]
        missing_lines = [x for x in old_lines if x not in b]
        retained = ordered_subsequence(a, b)
        if not retained:
            failures.append(f'{original.name}: original source is not an ordered subsequence')
        # TeX delimiters, environment nesting, comments ignored.
        plain = re.sub(r'(?<!\\)%[^\n]*', '', b)
        if '$$' in plain or '\\(' in plain or '\\)' in plain:
            failures.append(f'{original.name}: wrong math delimiter for TeX')
        stack = []
        for m in re.finditer(r'\\(begin|end)\{([^{}]+)\}', plain):
            if m[1] == 'begin':
                stack.append(m[2])
            elif not stack or stack.pop() != m[2]:
                failures.append(f'{original.name}: environment nesting mismatch near {m[0]}')
        if stack:
            failures.append(f'{original.name}: unclosed environments {stack}')
        if len(re.findall(r'(?<!\\)\\\[', plain)) != len(re.findall(r'(?<!\\)\\\]', plain)):
            failures.append(f'{original.name}: display delimiter mismatch')
        sections.append(dict(file=original.name, original_lines=len(a.splitlines()),
                             expanded_lines=len(b.splitlines()), original_display_blocks=len(original_displays),
                             expanded_display_blocks=len(displays(b)), original_labels=len(old_labels),
                             original_prose_retained_in_order=retained,
                             original_lines_split_by_insertions=len(missing_lines)))
    duplicates = [x for x, n in collections.Counter(all_labels).items() if n > 1]
    if duplicates:
        failures.append(f'Duplicate labels: {duplicates}')
    index = (ROOT / 'chapters/10-formula-index.tex').read_text(encoding='utf-8')
    indexed = re.findall(r'\\eqref\{([^{}]+)\}', index)
    equations = [x for x in all_labels if x.startswith('eq:')]
    if indexed != equations:
        failures.append('Formula index does not match the current source equation sequence')
    # Read every produced TeX file directly, including the guide, index and
    # derivative supplement. Auxiliary template .cls/.def files are unchanged.
    checked_tex = [ROOT / 'hmb-explanation.tex', *sorted((ROOT / 'chapters').glob('*.tex'))]
    for path in checked_tex:
        plain = re.sub(r'(?<!\\)%[^\n]*', '', path.read_text(encoding='utf-8'))
        if '$$' in plain or '\\(' in plain or '\\)' in plain:
            failures.append(f'{path.name}: invalid TeX math delimiter')
        depth = 0
        for token in re.finditer(r'(?<!\\)[{}]', plain):
            depth += 1 if token[0] == '{' else -1
            if depth < 0:
                failures.append(f'{path.name}: closing brace without opening brace')
                break
        if depth:
            failures.append(f'{path.name}: brace depth {depth}')
        if len(re.findall(r'(?<!\\)\$', plain)) % 2:
            failures.append(f'{path.name}: odd inline math delimiter count')
    number_check = None
    if len(sys.argv) > 1:
        old = json.loads((ROOT / 'checks/source-reference.json').read_text(encoding='utf-8'))
        new = numbers(Path(sys.argv[1]))
        mismatch = {k: {'original': v, 'expanded': new.get(k)} for k, v in old.items()
                    if k in all_labels or k.startswith('ch:') if new.get(k) != v}
        if mismatch:
            failures.append(f'Original reference numbers changed: {mismatch}')
        number_check = dict(original_labels_compared=len([k for k in old if k in all_labels or k.startswith('ch:')]), mismatches=mismatch)
    report = dict(source=str(SOURCE), sections=sections, reference_number_check=number_check,
                  indexed_equations=len(indexed),
                  tex_files_checked=[p.name for p in checked_tex],
                  source_matches_sync_baseline=not any('Synchronized source baseline changed' in s for s in failures), failures=failures)
    (ROOT / 'checks/source-coverage.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=True, indent=2))
    return bool(failures)


if __name__ == '__main__':
    raise SystemExit(main())
