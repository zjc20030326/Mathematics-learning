"""Static checks for the authored thesis, not a substitute for TeX compilation."""
from pathlib import Path
import collections
import json
import re

ROOT = Path(__file__).resolve().parents[1]
files = [ROOT / "fduthesis.tex"] + sorted((ROOT / "section").glob("0[0-7]-*.tex"))
issues = []
labels = []
references = []
citations = []
inline_operator_review = []
for path in files:
    raw = path.read_text(encoding="utf-8")
    text = re.sub(r"(?<!\\)%[^\n]*", "", raw)
    name = path.relative_to(ROOT).as_posix()
    if "$$" in text or "\x60\x60\x60" in text:
        issues.append(f"{name}: Markdown syntax in LaTeX source")
    balance = 0
    for token in re.findall(r"\\[{}]|[{}]", text):
        if token == "{":
            balance += 1
        elif token == "}":
            balance -= 1
        if balance < 0:
            issues.append(f"{name}: premature closing brace")
            break
    if balance:
        issues.append(f"{name}: brace balance {balance}")
    stack = []
    for match in re.finditer(r"\\(begin|end)\{([^}]+)\}", text):
        kind, env = match.groups()
        if kind == "begin":
            stack.append(env)
        elif not stack or stack.pop() != env:
            issues.append(f"{name}: environment nesting mismatch at {env}")
    if stack:
        issues.append(f"{name}: unclosed environments {stack}")
    if len(re.findall(r"(?<!\\)\$", text)) % 2:
        issues.append(f"{name}: odd inline dollar count")
    for left, right in [(r"\[", r"\]"), (r"\left", r"\right")]:
        pattern_left = re.escape(left) + (r"(?![A-Za-z])" if left == r"\left" else "")
        pattern_right = re.escape(right) + (r"(?![A-Za-z])" if right == r"\right" else "")
        if len(re.findall(pattern_left, text)) != len(re.findall(pattern_right, text)):
            issues.append(f"{name}: count mismatch {left}/{right}")
    labels.extend(re.findall(r"\\label\{([^}]+)\}", text))
    references.extend(re.findall(r"\\(?:eqref|ref)\{([^}]+)\}", text))
    for group in re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]+)\}", text):
        citations.extend(x.strip() for x in group.split(","))
    for match in re.finditer(r"(?<!\\)\$([^$]*)\$", text):
        formula = match.group(1)
        if re.search(r"\\(?:sum|lim|int|prod|bigcup|bigcap)\b", formula) and r"\displaystyle" not in formula:
            line = text[:match.start()].count("\n") + 1
            inline_operator_review.append({"file": name, "line": line, "formula": formula})
    if re.search(r"[\x00-\x08\x0b\x0c\x0e-\x1f]", raw):
        issues.append(f"{name}: unexpected control character")

bib = (ROOT / "citelist.bib").read_text(encoding="utf-8")
keys = set(re.findall(r"@\w+\s*\{\s*([^,\s]+)", bib))
for label, count in collections.Counter(labels).items():
    if count > 1:
        issues.append(f"Duplicate label {label}")
for ref in sorted(set(references) - set(labels)):
    issues.append(f"Missing label {ref}")
for cite in sorted(set(citations) - keys):
    issues.append(f"Missing bibliography key {cite}")
for report in ["README.md", "BUILD.md", "RESEARCH_STATUS.md", "REFERENCES.md"]:
    path = ROOT / report
    if not path.exists():
        issues.append(f"Missing report {report}")
        continue
    text = path.read_text(encoding="utf-8")
    # Exclude fenced source examples when checking Markdown math syntax.
    text = re.sub(r"\x60\x60\x60[\s\S]*?\x60\x60\x60", "", text)
    if any(delimiter in text for delimiter in [r"\(", r"\)", r"\[", r"\]"]):
        issues.append(f"{report}: LaTeX-only math delimiter in Markdown prose")
    for target in re.findall(r"\]\(([^)]+)\)", text):
        if "://" not in target and not target.startswith("#"):
            if target == "evidence/source-audit.json":
                continue  # This script writes the linked report below.
            if not (ROOT / target.split("#")[0]).exists():
                issues.append(f"{report}: missing linked file {target}")

result = {
    "authored_tex_files": [p.relative_to(ROOT).as_posix() for p in files],
    "labels": len(labels),
    "reference_commands": len(references),
    "cited_entries": sorted(set(citations)),
    "bibliography_entries": len(keys),
    "issues": issues,
    "inline_large_operators_needing_manual_review": inline_operator_review,
    "scope": "Static checks only; final XeLaTeX/BibTeX logs and rendered PDF must also be inspected.",
}
(ROOT / "evidence" / "source-audit.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(result, ensure_ascii=False, indent=2))
assert not issues, issues
