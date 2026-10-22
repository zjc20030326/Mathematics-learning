"""Inspect the compiled PDF and logs; render pages for separate human review."""
from pathlib import Path
import hashlib
import json
import re
import pymupdf
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence"
PREVIEW = EVIDENCE / "preview"
PREVIEW.mkdir(exist_ok=True)
pdf_path = ROOT / "fduthesis.pdf"
doc = pymupdf.open(pdf_path)
log = (EVIDENCE / "build" / "fduthesis.log").read_text(encoding="utf-8", errors="replace")
blg = (EVIDENCE / "build" / "fduthesis.blg").read_text(encoding="utf-8", errors="replace")
texts = [page.get_text() for page in doc]
(EVIDENCE / "thesis-extracted.txt").write_text(
    "\n".join(f"\n--- PDF PAGE {i+1} ---\n{text}" for i, text in enumerate(texts)),
    encoding="utf-8"
)
font_refs = {font[0]: font for page in doc for font in page.get_fonts(full=True)}
fonts = []
for xref, font in font_refs.items():
    font_bytes = len(doc.extract_font(xref)[3])
    glyph_bytes = 0
    if font[2] == "Type3":
        # Matplotlib embeds Type3 glyph programs, not an extractable font file.
        kind, value = doc.xref_get_key(xref, "CharProcs")
        if kind == "xref":
            dictionary = doc.xref_object(int(value.split()[0]))
            for target in re.findall(r"(\d+)\s+0\s+R", dictionary):
                glyph_bytes += len(doc.xref_stream(int(target)) or b"")
    fonts.append({
        "name": font[3], "type": font[2],
        "font_file_bytes": font_bytes, "type3_glyph_program_bytes": glyph_bytes,
        "embedded": font_bytes > 0 or glyph_bytes > 0,
    })
source_preservation = []
for entry in json.loads((EVIDENCE / "local-papers.json").read_text(encoding="utf-8")):
    original = ROOT.parent / "参考论文" / entry["file"]
    exists = original.is_file()
    actual = hashlib.sha256(original.read_bytes()).hexdigest() if exists else None
    source_preservation.append({
        "file": entry["file"], "exists": exists, "sha256": actual,
        "unchanged": actual == entry["sha256"],
    })
warnings = [line.strip() for line in log.splitlines() if re.search(r"(?:Package|LaTeX|Class).*Warning:", line)]
counts = {
    "overfull": log.count("Overfull"),
    "underfull": log.count("Underfull"),
    "missing_characters": log.count("Missing character:"),
    "undefined_references_or_citations": len(re.findall(r"(?:Reference|Citation).*undefined|There were undefined", log)),
    "bibtex_warnings": len(re.findall(r"Warning--", blg)),
}
result = {
    "validation_date": "2026-10-22",
    "pdf": pdf_path.name,
    "pdf_pages": len(doc),
    "pdf_bytes": pdf_path.stat().st_size,
    "pdf_sha256": hashlib.sha256(pdf_path.read_bytes()).hexdigest(),
    "log_counts": counts,
    "warnings": warnings,
    "embedded_fonts": fonts,
    "unresolved_question_mark_pairs_in_pdf": sum(t.count("??") for t in texts),
    "original_reference_pdfs": source_preservation,
    "rendered_review_note": "Contact sheets and selected full pages are for separate visual inspection; their generation does not itself certify layout.",
}
(EVIDENCE / "delivery-audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
for i in [0, 2, 3, 4, 7, 10, 11, 16, 21, 27, 28, 30, 31, 36, 37]:
    if i < len(doc):
        doc[i].get_pixmap(matrix=pymupdf.Matrix(1.35, 1.35)).save(PREVIEW / f"page-{i+1:02d}.png")
for start in range(0, len(doc), 12):
    sheet = Image.new("RGB", (1200, 4*440), "white")
    draw = ImageDraw.Draw(sheet)
    for offset, i in enumerate(range(start, min(start+12, len(doc)))):
        pix = doc[i].get_pixmap(matrix=pymupdf.Matrix(0.48, 0.48))
        page_image = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        xx = (offset % 3)*400 + 45
        yy = (offset // 3)*440 + 25
        sheet.paste(page_image, (xx, yy))
        draw.text((xx, yy-18), f"PDF {i+1}", fill="black")
    sheet.save(PREVIEW / f"contact-{start+1:02d}.png")
print(json.dumps(result, ensure_ascii=False, indent=2))
assert all(item["unchanged"] for item in source_preservation), "An original reference PDF differs."
assert not counts["overfull"] and not counts["missing_characters"]
assert not counts["undefined_references_or_citations"] and not counts["bibtex_warnings"]
assert all(font["embedded"] for font in fonts), "A font is not embedded."
