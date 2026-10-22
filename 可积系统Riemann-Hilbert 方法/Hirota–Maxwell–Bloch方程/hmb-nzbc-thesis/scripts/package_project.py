"""Create a self-contained source archive, without bundling reference-paper text."""
from pathlib import Path
import hashlib
import json
import zipfile

ROOT = Path(__file__).resolve().parents[1]
archive = ROOT.parent / "hmb-nzbc-thesis-source.zip"
selected = []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file():
        continue
    relative = path.relative_to(ROOT)
    if "__pycache__" in relative.parts or path.suffix == ".pyc":
        continue
    if relative.parts[0] == "evidence":
        if relative.parts[1] == "preview":
            continue
        if relative.parts[1] == "build":
            selected.append(path)
            continue
        if path.suffix == ".json" or path.name.startswith("template-") or path.name == "VISUAL_REVIEW.md":
            selected.append(path)
        continue
    if relative == Path("fduthesis.pdf"):
        continue  # Delivered as a separate directly readable artifact.
    selected.append(path)
manifest = [{
    "path": p.relative_to(ROOT).as_posix(),
    "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
    "bytes": p.stat().st_size,
} for p in selected]
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
    for path in selected:
        bundle.write(path, ROOT.name + "/" + path.relative_to(ROOT).as_posix())
    bundle.writestr(ROOT.name + "/SOURCE_MANIFEST.json", json.dumps(manifest, ensure_ascii=False, indent=2))
with zipfile.ZipFile(archive) as bundle:
    assert bundle.testzip() is None
    for item in manifest:
        content = bundle.read(ROOT.name + "/" + item["path"])
        assert hashlib.sha256(content).hexdigest() == item["sha256"]
report = {
    "archive": str(archive), "file_count": len(selected) + 1,
    "bytes": archive.stat().st_size,
    "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
    "zip_integrity": "passed", "per_file_hashes": "passed",
    "exclusions": ["compiled thesis PDF (delivered separately)", "reference-paper extracted text", "rendered preview images", "Python caches"],
}
(ROOT.parent / "hmb-nzbc-thesis-package.json").write_text(
    json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(report, ensure_ascii=False, indent=2))
