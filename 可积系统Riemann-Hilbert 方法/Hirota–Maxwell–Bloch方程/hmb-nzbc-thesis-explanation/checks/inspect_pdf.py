"""Render every page into contact sheets and report PDF geometry."""
from pathlib import Path
import json
import math
import sys
try:
    import pymupdf as fitz
except ImportError:
    import fitz
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]


def main():
    pdf = ROOT / 'hmb-explanation.pdf'
    output = Path(sys.argv[1])
    output.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(pdf)
    thumbnails, outside = [], []
    headings = []
    for index, page in enumerate(doc):
        pix = page.get_pixmap(matrix=fitz.Matrix(.45, .45), alpha=False)
        im = Image.frombytes('RGB', (pix.width, pix.height), pix.samples)
        canvas = Image.new('RGB', (im.width+12, im.height+30), '#dddddd')
        canvas.paste(im, (6, 24))
        ImageDraw.Draw(canvas).text((6, 6), f'PDF page {index+1}', fill='black')
        thumbnails.append(canvas)
        for block in page.get_text('dict')['blocks']:
            if block['type'] != 0:
                continue
            for line in block['lines']:
                for span in line['spans']:
                    x0, y0, x1, y1 = span['bbox']
                    if x0 < -1 or y0 < -1 or x1 > page.rect.width+1 or y1 > page.rect.height+1:
                        outside.append({'pdf_page': index+1, 'text': span['text'], 'bbox': span['bbox']})
        text = page.get_text()
        if '\ufffd' in text:
            outside.append({'pdf_page': index+1, 'error': 'replacement character in extracted text'})
    cols, rows = 4, 4
    paths = []
    for start in range(0, len(thumbnails), cols*rows):
        group = thumbnails[start:start+cols*rows]
        width, height = group[0].size
        sheet = Image.new('RGB', (cols*width, math.ceil(len(group)/cols)*height), 'white')
        for n, im in enumerate(group):
            sheet.paste(im, ((n%cols)*width, (n//cols)*height))
        path = output / f'contact-{start//(cols*rows)+1:02d}.png'
        sheet.save(path)
        paths.append(str(path))
    bookmarks = doc.get_toc()
    report = {'pdf': str(pdf), 'pages': len(doc), 'bytes': pdf.stat().st_size,
              'page_geometry': [doc[0].rect.width, doc[0].rect.height],
              'out_of_page_text_spans': outside, 'bookmarks': bookmarks,
              'contact_sheets': paths, 'contact_pages_covered': len(thumbnails)}
    (ROOT / 'checks/pdf-layout.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'pages':len(doc),'out_of_page_text_spans':outside,'contact_sheets':paths}, ensure_ascii=True, indent=2))


if __name__ == '__main__':
    main()
