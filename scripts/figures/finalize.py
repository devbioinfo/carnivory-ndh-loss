"""Fit every figure within Nature's print area (max 183 mm wide x 247 mm deep), 600 dpi RGB TIFF (LZW) + vector PDF."""
import glob, os
from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation, PageObject
MM = 72 / 25.4; WMAX, HMAX = 183, 247
files = sorted(glob.glob('Fig*.tif')) + sorted(glob.glob('ed/ED*.tif'), key=lambda x: int(x.split('ED')[1].split('_')[0]))
rows = []
for f in files:
    base = f[:-4]; name = os.path.basename(base)
    im = Image.open(f).convert('RGB'); w, h = im.size; dpi = 600
    wmm, hmm = w / dpi * 25.4, h / dpi * 25.4; s = min(1, WMAX / wmm, HMAX / hmm)
    if s < 1: im = im.resize((round(w * s), round(h * s)), Image.LANCZOS)
    im.save(f'final/{name}.tif', dpi=(dpi, dpi), compression='tiff_lzw')
    r = PdfReader(base + '.pdf'); p = r.pages[0]; pw, ph = float(p.mediabox.width), float(p.mediabox.height)
    sp = min(1, WMAX * MM / pw, HMAX * MM / ph)
    p.add_transformation(Transformation().scale(sp, sp)); p.mediabox.upper_right = (pw * sp, ph * sp); p.cropbox = p.mediabox
    wr = PdfWriter(); wr.add_page(p); wr.write(f'final/{name}.pdf')
    rows.append((name, im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4, s, 6 * s))
for n, a, b, s, fs in rows: print(f'{n:32s} {a:5.0f} x {b:5.0f} mm  scale {s:.2f}  min text {fs:.1f} pt')
