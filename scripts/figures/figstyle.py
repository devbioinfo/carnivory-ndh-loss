"""Shared figure style (Nature-family journals): 180 mm double column, Liberation Sans (Arial metrics),
7 pt text, thin recessive axes. Strategy colours validated with the dataviz palette checker (3 categorical
slots pass all-pairs CVD separation; autotrophy is the neutral grey baseline)."""
import matplotlib as mpl
import matplotlib.pyplot as plt
MM = 1 / 25.4
W2 = 183 * MM          # double-column width (in)
W1 = 88 * MM
# v3 palette (validated all-pairs: CVD worst ΔE 13.1 deutan; normal-vision ΔE ≥ 20)
C = {'autotrophic': '#1696a8', 'carnivorous': '#d4473a', 'parasitic': '#7a4cb3', 'hemiparasitic': '#7a4cb3',
     'holoparasitic': '#7a4cb3', 'mycoheterotrophic': '#d99d22', 'full_mycoheterotroph': '#d99d22',
     'partial_mycoheterotroph': '#d99d22', 'CAM': '#52514e'}
DARK = {'autotrophic': '#0e6b78', 'carnivorous': '#a8301f', 'parasitic': '#57338a', 'mycoheterotrophic': '#9c6c0c'}
INK, INK2, INK3, GRID = '#0b0b0b', '#3b3a37', '#8a8984', '#e4e3df'
LABEL = {'autotrophic': 'Autotrophic', 'carnivorous': 'Carnivorous', 'hemiparasitic': 'Hemiparasitic',
         'holoparasitic': 'Holoparasitic', 'full_mycoheterotroph': 'Full mycoheterotroph',
         'partial_mycoheterotroph': 'Partial mycoheterotroph'}
def setup():
    mpl.rcParams.update({
        'font.family': 'Liberation Sans', 'font.size': 7, 'axes.titlesize': 7.2, 'axes.labelsize': 7.2,
        'xtick.labelsize': 6.8, 'ytick.labelsize': 6.8, 'legend.fontsize': 6.6, 'axes.linewidth': 0.5,
        'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
        'text.color': INK, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.dpi': 600, 'font.size': 7, 'mathtext.fontset': 'custom',
        'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic', 'mathtext.bf': 'Liberation Sans:bold'})
def panel(ax, letter, x=-0.14, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom', ha='left')
MINFS = 6.0
def legible(fig, minfs=MINFS):
    """Enforce a minimum type size and dark text for readability (Nature: 5–7 pt; we use >= 6 pt)."""
    from matplotlib.text import Text
    for t in fig.findobj(Text):
        if not t.get_text().strip(): continue
        if t.get_fontsize() < minfs: t.set_fontsize(minfs)
        c = t.get_color()
        if isinstance(c, str) and c.lower() in (INK3.lower(), '#8a8984', '#898781'): t.set_color(INK2)
def fit(fig, wmax=183, hmax=247):
    """Shrink the canvas (not the text) until the tight bounding box fits Nature's print area."""
    for _ in range(4):
        fig.canvas.draw(); bb = fig.get_tightbbox(fig.canvas.get_renderer())
        bw, bh = bb.width * 25.4, bb.height * 25.4
        if bw <= wmax + 0.5 and bh <= hmax + 0.5: break
        w, h = fig.get_size_inches()
        fig.set_size_inches(w * min(1, wmax / bw), h * min(1, hmax / bh))
def save(fig, name, dpi=600):
    import os
    legible(fig); fit(fig)
    base = f'/home/claude/work/pan/figures/{name}'
    fig.savefig(base + '.pdf', bbox_inches='tight', pad_inches=0.02)
    fig.savefig(base + '.png', bbox_inches='tight', pad_inches=0.02, dpi=200)
    if os.environ.get('HIRES'):
        fig.savefig(base + '.tif', bbox_inches='tight', pad_inches=0.02, dpi=dpi, pil_kwargs={'compression': 'tiff_lzw'})

# ---- unified visual language (v2) -------------------------------------------------------------
# ndh gene state: one ordinal ramp (lightness carries the order), hue kept clear of strategy colours
GS = {'F': '#2b2d42', 'P': '#a3a7b8', '-': '#f2f1ed', '?': '#ffffff'}
GS_LABEL = {'F': 'Functional', 'P': 'Pseudogene', '-': 'Absent', '?': 'Unresolved'}
LOSS = '#0b0b0b'          # inferred ndh-loss branches: ink, not a hue
STRAT = ['autotrophic', 'carnivorous', 'parasitic', 'mycoheterotrophic']
def light(hexc, f=0.8):
    import matplotlib.colors as mc
    r = np.array(mc.to_rgb(hexc)); return mc.to_hex(r + (1 - r) * f)
# geological time (ICS International Chronostratigraphic Chart v2024/12; Ma)
GEO = [('Cretaceous', 'K', 145.0, 66.0, '#7fc64e'), ('Paleogene', 'Pg', 66.0, 23.03, '#fd9a52'),
       ('Neogene', 'N', 23.03, 2.58, '#ffe619'), ('Quaternary', 'Q', 2.58, 0.0, '#f9f97f')]
EPOCH = [('Paleocene', 66.0, 56.0), ('Eocene', 56.0, 33.9), ('Oligocene', 33.9, 23.03), ('Miocene', 23.03, 5.333),
         ('Pliocene', 5.333, 2.58)]
EVENTS = [(66.0, 'K–Pg'), (33.9, 'EOT')]
def geobands(ax, tmax, axis='x', labels=True, ypos=1.0, alpha_f=0.78, fs=5.4, short_below=12):
    """Pale ICS-coloured period bands on an absolute-age axis (Ma; axis must already be in Ma)."""
    import matplotlib.transforms as mtr
    tr = mtr.blended_transform_factory(ax.transData, ax.transAxes) if axis == 'x' else mtr.blended_transform_factory(ax.transAxes, ax.transData)
    for name, ab, t0, t1, col in GEO:
        a, b = min(t0, tmax), t1
        if a <= b: continue
        (ax.axvspan if axis == 'x' else ax.axhspan)(b, a, color=light(col, alpha_f), lw=0, zorder=0)
        if labels and (a - b) >= 4:
            lab = name if (a - b) >= short_below else ab
            if axis == 'x': ax.text((a + b) / 2, ypos, lab, transform=tr, ha='center', va='bottom', fontsize=fs, color=INK2, clip_on=False)
    for t, lab in EVENTS:
        if t < tmax: (ax.axvline if axis == 'x' else ax.axhline)(t, color=INK3, lw=0.5, ls=(0, (2, 1.5)), zorder=0.5)
import numpy as np
