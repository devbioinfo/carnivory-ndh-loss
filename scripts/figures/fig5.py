"""Figure 5 | Mechanistic context: cyclic electron flow routes (established), plastid–nuclear concordance (this study),
evolutionary model (observations solid, hypotheses dashed) and the gene-retention model."""
import csv, numpy as np, collections
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Ellipse, Rectangle, Circle, Polygon
from matplotlib.colors import LinearSegmentedColormap
import figstyle as S, icons as I
S.setup()
P = '/home/claude/work/pan/'
LIN, NDHC, PGRC, HC = '#2b2a27', S.GS['F'], '#2e8b57', '#6d87ad'   # linear flow; NDH route (ndh = indigo); PGR5 route; protons
fig = plt.figure(figsize=(S.W2, 247 * S.MM))
gs = fig.add_gridspec(3, 2, width_ratios=[1.05, 1], height_ratios=[0.92, 1.12, 1.08], hspace=0.3, wspace=0.28, left=0.03, right=0.985, top=0.98, bottom=0.04)

# ---- A: thylakoid membrane ------------------------------------------------------------------------
ax = fig.add_subplot(gs[0, :]); S.panel(ax, 'a', x=-0.005, y=1.0); ax.axis('off'); ax.set_xlim(0, 100); ax.set_ylim(-11.5, 64); ax.set_aspect('auto')
ax.add_patch(Rectangle((0, 37), 100, 23, fc='#f4f7f2', ec='none', zorder=0)); ax.add_patch(Rectangle((0, -4.5), 100, 18.5, fc='#eef1f7', ec='none', zorder=0))
ax.text(1, 57.5, 'STROMA', fontsize=5.6, color=S.INK3, fontweight='bold'); ax.text(1, 6.0, 'LUMEN', fontsize=5.6, color=S.INK3, fontweight='bold')
# lipid bilayer: two leaflets of head groups with tails
for yh, d in [(36.0, -1), (15.0, 1)]:
    for x in np.arange(0.8, 100, 1.6):
        ax.plot([x - 0.25, x - 0.25], [yh + d * 0.9, yh + d * 9.0], color='#d8cfb8', lw=0.35, zorder=1)
        ax.plot([x + 0.25, x + 0.25], [yh + d * 0.9, yh + d * 9.0], color='#d8cfb8', lw=0.35, zorder=1)
        ax.add_patch(Ellipse((x, yh), 1.35, 1.9, fc='#e8c98f', ec='#c9a865', lw=0.25, zorder=1.2))
def grad(patch, c_top, c_bot, z=3, x0=0, x1=100, y0=-4, y1=64):
    ax.add_patch(patch); patch.set_zorder(z)
    g = np.linspace(1, 0, 256)[:, None]
    im = ax.imshow(g, extent=(x0, x1, y0, y1), cmap=LinearSegmentedColormap.from_list('g', [c_bot, c_top]), aspect='auto', zorder=z - 0.01, interpolation='bicubic')
    im.set_clip_path(patch)
def cplx(x, y, w, h, lab, ctop, cbot, ec, fs=6.0, z=3, sub=None, tc=S.INK):
    p = FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0,rounding_size=2.4', fc='none', ec=ec, lw=0.7)
    grad(p, ctop, cbot, z=z); ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0,rounding_size=2.4', fc='none', ec=ec, lw=0.7, zorder=z + 0.1))
    ax.text(x, y + (1.2 if sub else 0), lab, ha='center', va='center', fontsize=fs, fontweight='bold', zorder=z + 0.2, color=tc)
    if sub: ax.text(x, y - 2.6, sub, ha='center', va='center', fontsize=4.5, zorder=z + 0.2, color=S.INK2)
# complexes (established composition; shapes schematic)
cplx(11, 25.5, 11, 25, 'PSII', '#cfe3c4', '#8fb383', '#5d7f52', sub='P680')
cplx(31, 25.5, 10, 23, 'Cyt $b_6f$', '#efe2c4', '#c9ad73', '#8f7642')
cplx(51, 25.5, 11, 25, 'PSI', '#cfe3c4', '#8fb383', '#5d7f52', sub='P700')
cplx(69, 27.5, 15, 21, 'NDH', '#c9ccd8', '#5a5f7a', NDHC, sub='plastid ndhA–K\n+ nuclear subunits', tc='white')
ax.add_patch(FancyBboxPatch((58.8, 17.5), 4.0, 15, boxstyle='round,pad=0,rounding_size=1.5', fc='#a9c79c', ec='#5d7f52', lw=0.5, zorder=2.6))
ax.text(60.8, 12.0, 'Lhca5/6\nlinkers', fontsize=4.0, ha='center', va='top', color=S.INK2, zorder=6)
# ATP synthase: CF0 ring + stalk + CF1 head
cplx(87, 25.5, 9, 17, 'CF$_0$', '#dfe3ea', '#a9b3c4', '#6c7a91', fs=5.2)
ax.add_patch(Rectangle((86.3, 34), 1.4, 6, fc='#a9b3c4', ec='#6c7a91', lw=0.5, zorder=3))
grad(Ellipse((87, 45), 11, 9, fc='none', ec='#6c7a91', lw=0.7), '#e8ebf1', '#a9b3c4', z=3)
ax.add_patch(Ellipse((87, 45), 11, 9, fc='none', ec='#6c7a91', lw=0.7, zorder=3.2)); ax.text(87, 45, 'CF$_1$', fontsize=5.2, ha='center', va='center', zorder=4, fontweight='bold')
ax.text(91.5, 11.5, 'c$_{14}$ ring', fontsize=4.4, ha='left', va='top', color=S.INK2, zorder=6)
# mobile carriers
ax.add_patch(Ellipse((21, 30), 5, 3.4, fc='#f2df9b', ec='#a88f2f', lw=0.5, zorder=5)); ax.text(21, 30, 'PQ', fontsize=4.8, ha='center', va='center', zorder=6)
ax.add_patch(Ellipse((41, 9.5), 5, 3.4, fc='#b8d0ea', ec='#5c7fa8', lw=0.5, zorder=5)); ax.text(41, 9.5, 'PC', fontsize=4.8, ha='center', va='center', zorder=6)
ax.add_patch(Ellipse((55, 46), 5.4, 3.6, fc='#e7c3a8', ec='#a8714a', lw=0.5, zorder=5)); ax.text(55, 46, 'Fd', fontsize=4.8, ha='center', va='center', zorder=6)
ax.add_patch(FancyBboxPatch((62, 51), 8.5, 4.2, boxstyle='round,pad=0,rounding_size=1.2', fc='#e7e5df', ec=S.INK3, lw=0.5, zorder=5)); ax.text(66.25, 53.1, 'FNR', fontsize=4.8, ha='center', va='center', zorder=6)
ax.add_patch(FancyBboxPatch((31, 47.5), 15, 5, boxstyle='round,pad=0,rounding_size=1.4', fc='#d3ecdc', ec=PGRC, lw=0.7, zorder=5)); ax.text(38.5, 50, 'PGR5 / PGRL1', fontsize=5.0, ha='center', va='center', zorder=6, fontweight='bold')
def arr(p, q, col=LIN, lw=0.9, cs='arc3,rad=0', ls='-', z=7, ms=6):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=ms, lw=lw, color=col, connectionstyle=cs, ls=ls, zorder=z, shrinkA=0, shrinkB=0))
# linear electron flow (ink)
ax.text(5, 50.5, r'2H$_2$O $\rightarrow$ O$_2$ + 4H$^+$', fontsize=4.9, zorder=7)
arr((9.5, 49.5), (10, 38.5), lw=0.6, col=S.INK2)
arr((15.5, 31), (18.8, 30.3)); arr((23.5, 30), (26.5, 28.5)); arr((33, 14), (38.6, 10.2)); arr((43.3, 10.4), (47, 13.5))
arr((52.5, 38), (53.8, 44.2)); arr((57.8, 47), (62, 51.5)); arr((70.7, 53.1), (76, 53.1)); ax.text(76.5, 53.1, r'NADPH', fontsize=5.4, va='center', fontweight='bold')
# cyclic routes
arr((52.3, 46.6), (46.2, 49.5), col=PGRC, lw=1.4); arr((34, 47.3), (22.5, 32), col=PGRC, lw=1.4, cs='arc3,rad=0.25')
arr((57.6, 44.8), (65, 38.5), col=NDHC, lw=1.4); arr((61.5, 22), (23.3, 28.2), col=NDHC, lw=1.4, cs='arc3,rad=-0.2')
# protons
def hp(x, y, col=HC):
    ax.add_patch(Circle((x, y), 1.25, fc='white', ec=col, lw=0.5, zorder=7)); ax.text(x, y - 0.05, 'H$^+$', fontsize=3.7, ha='center', va='center', color=col, zorder=8)
for x, lab in [(31, '2 H$^+$/e$^-$ (Q cycle)'), (73, '+2 H$^+$/e$^-$')]:
    arr((x, 41.5), (x, 36.5), col=HC, lw=0.5, ms=4); hp(x, 43.4); arr((x, 14.5), (x, 4.6), col=HC, lw=0.5, ms=4); hp(x, 2.6)
    ax.text(x + 1.8, 2.6, lab, fontsize=4.4, va='center', color=HC, zorder=8)
hp(11, 2.6); arr((11, 14.5), (11, 4.6), col=HC, lw=0.5, ms=4)
arr((87, 4.6), (87, 15), col=HC, lw=0.6, ms=4); hp(87, 2.6)
arr((92.5, 45), (96, 45), lw=0.7); ax.text(96.3, 45, 'ATP', fontsize=5.4, va='center', fontweight='bold'); ax.text(87, 51.2, 'ADP + P$_i$', fontsize=4.4, ha='center', color=S.INK2)
# legend of routes + key numbers
ly = 61.6
for k, (col, lab) in enumerate([(LIN, 'Linear electron flow'), (PGRC, 'PGR5/PGRL1-dependent cyclic flow'), (NDHC, 'NDH-dependent cyclic flow')]):
    xx = [1, 26, 64][k]; ax.add_patch(FancyArrowPatch((xx, ly), (xx + 4.5, ly), arrowstyle='-|>', mutation_scale=5, lw=1.3, color=col)); ax.text(xx + 5.5, ly, lab, fontsize=5.0, va='center')
ax.add_patch(FancyBboxPatch((0.6, -11.2), 70, 6.2, boxstyle='round,pad=0,rounding_size=0.8', fc='white', ec=S.GRID, lw=0.5, zorder=8))
ax.text(1.4, -8.1, r'Linear flow: ATP/NADPH = $\frac{3\ \mathrm{H^+/e^-}\times 2}{14/3\ \mathrm{H^+/ATP}}\approx$ 1.29; Calvin–Benson cycle needs 1.5.' + '\n' + 'Cyclic flow makes up the deficit: PGR5 route 2 H$^+$/e$^-$, NDH route 4 H$^+$/e$^-$.', fontsize=4.6, va='center', zorder=9, linespacing=1.35)
ax.text(100, 57.5, 'Established biochemistry (Kramer & Evans 2011; Shikanai 2016; Strand et al. 2017); shapes schematic', fontsize=4.6, color=S.INK3, ha='right', style='italic')

# ---- B: plastid–nuclear concordance ------------------------------------------------------------
ax = fig.add_subplot(gs[1, 0]); S.panel(ax, 'b', x=0.0, y=1.0)
M = list(csv.DictReader(open(P + 'nuclear/nuclear_matrix.tsv'), delimiter='\t'))
mods = ['NDH_subA', 'NDH_subED', 'NDH_subB', 'NDH_lumen', 'NDH_PSI_link', 'NDH_assembly', 'NDH_expression', 'PGR5_route']
mlab = ['Subcomplex A', 'Subcomplex E/D', 'Subcomplex B', 'Lumen subcomplex', 'PSI linkers (Lhca5/6)', 'Assembly factors', 'Plastid RNA factors', 'PGR5 / PGRL1']
frac = collections.defaultdict(lambda: [0, 0])
for r in M:
    if r['module'] == 'control': continue
    k = (r['species'], r['module']); frac[k][1] += 1; frac[k][0] += (r['call'] == 'present') + 0.5 * (r['call'] == 'partial')
plast = {'Drosera_rotundifolia': '-', 'Drosera_capensis': '-', 'Dionaea_muscipula': '-', 'Cephalotus_follicularis': '-', 'Sarracenia_purpurea': '-', 'Genlisea_aurea': '-',
         'Pinguicula_primuliflora': 'P', 'Drosophyllum_lusitanicum': 'P', 'Striga_asiatica': '-', 'Cuscuta_australis': '-', 'Nepenthes_gracilis': 'F', 'Triphyophyllum_peltatum': 'F',
         'Utricularia_aurea': 'F', 'Utricularia_gibba': 'F', 'Beta_vulgaris': 'F', 'Spinacia_oleracea': 'F', 'Fagopyrum_tataricum': 'F', 'Actinidia_chinensis': 'F', 'Erythranthe_guttata': 'F', 'Ananas_comosus': 'F'}
plab = {'-': 'lost', 'P': 'pseudogenised', 'F': 'intact'}
groups = [('Carnivorous, plastid ndh lost', 'carnivorous', 'sundew', ['Drosera_rotundifolia', 'Drosera_capensis', 'Dionaea_muscipula', 'Cephalotus_follicularis', 'Sarracenia_purpurea', 'Genlisea_aurea', 'Pinguicula_primuliflora']),
          ('Parasitic', 'parasitic', 'broomrape', ['Striga_asiatica', 'Cuscuta_australis']),
          ('Carnivorous, transitional', 'carnivorous', None, ['Drosophyllum_lusitanicum']),
          ('Carnivorous, plastid ndh intact', 'carnivorous', 'nepenthes', ['Nepenthes_gracilis', 'Triphyophyllum_peltatum', 'Utricularia_aurea', 'Utricularia_gibba']),
          ('Autotrophic relatives', 'autotrophic', 'autotroph', ['Beta_vulgaris', 'Spinacia_oleracea', 'Fagopyrum_tataricum', 'Actinidia_chinensis', 'Erythranthe_guttata', 'Ananas_comosus'])]
bins = [(1.0, 1.01, NDHC, 'complete'), (0.5, 1.0, '#6e7290', 'partial ≥50%'), (0.01, 0.5, '#c7c9d4', 'partial <50%'), (-0.01, 0.01, '#f7f6f2', 'none detected')]
def bcol(f):
    for lo, hi, c, _ in bins:
        if lo <= f < hi: return c
    return '#ffffff'
y = 0; ypos = {}
rows_total = sum(len(g[3]) for g in groups) + 1.3 * len(groups)
y = rows_total
for gname, strat, ik, sps in groups:
    y -= 1.3
    ax.text(-7.4, y + 0.15, gname, fontsize=6.5, fontweight='bold', color=S.C[strat] if strat != 'autotrophic' else S.INK2, va='bottom')
    y0 = y
    for sp in sps:
        y -= 1
        for j, m in enumerate(mods):
            a, b = frac[(sp, m)]; f = a / b if b else np.nan
            ax.add_patch(Rectangle((j, y + 0.06), 0.92, 0.88, fc=bcol(f), ec=S.INK3 if f < 0.01 else 'none', lw=0.2))
        st = plast[sp]; ax.add_patch(Rectangle((8.5, y + 0.06), 0.92, 0.88, fc=S.GS[st], ec=S.INK3 if st == '-' else 'none', lw=0.2))
        ax.text(9.6, y + 0.5, plab[st], fontsize=4.6, va='center', color=S.INK2)
        ax.text(-0.25, y + 0.5, sp.replace('_', ' '), ha='right', va='center', fontsize=4.5, style='italic')
    ax.add_patch(Rectangle((-7.45, y + 0.05), 0.15, y0 - y - 0.1, fc=S.C[strat] if strat != 'autotrophic' else S.INK3, ec='none', clip_on=False))
    if ik: I.icon_in(ax, ik, [(12.0 + 7.5) / 21.1, (y + 0.1 - (-1.5)) / (rows_total + 5.2), 0.075, 0.075], S.C[strat])
ax.plot([7.96, 7.96], [0, rows_total], color='white', lw=0)
for j, l in enumerate(mlab): ax.text(j + 0.46, rows_total + 0.3, l, rotation=55, ha='left', va='bottom', fontsize=4.9)
ax.text(8.96, rows_total + 0.3, 'Plastid ndh', rotation=55, ha='left', va='bottom', fontsize=4.9, fontweight='bold')
ax.set_xlim(-7.5, 13.6); ax.set_ylim(-1.5, rows_total + 5.2); ax.axis('off')
for k, (lo, hi, c, l) in enumerate(bins):
    ax.add_patch(Rectangle((-7.3 + k * 5.2, -1.35), 0.6, 0.6, fc=c, ec=S.INK3 if c == '#f7f6f2' else 'none', lw=0.2, clip_on=False)); ax.text(-6.5 + k * 5.2, -1.05, l, fontsize=4.6, va='center')
ax.text(-7.3, -2.0, 'Columns 1–8: nuclear modules, fraction of genes retained (reciprocal\nbest hit; absences checked in the genome). Plastid column: gene-state\ncolours of Figs 1–2. Concordance, 18 species: Fisher P = 2.3 × 10$^{-5}$.', fontsize=4.6, color=S.INK2, va='top', clip_on=False)

# ---- C: model --------------------------------------------------------------------------------------
ax = fig.add_subplot(gs[2, :]); S.panel(ax, 'c', x=-0.005, y=1.0); ax.axis('off'); ax.set_xlim(0, 100); ax.set_ylim(0, 62)
def box(x, y, w, h, t, ec=S.INK2, ls='-', fc='white', fs=5.0, bold=False, tc=S.INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=1.6', fc=fc, ec=ec, lw=0.8, ls=ls))
    ax.text(x + w / 2, y + h / 2, t, ha='center', va='center', fontsize=fs, fontweight='bold' if bold else 'normal', linespacing=1.22, color=tc)
# column headers
for x, t in [(1, 'Strategy'), (26, 'Chloroplast\nconsequence'), (52, 'Selective regime\n(hypothesis)'), (76, 'Genomic outcome\n(observed, this study)')]:
    ax.text(x, 60.5, t, fontsize=5.0, fontweight='bold', color=S.INK2, va='top')
strat = [('Carnivory', 'prey N and P', S.C['carnivorous'], 'trumpet', 41), ('Parasitism', 'host C and N', S.C['parasitic'], 'broomrape', 24.5), ('CAM', 'night CO$_2$ uptake', S.INK, 'succulent', 8)]
for name, sub, col, ik, yy in strat:
    box(1, yy, 21, 12.5, '', ec=col, fc=S.light(col, 0.88) if col != S.INK else '#efeeea')
    ax.text(9.5, yy + 8.2, name, fontsize=5.8, fontweight='bold', color=col if col != S.INK else S.INK, ha='left'); ax.text(9.5, yy + 4.0, sub, fontsize=3.9, color=S.INK2, ha='left')
    I.icon_in(ax, ik, [0.015, (yy + 1.5) / 62, 0.07, 9.5 / 62], col)
box(26, 41, 22, 12.5, 'Prey supplies N and P;\nlower photosynthetic\ncapacity (documented)', fc='white', ec=S.INK3, fs=4.5)
box(26, 24.5, 22, 12.5, 'Host supplies C; reduced\nor no photosynthesis\n(documented)', fc='white', ec=S.INK3, fs=4.5)
box(26, 8, 22, 12.5, 'Altered ATP/NADPH\ndemand by day\n(documented)', fc='white', ec=S.INK3, fs=4.5)
box(52, 18, 20, 30, 'Relaxed purifying\nselection on\nNDH-dependent\ncyclic flow\n\n(ATP balancing and\nphotoprotection in\nfluctuating light\nless critical)', ls=(0, (3, 2)), ec=NDHC, fc='#eeeff3', fs=4.5)
box(76, 43, 23, 13, '', ec=NDHC, fc='white')
ax.text(87.5, 52.6, '1  Plastid ndh decay', fontsize=4.9, ha='center', fontweight='bold', color=NDHC)
for j in range(11): ax.add_patch(Rectangle((78.3 + j * 1.7, 46.2), 1.4, 3.2, fc=S.GS['P' if j % 3 else '-'], ec=S.INK3 if j % 3 == 0 else 'none', lw=0.2))
ax.text(87.5, 44.2, 'rate ×21 (Fig. 3)', fontsize=4.2, ha='center', color=S.INK2)
box(76, 25, 23, 13, '', ec=NDHC, fc='white')
ax.text(87.5, 34.6, '2  Nuclear NDH loss', fontsize=4.9, ha='center', fontweight='bold', color=NDHC)
ax.text(87.5, 29.0, 'follows plastid loss (P = 2.3 × 10$^{-5}$)\nDrosophyllum: plastid decayed,\nnuclear kept', fontsize=4.2, ha='center', va='center', color=S.INK2, linespacing=1.2)
box(76, 6, 23, 14, '', ec=PGRC, fc='#eef7f1')
ax.text(87.5, 16.7, '3  PGR5 route kept', fontsize=4.9, ha='center', fontweight='bold', color=PGRC)
ax.text(87.5, 11.2, 'PGR5: all 22 genomes;\nPGRL1: all carnivores.\nCompensation: hypothesis', fontsize=4.2, ha='center', va='center', color=S.INK2, linespacing=1.2)
def a2(p, q, ls='-', col=S.INK2, cs='arc3,rad=0'):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=6, lw=0.8, color=col, ls=ls, connectionstyle=cs))
for yy in (47.25, 30.75, 14.25): a2((22.3, yy), (25.7, yy))
a2((48.3, 47.25), (51.7, 40), ls=(0, (3, 2))); a2((48.3, 30.75), (51.7, 32), ls=(0, (3, 2))); a2((48.3, 14.25), (51.7, 24), ls=(0, (3, 2)))
a2((72.3, 40), (75.7, 47), ls=(0, (3, 2)), col=NDHC); a2((87.5, 42.7), (87.5, 38.3), col=NDHC); a2((72.3, 26), (75.7, 13), ls=(0, (3, 2)), col=PGRC)
ax.plot([0.5, 6], [2.2, 2.2], color=S.INK2, lw=0.9); ax.text(7, 2.2, 'documented observation', fontsize=4.7, va='center')
ax.plot([30, 35.5], [2.2, 2.2], color=S.INK2, lw=0.9, ls=(0, (3, 2))); ax.text(36.5, 2.2, 'hypothesis or model inference', fontsize=4.7, va='center')
ax.text(99, 2.2, 'CAM: odds ratio ≈ 97, P = 0.019 (Fig. 4c)', fontsize=4.6, va='center', ha='right', color=S.INK2)

# ---- D: retention model ----------------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 1]); S.panel(ax, 'd', x=-0.18, y=1.0)
R = list(csv.DictReader(open(P + 'phaseC/f311_GBM_t3_rates_per_map.csv')))
tt = np.linspace(0, 140, 300)
for reg, lab, col in [('autotrophic', 'Autotrophic', S.C['autotrophic']), ('hemiparasitic', 'Hemiparasitic', S.C['parasitic']), ('carnivorous', 'Carnivorous', S.C['carnivorous'])]:
    lam = np.array([float(r['rate_per_Myr']) for r in R if r['regime'] == reg and r['rate_per_Myr'] not in ('NA', '')])
    curves = np.exp(-np.outer(lam, tt)); lo, md, hi = np.percentile(curves, [2.5, 50, 97.5], axis=0)
    ax.fill_between(tt, lo, hi, color=col, alpha=0.2, lw=0, zorder=1); ax.plot(tt, lo, color=col, lw=0.3, alpha=0.6); ax.plot(tt, hi, color=col, lw=0.3, alpha=0.6)
    ax.plot(tt, md, color='white', lw=2.6, zorder=2); ax.plot(tt, md, color=col, lw=1.5, zorder=3)
    th = np.log(2) / np.median(lam); thlo, thhi = np.log(2) / np.percentile(lam, [97.5, 2.5])
    lp = {'autotrophic': (60, 1.0, 'left', 'bottom'), 'carnivorous': (138, 0.14, 'right', 'bottom'), 'hemiparasitic': (2, 0.2, 'left', 'bottom')}[reg]
    ax.text(lp[0], lp[1], f'{lab}\n$t_{{1/2}}$ = {th:.0f} Myr ({thlo:.0f}–{thhi:.0f})', fontsize=5.4, color=col if reg != 'autotrophic' else S.INK2, ha=lp[2], va=lp[3], fontweight='bold', linespacing=1.15, zorder=6)
T = collections.defaultdict(list)
for r in csv.DictReader(open(P + 'phaseC/f311_GBM_t3_timing_per_map.csv')):
    if r['carnivory_origin_Ma'] not in ('NA', ''): T[r['lineage']].append(float(r['carnivory_origin_Ma']))
retained = {'Lentibulariaceae': 0.5, 'Caryophyllales_carnivores': 0.5, 'Brocchinia': 1, 'Sarraceniaceae': 0, 'Cephalotus': 0, 'Roridula': 0, 'Triantha': 0}
offs = {'Lentibulariaceae': (0, 0.06, 'center'), 'Caryophyllales_carnivores': (0, -0.075, 'center'), 'Brocchinia': (3, -0.06, 'left'), 'Sarraceniaceae': (0, 0.05, 'center'),
        'Cephalotus': (-2.5, 0.05, 'right'), 'Roridula': (2.5, 0.05, 'left'), 'Triantha': (0, -0.07, 'center')}
for L, v in retained.items():
    a = np.median(T[L]); lo_, hi_ = np.percentile(T[L], [2.5, 97.5])
    ax.plot([lo_, hi_], [v, v], color=S.C['carnivorous'], lw=0.6, alpha=0.5, zorder=4)
    ax.plot(a, v, 'D', ms=4.2, mfc='white', mec=S.DARK['carnivorous'], mew=0.9, zorder=5)
    nm = {'Caryophyllales_carnivores': 'Caryophyllales', 'Lentibulariaceae': 'Lentibulariaceae'}.get(L, L); o = offs[L]
    ax.text(a + o[0], v + o[1], nm, fontsize=4.9, ha=o[2], va='center', style='italic' if L in ('Cephalotus', 'Roridula', 'Triantha', 'Brocchinia') else 'normal', zorder=6)
ax.set_xlim(0, 140); ax.set_ylim(-0.12, 1.16); ax.set_xlabel('Time spent in the regime, t (Myr)'); ax.set_ylabel('P(ndh complex intact)')
ax.text(0.98, 0.985, r'$S_r(t)=e^{-\lambda_r t},\ \ t_{1/2}=\ln 2/\lambda_r$', transform=ax.transAxes, ha='right', va='top', fontsize=6.4)
ax.plot([], [], color=S.INK2, lw=1.5, label='model: posterior median, 95% band'); ax.plot([], [], 'D', ms=4.2, mfc='white', mec=S.DARK['carnivorous'], mew=0.9, ls='none', label='observed lineage (1 intact, 0.5 mixed, 0 lost)\nat posterior age of carnivory (95% interval)')
ax.legend(loc='center right', bbox_to_anchor=(1.0, 0.7), fontsize=4.9, handlelength=1.4, borderaxespad=0.2)
ax.text(0.0, -0.16, 't is time spent in a nutritional regime, not geological age, so no geological timescale is shown.', transform=ax.transAxes, fontsize=4.7, color=S.INK3, style='italic', va='top')
ax.yaxis.grid(True, color=S.GRID, lw=0.4); ax.set_axisbelow(True)
S.save(fig, 'Fig5_mechanism')
print('ok')
