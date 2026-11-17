"""Figure 1 | Analytical framework, dated 311-plastome phylogeny, topological uncertainty and clock-model choice (v3)."""
import csv, re, io, json, numpy as np, collections
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from Bio import Phylo
import figstyle as S
S.setup(); P = '/home/claude/work/pan/'
FS = 6.5   # base text size for this figure
fig = plt.figure(figsize=(S.W2, 222 * S.MM))

# ---- a: workflow -----------------------------------------------------------------------------
ax = fig.add_axes([0.02, 0.79, 0.98, 0.2]); ax.axis('off'); ax.set_xlim(0, 100); ax.set_ylim(0, 40)
fig.text(0.0, 0.997, 'a', fontsize=9, fontweight='bold', va='top')
steps = [('Data', '311 plastomes\n50 families\nnew: Triantha,\nByblis, Roridula', r'depth $\geq160\times$'),
 ('Gene scoring', '79 genes\nF / P / − calls\nORF + tblastn\nread-verified', r'lost: $\geq3$ of 11'),
 ('Phylogeny', 'partitioned ML\n5 topologies\nAU, ELW\ngCF, sCF', r'$\ln L=\sum_k \ln P(D_k)$'),
 ('Dating', 'MCMCTree\n124 taxa, 5 fossils\nGaussian branch-\nlength likelihood', r'$\ell\simeq\hat\ell+\mathbf{g}^{\top}\Delta+\frac{1}{2}\Delta^{\top}\mathbf{H}\Delta$'),
 ('Clock choice', 'power posteriors\nGauss–Legendre\nILN vs GBM', r'$\ln Z=\int_0^1 \mathrm{E}_\beta[\ln L]d\beta$'),
 ('Trait models', 'Mk, stochastic\nmaps, Pagel,\nphylo. logistic,\njackknife', r'$\mathbf{P}(t)=e^{\mathbf{Q}t}$'),
 ('Nuclear survey', '22 genomes\n36 genes\nreciprocal best\nhits + tblastn', r'RBH: $a\rightarrow b\rightarrow a$')]
stage = [0, 0, 1, 1, 1, 2, 2]; stc = ['#fbeee9', '#e7f3f5', '#f1ebf8']; ste = [S.C['carnivorous'], S.C['autotrophic'], S.C['parasitic']]
stl = ['Data and annotation', 'Phylogeny and time', 'Comparative inference']
n = len(steps); w = 13.4; gap = (99.6 - n * w) / (n - 1)
for g in range(3):
    ii = [i for i in range(n) if stage[i] == g]; x0 = 0.2 + ii[0] * (w + gap); x1 = 0.2 + ii[-1] * (w + gap) + w
    ax.add_patch(FancyBboxPatch((x0, 35.2), x1 - x0, 2.6, boxstyle='round,pad=0,rounding_size=1', fc=ste[g], ec='none'))
    ax.text((x0 + x1) / 2, 36.5, stl[g], ha='center', va='center', fontsize=FS, color='white', fontweight='bold')
for i, (t, body, eq) in enumerate(steps):
    x = 0.2 + i * (w + gap); g = stage[i]
    ax.add_patch(FancyBboxPatch((x, 1), w, 32.5, boxstyle='round,pad=0,rounding_size=1.2', fc=stc[g], ec=ste[g], lw=0.8))
    ax.text(x + w / 2, 30.6, f'{i + 1}  {t}', ha='center', va='center', fontsize=6.8, fontweight='bold')
    ax.text(x + w / 2, 27.4, body, ha='center', va='top', fontsize=6.6, color=S.INK, linespacing=1.2)
    ax.add_patch(FancyBboxPatch((x + 0.6, 2.2), w - 1.2, 6.2, boxstyle='round,pad=0,rounding_size=0.6', fc='white', ec='none'))
    ax.text(x + w / 2, 5.3, eq, ha='center', va='center', fontsize=6.6 if i != 3 else 6.0)
    if i < n - 1: ax.add_patch(FancyArrowPatch((x + w + 0.1, 17), (x + w + gap - 0.1, 17), arrowstyle='-|>', mutation_scale=7, lw=0.8, color=S.INK2))

# ---- b: circular chronogram ---------------------------------------------------------------------
ax = fig.add_axes([0.0, 0.215, 0.62, 0.62 * 183 / 222], projection='polar')
fig.text(0.005, 0.765, 'b', fontsize=9, fontweight='bold', va='top')
t = Phylo.read(P + 'phaseC/f311_GBM_t3_chronogram.nwk', 'newick')
tv = {r['tip']: r for r in csv.DictReader(open(P + 'phaseC/phaseC_traits_cov311.csv'))}
tips = t.get_terminals(); N = len(tips); dep = t.depths(); H = max(dep[x] for x in tips)
ang = {}; span = 2 * np.pi * 0.955
for i, tp in enumerate(tips): ang[tp] = i / (N - 1) * span
for c in t.get_nonterminals(order='postorder'): ang[c] = (ang[c.clades[0]] + ang[c.clades[-1]]) / 2
ax.set_theta_offset(np.pi / 2); ax.set_theta_direction(-1); ax.axis('off')
for name, abv, t0, t1, col in S.GEO:
    a0, a1 = min(t0, H), t1
    ax.bar(span / 2, a0 - a1, width=span + 0.06, bottom=H - a0, color=S.light(col, 0.82), lw=0, zorder=0)
for c in t.find_clades():
    for ch in c.clades: ax.plot([ang[ch], ang[ch]], [dep[c], dep[ch]], color=S.INK2, lw=0.4, solid_capstyle='butt')
    if c.clades:
        a0, a1 = ang[c.clades[0]], ang[c.clades[-1]]; th = np.linspace(a0, a1, max(3, int(abs(a1 - a0) * 60)))
        ax.plot(th, [dep[c]] * len(th), color=S.INK2, lw=0.4)
for tm in (100, 66, 33.9): ax.plot(np.linspace(0, span, 300), [H - tm] * 300, color=S.INK3, lw=0.5, ls=(0, (2, 1.5)), zorder=0.5)
gx = span + (2 * np.pi - span) / 2
for tm, lab in [(125, '125'), (100, '100'), (66, 'K–Pg'), (33.9, 'EOT'), (0, '0 Ma')]:
    ax.text(gx, H - tm, lab, fontsize=6, color=S.INK, ha='center', va='center', bbox=dict(fc='white', ec='none', pad=0.4, alpha=0.9))
name2tip = {x.name: x for x in tips}
for nut_, k_, pp_, tl_ in json.load(open(P + 'phaseC/loss_clades.json')):
    cl = t.common_ancestor([name2tip[x] for x in tl_]) if len(tl_) > 1 else name2tip[tl_[0]]
    par = t.get_path(cl); p = par[-2] if len(par) > 1 else t.root
    ax.plot([ang[cl]] * 2, [dep[p], dep[cl]], color='white', lw=2.8, solid_capstyle='butt', zorder=3)
    ax.plot([ang[cl]] * 2, [dep[p], dep[cl]], color=S.LOSS, lw=1.6, solid_capstyle='butt', zorder=3.1)
R0 = H + 4; dth = span / (N - 1)
grp = lambda n: 'carnivorous' if n == 'carnivorous' else 'parasitic' if 'parasit' in n else 'mycoheterotrophic' if 'myco' in n else 'autotrophic'
AUT = S.light(S.C['autotrophic'], 0.55)
for tp in tips:
    r = tv[tp.name]; a = ang[tp]; g = grp(r['nutrition'])
    ax.bar(a, 5.5, width=dth, bottom=R0, color=S.C[g] if g != 'autotrophic' else AUT, lw=0)
    ax.bar(a, 1.7, width=dth, bottom=R0 + 6.1, color=S.INK if r['CAM'] == '1' else '#f2f1ed', lw=0)
    ax.bar(a, 1.7, width=dth, bottom=R0 + 8.1, color='#5b8db8' if r['aquatic'] == '1' else '#f2f1ed', lw=0)
    calls = r['ndh_calls'] if tp.name != 'Byblis_gigantea' else ''.join('F' if g_ in 'CGIJK' else '?' for g_ in 'ABCDEFGHIJK')
    for j, ch in enumerate(calls): ax.bar(a, 1.0, width=dth, bottom=R0 + 10.4 + j * 1.1, color=S.GS.get(ch, '#ffffff'), lw=0)
def arc(names, label, rr, it=False, col=S.INK, fs=6.6):
    tt = [name2tip[x] for x in names if x in name2tip]
    a = sorted(ang[x] for x in tt); a0, a1 = a[0] - dth / 2, a[-1] + dth / 2
    th = np.linspace(a0, a1, 30); ax.plot(th, [rr] * 30, color=col, lw=1.0)
    am = (a0 + a1) / 2; deg = np.degrees(am); rot = 90 - deg if deg <= 180 else 270 - deg
    ax.text(am, rr + 3, label, rotation=rot, ha='left' if deg <= 180 else 'right', va='center', fontsize=fs, style='italic' if it else 'normal',
            rotation_mode='anchor', color=col if col != S.INK2 else S.INK, fontweight='bold' if it else 'normal')
fam = collections.defaultdict(list)
for tp in tips: fam[tv[tp.name]['family']].append(tp.name)
RL = R0 + 24.5
for f, lab in [('Lentibulariaceae', 'Lentibulariaceae'), ('Droseraceae', 'Droseraceae'), ('Nepenthaceae', 'Nepenthaceae'), ('Sarraceniaceae', 'Sarraceniaceae'),
               ('Orobanchaceae', 'Orobanchaceae'), ('Ericaceae', 'Ericaceae'), ('Bromeliaceae', 'Bromeliaceae'), ('Cactaceae', 'Cactaceae'), ('Tofieldiaceae', 'Tofieldiaceae (incl. Triantha)')]:
    arc(fam[f], lab, RL, col=S.INK2)
for nm, lab, dr in [('Cephalotus_follicularis', 'Cephalotus', 0), ('Roridula_gorgonias', 'Roridula', 0), ('Byblis_gigantea', 'Byblis', 50),
                    ('Drosophyllum_lusitanicum', 'Drosophyllum', 0), ('Triphyophyllum_peltatum', 'Triphyophyllum', 0)]:
    arc([nm], lab, RL + dr, it=True, col=S.DARK['carnivorous'])
ax.set_ylim(0, RL + 50)
# legend for b
lx, ly = 0.015, 0.175
blocks = [('Strategy', [(S.C['carnivorous'], 'Carnivorous'), (S.C['parasitic'], 'Parasitic'), (S.C['mycoheterotrophic'], 'Mycoheterotrophic'), (AUT, 'Autotrophic')]),
          ('ndh gene state', [(S.GS['F'], 'Functional'), (S.GS['P'], 'Pseudogene'), (S.GS['-'], 'Absent')]),
          ('Period', [(S.light(S.GEO[0][4], 0.82), 'Cretaceous'), (S.light(S.GEO[1][4], 0.82), 'Paleogene'), (S.light(S.GEO[2][4], 0.82), 'Neogene')]),
          ('Annotation', [(S.INK, 'CAM'), ('#5b8db8', 'Submerged aquatic'), (None, 'ndh loss (≥50% maps)')])]
xs_ = [0, 0.165, 0.305, 0.42]
for bi, (tt, items) in enumerate(blocks):
    cx0 = lx + xs_[bi]; fig.text(cx0, ly + 0.012, tt, fontsize=6.5, fontweight='bold', va='bottom')
    for k, (col, lab) in enumerate(items):
        cy = ly - k * 0.017
        if col is None: fig.add_artist(plt.Line2D([cx0, cx0 + 0.014], [cy + 0.005] * 2, color=S.LOSS, lw=1.8, transform=fig.transFigure))
        else: fig.patches.append(FancyBboxPatch((cx0, cy), 0.014, 0.01, boxstyle='square,pad=0', fc=col, ec=S.INK3 if col == S.GS['-'] else 'none', lw=0.3, transform=fig.transFigure))
        fig.text(cx0 + 0.019, cy + 0.005, lab, fontsize=6.5, va='center')
fig.text(lx, ly - 0.075, 'Rings, inner to outer: strategy; CAM; aquatic; ndhA (inner) to ndhK (outer).\nDashed rings: 100 Ma, K–Pg (66.0 Ma), EOT (33.9 Ma). Periods: ICS 2024.\nNo Jurassic band: root (137.9 Ma) is younger than the Jurassic–Cretaceous boundary (145.0 Ma).',
         fontsize=6, color=S.INK2, va='top', linespacing=1.3)

# ---- c: topology scenarios ------------------------------------------------------------------------
gs = fig.add_gridspec(2, 2, left=0.665, right=0.985, top=0.735, bottom=0.045, hspace=0.3, height_ratios=[1, 1.1])
gc = gs[0, :].subgridspec(1, 2, width_ratios=[1.55, 1], wspace=0.08)
ax = fig.add_subplot(gc[0]); fig.text(0.6, 0.765, 'c', fontsize=9, fontweight='bold', va='top')
sc = [('ML', 'Ericaceae + (Roridula,\nActinidiaceae)', 0, 0.24, 1.0), ('SR', 'Sarraceniaceae + (Roridula,\nActinidiaceae)', 8.69, 0.375, 0.437),
      ('H4', 'Nepenthes + (Drosophyllum,\nDioncophyllaceae)', 21.02, 0.0049, 0.0429), ('H5', 'Drosera monophyletic', 8.36, 0.154, 0.342),
      ('H6', 'Cephalotus outside\nBrunelliaceae + Elaeocarpaceae', 2.26, 0.225, 0.495)]
TOP = '#4f6d8a'
for i, (k, d, dl, elw, au) in enumerate(sc):
    y = len(sc) - 1 - i; rej = au < 0.05
    if rej: ax.barh(y + 0.15, dl, height=0.38, color='white', ec=S.C['carnivorous'], hatch='////', lw=0.6)
    else: ax.barh(y + 0.15, max(dl, 0.2), height=0.38, color=TOP, lw=0)
    ax.text(-0.8, y + 0.15, k, ha='right', va='center', fontsize=7, fontweight='bold', color=S.DARK['carnivorous'] if rej else S.INK)
    ax.text(0.6 + dl, y + 0.15, f'{dl:.1f}' if dl else 'ref.', va='center', fontsize=6.2)
    ax.text(0, y - 0.3, d, fontsize=5.9, color=S.INK2, va='center', linespacing=1.05)
ax.set_yticks([]); ax.spines['left'].set_visible(False); ax.set_xlim(0, 26); ax.set_ylim(-0.8, len(sc) - 0.4)
ax.set_xlabel(r'$\Delta\ln L$ vs ML tree', fontsize=6.8); ax.tick_params(labelsize=6.2)
ax2 = fig.add_subplot(gc[1], sharey=ax); ax2.patch.set_alpha(0)
for i, (k, d, dl, elw, au) in enumerate(sc):
    y = len(sc) - 1 - i + 0.15; rej = au < 0.05
    ax2.plot([0.01, au], [y, y], color=S.GRID, lw=0.9); ax2.plot(au, y, 'o', ms=5, mfc='white' if rej else TOP, mec=S.DARK['carnivorous'] if rej else TOP, mew=0.9)
    ax2.text(1.25, y + 0.36, f'AU {au:.3f} · ELW {elw:.2f}', fontsize=6, color=S.INK, va='center', ha='right', bbox=dict(fc='white', ec='none', pad=0.2))
ax2.axvline(0.05, color=S.INK2, lw=0.7, ls=(0, (2, 1.5))); ax2.text(0.055, -0.7, r'$\alpha$ = 0.05', fontsize=6, ha='left')
ax2.set_xscale('log'); ax2.set_xlim(0.01, 1.3); ax2.set_xlabel('AU test P', fontsize=6.8); ax2.tick_params(axis='y', left=False, labelleft=False); ax2.tick_params(labelsize=6.2)
ax2.spines['left'].set_visible(False)
ax.text(0, len(sc) - 0.15, 'Carnivore placements (295 taxa); hatched = rejected', fontsize=6.3, va='bottom', fontweight='bold')

# ---- d: clock-model choice -------------------------------------------------------------------------
ax = fig.add_subplot(gs[1, :]); fig.text(0.6, 0.425, 'd', fontsize=9, fontweight='bold', va='top')
def read_ft(p):
    s = open(p).read(); i = s.index('(', s.index('UTREE')); s = s[i:s.index(';', i) + 1]
    s = re.sub(r'\)\s*\[&95%(?:HPD)?=\{\s*([\d.eE+-]+)\s*,\s*([\d.eE+-]+)\s*\}\]', lambda m: f')HPD_{m.group(1)}_{m.group(2)}', s)
    s = re.sub(r'\[[^\]]*\]', '', s); tr = Phylo.read(io.StringIO(s), 'newick'); out = {}
    for c in tr.get_nonterminals():
        key = frozenset(x.name for x in c.get_terminals()); age = tr.distance(c, c.get_terminals()[0]) * 100
        lo, hi = (float(v) * 100 for v in c.name.split('_')[1:]) if c.name else (age, age); out[key] = (age, lo, hi)
    return out
R = P + 'results/dating_final/runs/'
iln, gbm = read_ft(R + 'post_ML_r1/FigTree.tre'), read_ft(R + 'sens_GBM_r1/FigTree.tre')
keys = [k for k in iln if k in gbm]; xs = np.array([iln[k][0] for k in keys]); ys = np.array([gbm[k][0] for k in keys])
carn = {x for x in tv if tv[x]['nutrition'] == 'carnivorous'}
isc = np.array([len(k & carn) > 0 and len(k - carn) == 0 for k in keys])
ax.fill_between([0, 140], [0, 140], [140, 140], color='#f4f3ef', lw=0, zorder=0)
for k, x, y in zip(keys, xs, ys):
    c = len(k & carn) > 0 and len(k - carn) == 0; cc = S.light(S.C['carnivorous'], 0.55) if c else '#dcdbd6'
    ax.plot([iln[k][1], iln[k][2]], [y, y], color=cc, lw=0.5, zorder=1); ax.plot([x, x], [gbm[k][1], gbm[k][2]], color=cc, lw=0.5, zorder=1)
ax.scatter(xs[~isc], ys[~isc], s=11, color=S.C['autotrophic'], ec='white', lw=0.3, label='other nodes', zorder=3)
ax.scatter(xs[isc], ys[isc], s=16, color=S.C['carnivorous'], ec='white', lw=0.4, label='within carnivorous clades', zorder=4)
ax.plot([0, 140], [0, 140], color=S.INK2, lw=0.7, ls=(0, (3, 2)))
ax.set_xlim(0, 140); ax.set_ylim(0, 140); ax.set_xlabel('Node age, independent rates (Ma)', fontsize=6.8); ax.set_ylabel('Node age, autocorrelated (Ma)', fontsize=6.8)
ax.tick_params(labelsize=6.2); ax.legend(loc='lower right', fontsize=6.2, handletextpad=0.2, borderaxespad=0.3)
ax.text(4, 136, r'$\ln Z_{\mathrm{GBM}}-\ln Z_{\mathrm{ILN}}=48.8\pm0.8$', fontsize=7, va='top', fontweight='bold')
ax.text(4, 121, r'$\ln Z\approx\sum_{i=1}^{8}\frac{w_i}{2}\overline{\ln L}_{\beta_i}$', fontsize=6.6, va='top')
ax.text(4, 103, 'P(GBM | data) > 0.9999\nbars: 95% HPD', fontsize=6.2, va='top', color=S.INK2)
ax.text(96, 70, 'older under ILN', rotation=45, fontsize=6, color=S.INK3, ha='center', va='center')
S.save(fig, 'Fig1_framework_phylogeny'); print('ok', N)
