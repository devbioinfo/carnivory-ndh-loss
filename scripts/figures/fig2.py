"""Figure 2 | Genomic evidence of plastid ndh degradation in carnivorous lineages."""
import csv, json, sys, numpy as np, collections
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import figstyle as S
sys.path.insert(0, '/home/claude/work/pan')
from genbank_utils import read_genbank, _segments
S.setup()
P = '/home/claude/work/pan/'
CALLC = dict(S.GS)
UNRES = dict(fc='white', ec='#a5a39c', hatch='//////', lw=0.3)
import icons as I
byb_orf = {r['gene']: r['call'] for r in csv.DictReader(open(P + 'results/phaseA/byblis_ndh_orfcalls.tsv'), delimiter='\t')}
GENES = ['ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'ndhK']
fig = plt.figure(figsize=(S.W2, 225 * S.MM))
gs = fig.add_gridspec(3, 2, width_ratios=[0.95, 1.4], height_ratios=[1, 0.95, 0.95], hspace=0.5, wspace=0.32,
                      left=0.17, right=0.985, top=0.955, bottom=0.09)

# ---- A: ndh matrix across carnivorous lineages and their relatives ------------------------------
ax = fig.add_subplot(gs[:, 0]); S.panel(ax, 'a', x=-0.55, y=1.0)
tv = list(csv.DictReader(open(P + 'phaseC/phaseC_traits_cov311.csv')))
T = {r['tip']: r for r in tv}
blocks = [('Lentibulariaceae', lambda r: r['family'] == 'Lentibulariaceae'),
          ('Byblidaceae', lambda r: r['family'] == 'Byblidaceae'),
          ('Lamiales relatives', lambda r: r['tip'] in ('Erythranthe_guttata', 'Callitriche_palustris', 'Hippuris_montana', 'Limnophila_sessiliflora', 'Bacopa_monnieri', 'Lindenbergia_philippensis')),
          ('Droseraceae', lambda r: r['family'] == 'Droseraceae'), ('Drosophyllaceae', lambda r: r['family'] == 'Drosophyllaceae'),
          ('Nepenthaceae', lambda r: r['family'] == 'Nepenthaceae'), ('Dioncophyllaceae', lambda r: r['family'] == 'Dioncophyllaceae'),
          ('Ancistrocladaceae', lambda r: r['family'] == 'Ancistrocladaceae'),
          ('Sarraceniaceae', lambda r: r['family'] == 'Sarraceniaceae'), ('Roridulaceae', lambda r: r['family'] == 'Roridulaceae'),
          ('Actinidiaceae', lambda r: r['family'] == 'Actinidiaceae'),
          ('Cephalotaceae', lambda r: r['family'] == 'Cephalotaceae'), ('Brunelliaceae', lambda r: r['family'] == 'Brunelliaceae'),
          ('Bromeliaceae (part)', lambda r: r['tip'].startswith(('Brocchinia', 'Catopsis'))),
          ('Tofieldiaceae', lambda r: r['family'] == 'Tofieldiaceae')]
AUT = S.light(S.C['autotrophic'], 0.55)
rows = []
for name, f in blocks:
    sel = [r for r in tv if f(r)]
    grpd = collections.OrderedDict()
    for r in sorted(sel, key=lambda r: (r['nutrition'] != 'carnivorous', r['tip'])):
        calls = r['ndh_calls'] if r['tip'] != 'Byblis_gigantea' else ''.join('F' if byb_orf.get(g) == 'intact' else '?' for g in GENES)
        key = (r['tip'].split('_')[0], r['nutrition'], calls, r['evidence'])
        grpd.setdefault(key, []).append(r['tip'])
    rows.append(('__gap', None)); rows.append(('__' + name, None))
    for (gen, nut, calls, evd), tl in grpd.items(): rows.append((tl, (gen, nut, calls, evd)))
n = len(rows); y = n
for tip, info in rows:
    y -= 1
    if info is None:
        if tip != '__gap':
            ax.text(-0.25, y + 0.35, tip[2:], ha='right', va='center', fontsize=6.6, fontweight='bold')
            ax.plot([-6.0, 11.6], [y - 0.02, y - 0.02], color=S.INK3, lw=0.4)
        continue
    gen, nut, calls, evd = info
    if nut == 'carnivorous': ax.add_patch(Rectangle((-8.6, y), 20.8, 1, fc=S.light(S.C['carnivorous'], 0.9), ec='none', zorder=0))
    for j, ch in enumerate(calls):
        if ch == '?': ax.add_patch(Rectangle((j + 0.04, y + 0.1), 0.84, 0.8, **UNRES))
        else: ax.add_patch(Rectangle((j, y + 0.06), 0.92, 0.88, fc=CALLC.get(ch, '#fff'), ec=S.INK3 if ch == '-' else 'none', lw=0.2))
    nf = calls.count('F'); ax.text(12.45, y + 0.5, f'{nf}' if '?' not in calls else f'≥{nf}', fontsize=6.2, va='center', ha='right')
    nm = tip[0].replace('_', ' ') if len(tip) == 1 else f'{gen} ({len(tip)} spp.)'
    nm = nm if len(nm) < 27 else nm[:25] + '…'
    ax.text(-0.25, y + 0.5, nm, ha='right', va='center', fontsize=6.2, style='italic', color=S.INK if nut == 'carnivorous' else S.INK2)
    grp_ = 'carnivorous' if nut == 'carnivorous' else 'parasitic' if 'parasit' in nut else 'autotrophic'
    ax.add_patch(Rectangle((11.1, y + 0.1), 0.32, 0.8, fc=S.C[grp_] if grp_ != 'autotrophic' else AUT, ec='none'))
    if 'transcript' in evd: ax.text(12.7, y + 0.5, 'T', fontsize=6.2, va='center', fontweight='bold')
    if 'new assembly' in evd: ax.text(12.7, y + 0.5, 'A', fontsize=6.2, va='center', fontweight='bold')
for j, g in enumerate(GENES): ax.text(j + 0.46, n + 0.2, g, rotation=90, ha='center', va='bottom', fontsize=6.6, style='italic')
ax.text(12.1, n + 0.2, 'n F', rotation=90, ha='center', va='bottom', fontsize=6.4)
ax.set_xlim(-0.2, 12.4); ax.set_ylim(0, n + 3.0); ax.axis('off')
lg = [(CALLC['F'], 'functional'), (CALLC['P'], 'pseudogene'), (CALLC['-'], 'absent'), (None, 'unresolved'), (S.C['carnivorous'], 'carnivorous'), (S.C['parasitic'], 'parasitic'), (AUT, 'autotrophic')]
for k, (c, l) in enumerate(lg):
    xx_, yy_ = -8.4 + (k % 4) * 5.3, -2.6 - (k // 4) * 1.9
    if c is None: ax.add_patch(Rectangle((xx_, yy_), 0.7, 1.1, clip_on=False, **UNRES))
    else: ax.add_patch(Rectangle((xx_, yy_), 0.7, 1.1, fc=c, ec=S.INK3 if c == CALLC['-'] else 'none', lw=0.3, clip_on=False))
    ax.text(xx_ + 0.95, yy_ + 0.55, l, fontsize=6.4, va='center')
ax.text(-8.4, -6.6, 'Rows: species with identical ndh patterns within a genus are merged\n(n spp.; all 311 taxa in Extended Data Fig. 1). n F, functional genes.\nT, transcript-based (Byblis: unassembled or truncated genes shown\nas unresolved); A, new assembly. Shaded rows: carnivorous taxa.', fontsize=6.2, color=S.INK2, va='top', linespacing=1.25)

# ---- B: Triantha plastome, the two ndh regions ------------------------------------------------
rec = next(read_genbank(P + 'genbank/TRIOCC_ERR8268345.gb')); L = len(rec['seq'])
dep = np.load(P + 'figures/triantha_depth.npy')
calls = dict(zip(GENES, T['Triantha_occidentalis']['ndh_calls']))
supp = {r['gene']: r for r in csv.DictReader(open(P + 'results/phaseA/triantha_ndh_readsupport.tsv'), delimiter='\t')}
gsb = gs[0, 1].subgridspec(2, 2, width_ratios=[13, 3.2], height_ratios=[1, 1.15], wspace=0.08, hspace=0.05)
regions = [(53.6, 66.6), (112.6, 115.8)]
feats = []
for f in rec['features']:
    if f['type'] != 'CDS': continue
    g = f['qualifiers']['gene'][0]; segs = _segments(f['location'])
    feats.append((g, min(a for a, b, s_ in segs) / 1000, max(b for a, b, s_ in segs) / 1000, segs[0][2]))
for k, (x0, x1) in enumerate(regions):
    axd = fig.add_subplot(gsb[0, k]); axg = fig.add_subplot(gsb[1, k], sharex=axd)
    if k == 0: S.panel(axd, 'b', x=-0.12, y=1.08)
    i0, i1 = int(x0 * 1000), int(x1 * 1000); xx = np.arange(i0, i1) / 1000
    axd.fill_between(xx, 0, dep[i0:i1], color='#c9d6e8', lw=0); axd.plot(xx, dep[i0:i1], color='#6d87ad', lw=0.35)
    axd.axhline(160, color=S.INK3, lw=0.5, ls=(0, (2, 1.5))); axd.text(x1, 175, '160×', fontsize=4.6, color=S.INK3, ha='right', va='bottom') if k == 1 else None; axd.set_ylim(0, 700); axd.set_xlim(x0, x1)
    axd.tick_params(labelbottom=False, bottom=False); axd.spines['bottom'].set_visible(False)
    if k == 0: axd.set_ylabel('Depth', labelpad=1)
    else: axd.tick_params(labelleft=False); axd.spines['left'].set_visible(False); axd.tick_params(left=False)
    axg.set_ylim(-1.6, 1.4); axg.axhline(0, color=S.INK3, lw=0.4)
    for g, a0, a1, st in feats:
        if a1 < x0 or a0 > x1: continue
        isn = g.startswith('ndh'); yy = 0.25 if st == 1 else -0.25
        col = (CALLC['F'] if calls.get(g) == 'F' else CALLC['P']) if isn else '#dcdbd6'
        if isn and calls.get(g) != 'F':
            for axx in (axd, axg): axx.axvspan(a0, a1, color=S.light(CALLC['P'], 0.6), lw=0, zorder=0)
        axg.add_patch(Rectangle((a0, yy - 0.22), a1 - a0, 0.44, fc=col, ec=S.INK2 if (isn and calls.get(g) != 'F') else 'none', lw=0.3))
        if not isn: continue
        nshift = {'ndhD': -0.55, 'ndhE': 0.0, 'ndhI': 0.45, 'ndhC': -0.25, 'ndhK': -0.32, 'ndhJ': 0.42}.get(g, 0)
        lev = {'ndhE': 0.42, 'ndhK': 0.0, 'ndhC': 0.42}.get(g, 0)
        axg.text((a0 + a1) / 2 + nshift, yy + (0.35 + lev if st == 1 else -0.35 - lev), g, ha='center', va='bottom' if st == 1 else 'top',
                 fontsize=6.2, style='italic', color=S.INK, fontweight='bold' if calls.get(g) == 'F' else 'normal')
    if k == 0:
        xg = (61597 + 61776) / 2000
        axg.plot(xg, 0, marker='x', ms=4, mew=1.0, color=S.INK, zorder=5)
        axg.annotate('ndhG deleted (read-confirmed)', (xg, 0), xytext=(xg - 1.6, -1.45), fontsize=5, color=S.INK, ha='center',
                     arrowprops=dict(arrowstyle='-', lw=0.5, color=S.INK))
    axg.set_yticks([]); axg.spines['left'].set_visible(False)
    if k == 0:
        axg.set_xlabel('Position (kb)', labelpad=1)
        axg.text(0.0, -0.42, 'Top: read depth (dashed, 160×). Bottom: CDS on + / − strands;\ngrey shading: read-confirmed ndh lesions; bold: functional genes.', transform=axg.transAxes, fontsize=6.2, color=S.INK2, va='top')
    axg.tick_params(axis='x', labelsize=5.8)
md = np.median([float(r['mean_depth']) for r in supp.values() if r['gene'] in ('ndhC','ndhD','ndhF','ndhH','ndhI','ndhA','ndhE','ndhJ','ndhK')])
fig.axes[1].text(0.0, 1.08, f'Triantha occidentalis (new assembly, 154.9 kb): lesions in 6 of 11 ndh\ngenes; median depth {md:.0f}× at ndh loci; no site with <80% read agreement', transform=fig.axes[1].transAxes, fontsize=5.6, color=S.INK2, va='bottom')

# ---- C: Roridula ------------------------------------------------------------------------------
gsc = gs[1, 1].subgridspec(1, 2, width_ratios=[1.25, 1], wspace=0.45)
ax = fig.add_subplot(gsc[0]); S.panel(ax, 'c', x=-0.2)
gi = [l.rstrip('\n').split('\t') for l in open(P + 'results/Roridula_verification/gene_identity_fu_vs_chang.tsv') if '\t' in l][1:]
lsc = [float(r[6]) for r in gi if r[1] == 'LSC']; ssc = [(r[0], float(r[6])) for r in gi if r[1] == 'SSC']
rng = np.random.default_rng(1)
ax.scatter(rng.uniform(-0.15, 0.15, len(lsc)), lsc, s=11, color=S.light(S.INK3, 0.3), ec='white', lw=0.3)
ax.axhspan(60, 99, xmin=0.5, xmax=1, color='#f6f5f2', lw=0, zorder=0)
ax.scatter(1 + rng.uniform(-0.15, 0.15, len(ssc)), [v for _, v in ssc], s=12, color=S.INK, ec='white', lw=0.3)
ax.set_xticks([0, 1]); ax.set_xticklabels(['LSC genes', 'SSC genes']); ax.set_xlim(-0.5, 1.5)
ax.set_ylabel('Identity, NC_069066 vs\nPP273159.1 (%)', labelpad=1); ax.set_ylim(60, 101)
ax.text(0.86, min(v for _, v in ssc), 'ycf1', fontsize=6.2, style='italic', va='center', ha='right')
ax.set_title('Roridula gorgonias: 3 SSC ndh genes\nabsent from NC_069066', fontsize=6.4, loc='left', color=S.INK)
ax = fig.add_subplot(gsc[1])
k1, k2 = json.load(open(P + 'results/Roridula_verification/rc_kew.json')), json.load(open(P + 'results/Roridula_verification/rc_wgs2020.json'))
J = 'junction_spanning_reads(>=99%,±40bp)'
cats = ['LSC/IR', 'SSC/IRb', 'IRa/SSC']
def g(d, pre):
    j = d[pre + ':' + J]; return [sum(v for k, v in j.items() if k.startswith('LSC')), j.get('SSC/IRb', 0), j.get('IRa/SSC', 0)]
xx = np.arange(3); wbar = 0.36
for i, (lab_, pre, col) in enumerate([('NC_069066', 'fu', 'white'), ('PP273159.1', 'chang', S.GS['F'])]):
    v = np.array(g(k1, pre)) + np.array(g(k2, pre))
    ax.bar(xx + (i - 0.5) * wbar, v, width=wbar * 0.92, color=col, ec=S.INK2 if i == 0 else 'none', hatch='////' if i == 0 else None, lw=0.5, label=lab_)
    for x_, vv in zip(xx + (i - 0.5) * wbar, v):
        if vv == 0: ax.text(x_, 8, '0', ha='center', fontsize=5.5, color=S.INK, fontweight='bold')
ax.set_xticks(xx); ax.set_xticklabels(cats, fontsize=5.8); ax.set_ylabel('Junction-spanning reads\n(≥99% identity, 2 datasets)')
ax.legend(fontsize=5.3, loc='upper left', handlelength=1, bbox_to_anchor=(0.0, 1.12), ncol=2, columnspacing=0.6)

# ---- D: Byblis transcripts ---------------------------------------------------------------------
ax = fig.add_subplot(gs[2, 1]); S.panel(ax, 'd', x=-0.06)
seg = list(csv.DictReader(open(P + 'results/phaseA/byblis_ndh_segments.tsv'), delimiter='\t'))
reflen = {}
rec2 = next(read_genbank(P + 'genbank/PP316519.1.gb'))
from genbank_utils import cds_records
for gname, nt, p in cds_records(rec2):
    if gname.startswith('ndh') and gname not in reflen: reflen[gname] = len(nt)
byb = {r['gene']: r['call'] for r in csv.DictReader(open(P + 'results/phaseA/byblis_ndh_orfcalls.tsv'), delimiter='\t')}
for j, gname in enumerate(GENES):
    yy = len(GENES) - 1 - j; Lr = reflen.get(gname, 1)
    ax.add_patch(Rectangle((0, yy - 0.08), 1, 0.16, fc='#ecebe6', ec='none'))
    for r in seg:
        if r['gene'] != gname: continue
        a, b = int(r['ref_from']) / Lr, int(r['ref_to']) / Lr; ok = True
        cov = float(r['transcript_kmer_cov'])
        st0 = byb.get(gname, 'unresolved')
        fc_ = CALLC['F'] if st0 == 'intact' else '#8a8984'
        ax.add_patch(Rectangle((a, yy - 0.3), b - a, 0.6, fc=fc_, ec='white', lw=0.3, alpha=min(1, 0.3 + np.log10(cov + 1) / 2.5)))
    st = byb.get(gname, 'unresolved')
    if st != 'intact': ax.add_patch(Rectangle((-0.004, yy - 0.36), 1.008, 0.72, fc='none', ec='#a5a39c', hatch='//////', lw=0.3, zorder=0))
    ax.text(1.03, yy, {'intact': 'intact ORF', 'truncated': 'truncated*'}.get(st, 'unresolved'), fontsize=5.4, va='center',
            color=CALLC['F'] if st == 'intact' else S.INK2, fontweight='bold' if st == 'intact' else 'normal', style='normal' if st == 'intact' else 'italic')
ax.set_yticks(range(len(GENES))); ax.set_yticklabels(GENES[::-1], style='italic', fontsize=6.4); ax.set_ylim(-0.7, len(GENES) - 0.4)
ax.set_xlim(0, 1.3); ax.set_xticks([0, 0.25, 0.5, 0.75, 1]); ax.set_xticklabels(['0', '25', '50', '75', '100'])
ax.set_xlabel('Position along reference CDS (%)'); ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
ax.set_title('Byblis gigantea plastid transcripts (2 RNA-seq runs). Dark: segments of intact ORFs.\nGrey on hatched track: truncated* or unassembled CDS (not evidence of loss).\nOpacity scales with log k-mer coverage. *Not separable from fragmented assembly', loc='left', fontsize=5.6, color=S.INK2)
S.save(fig, 'Fig2_genomic_evidence')
print('ok', n)
