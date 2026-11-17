"""ED2 | MCMCTree convergence: node ages from replicate chains, and effective priors vs posteriors at calibrations."""
import re, io, csv, numpy as np, matplotlib.pyplot as plt
from Bio import Phylo
import figstyle as S
S.setup(); P = '/home/claude/work/pan/results/dating_final/'
def read_ft(p):
    s = open(p).read(); i = s.index('(', s.index('UTREE')); s = s[i:s.index(';', i) + 1]
    s = re.sub(r'\[[^\]]*\]', '', s); tr = Phylo.read(io.StringIO(s), 'newick'); out = {}
    for c in tr.get_nonterminals(): out[frozenset(x.name for x in c.get_terminals())] = tr.distance(c, c.get_terminals()[0]) * 100
    return out
conv = {r['analysis']: r for r in csv.DictReader(open(P + 'results/convergence.csv'))}
runs = [('post_ML', 'ML topology (ILN)'), ('post_SR', 'SR topology (ILN)'), ('post_H4', 'H4 topology (ILN)'), ('post_H5', 'H5 topology (ILN)'),
        ('post_H6', 'H6 topology (ILN)'), ('sens_GBM', 'ML topology (GBM, preferred)')]
fig = plt.figure(figsize=(S.W2, 150 * S.MM))
gs = fig.add_gridspec(2, 4, hspace=0.5, wspace=0.6, left=0.06, right=0.98, top=0.95, bottom=0.08)
for k, (run, lab) in enumerate(runs):
    ax = fig.add_subplot(gs[k // 3, k % 3]); a, b = read_ft(P + f'runs/{run}_r1/FigTree.tre'), read_ft(P + f'runs/{run}_r2/FigTree.tre')
    ks = [x for x in a if x in b]; xs = np.array([a[x] for x in ks]); ys = np.array([b[x] for x in ks])
    ax.plot([0, 145], [0, 145], color=S.INK3, lw=0.6, ls=(0, (2, 1.5))); ax.scatter(xs, ys, s=7, color='#3e5373' if run != 'sens_GBM' else S.INK, ec='white', lw=0.25)
    ax.set_xlim(0, 145); ax.set_ylim(0, 145); ax.set_aspect('equal'); ax.set_title(lab, fontsize=6.2, loc='left', color=S.INK2)
    c = conv.get(run)
    if c: ax.text(4, 140, f"r = {float(c['corr_node_means']):.5f}\nmax diff {c['max_rel_diff_%_nodes_>5Myr']}%\nESS≥200: {c['pct_nodes_ESS>=200_combined']}% nodes", fontsize=5.0, va='top')
    if k % 3 == 0: ax.set_ylabel('Chain 2 age (Ma)')
    if k // 3 == 1: ax.set_xlabel('Chain 1 age (Ma)')
    if k == 0: S.panel(ax, 'a', x=-0.3)
# b: effective prior vs posterior at calibrations
W = {r['node']: r for r in csv.DictReader(open(P + 'results/node_ages_wide.csv'))}
def parse(s): m = re.match(r'([\d.]+) \(([\d.]+)-([\d.]+)\)', s); return tuple(map(float, m.groups()))
ax = fig.add_subplot(gs[:, 3]); S.panel(ax, 'b', x=-0.55)
nodes = [n for n in W if re.search(r'\((C\d|X\d)\)', n)]
lab = {'root (C1)': 'C1 root', 'crown mesangiosperms (C2)': 'C2 mesangiosperms', 'crown core Ericales (C3)': 'C3 core Ericales', 'stem Roridula (C4)': 'C4 stem Roridula',
       'stem Aldrovanda (C5)': 'C5 stem Aldrovanda'}
for i, n in enumerate(nodes):
    y = len(nodes) - 1 - i
    if i % 2 == 0: ax.axhspan(y - 0.45, y + 0.45, color='#f6f5f2', lw=0, zorder=0)
    for k, (col_, c, mk) in enumerate([('prior_ML', '#b9c6da', 'o'), ('post_ML_r1', '#3e5373', 's'), ('sens_GBM_r1', S.INK, 'D')]):
        m, lo, hi = parse(W[n][col_]); yy = y + 0.22 - k * 0.22
        ax.plot([lo, hi], [yy, yy], color=c, lw=1.2 if k == 0 else 0.8); ax.plot(m, yy, mk, ms=3, color=c, mec='white', mew=0.3)
ax.set_yticks(range(len(nodes))); ax.set_yticklabels([lab.get(n, n.replace(' (X', '\n(X')) for n in nodes][::-1], fontsize=4.8)
ax.set_xlim(150, 0); ax.set_xlabel('Age (Ma)')
for c, l, mk in [('#b9c6da', 'effective prior', 'o'), ('#3e5373', 'posterior, ILN', 's'), (S.INK, 'posterior, GBM', 'D')]: ax.plot([], [], mk + '-', color=c, ms=3, lw=0.8, label=l)
ax.legend(loc='upper center', bbox_to_anchor=(0.4, -0.07), fontsize=5, ncol=1)
ax.set_title('Calibrations (C) and cross-validation\nfossils (X), mean and 95% HPD', fontsize=5.8, loc='left', color=S.INK2)
S.save(fig, 'ed/ED3_dating_convergence'); print('ok', len(nodes))
