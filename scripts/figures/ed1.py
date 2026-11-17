"""Extended Data Fig. 1 | Dated 311-taxon plastome phylogeny with tip names, strategy and ndhA–K states (GBM chronogram)."""
import csv, collections, numpy as np, matplotlib.pyplot as plt
from Bio import Phylo
from matplotlib.patches import Rectangle
import figstyle as S
S.setup(); P = '/home/claude/work/pan/'
t = Phylo.read(P + 'phaseC/f311_GBM_t3_chronogram.nwk', 'newick'); tv = {r['tip']: r for r in csv.DictReader(open(P + 'phaseC/phaseC_traits_cov311.csv'))}
tips = t.get_terminals(); N = len(tips); dep = t.depths(); H = max(dep[x] for x in tips)
y = {tp: i for i, tp in enumerate(tips)}
for c in t.get_nonterminals(order='postorder'): y[c] = (y[c.clades[0]] + y[c.clades[-1]]) / 2
import json
name2tip = {x.name: x for x in tips}; lossnodes = []
for nut, k, pp, tl in json.load(open(P + 'phaseC/loss_clades.json')):
    cl = t.common_ancestor([name2tip[x] for x in tl]) if len(tl) > 1 else name2tip[tl[0]]
    assert len(cl.get_terminals()) == k; lossnodes.append((cl, pp))
grp = lambda n: 'carnivorous' if n == 'carnivorous' else 'parasitic' if 'parasit' in n else 'mycoheterotrophic' if 'myco' in n else 'autotrophic'
NC = 4; per = -(-N // NC)
LEG = ('S, strategy: ' + 'carnivorous / parasitic / mycoheterotrophic / autotrophic (strategy colours as in Fig. 1); C, CAM; Aq, submerged aquatic. '
       'ndhA–K: dark functional, grey pseudogene, pale absent, hatched unresolved.\nBlack branches: ndh loss in ≥50% of 200 stochastic maps (posterior frequency above branch). Background: Cretaceous, Paleogene, Neogene (ICS 2024). Time axis in Ma.')
for page in range(2):
    fig = plt.figure(figsize=(S.W2, 247 * S.MM))
    for h in range(2):
        col = page * 2 + h
        ax = fig.add_axes([0.005 + h * 0.5, 0.04, 0.49, 0.9]); ax.axis('off')
        lo, hi = col * per - 0.5, min(N, (col + 1) * per) - 0.5
        S.geobands(ax, H, labels=False, alpha_f=0.88)
        ax.set_xlim(H, -130); ax.set_ylim(hi, lo)
        for c in t.find_clades():
            for ch in c.clades: ax.plot([H - dep[c], H - dep[ch]], [y[ch], y[ch]], color=S.INK2, lw=0.4, solid_capstyle='butt')
            if c.clades: ax.plot([H - dep[c]] * 2, [y[c.clades[0]], y[c.clades[-1]]], color=S.INK2, lw=0.4)
        for cl, pp in lossnodes:
            par = t.get_path(cl); p = par[-2] if len(par) > 1 else t.root
            ax.plot([H - dep[p], H - dep[cl]], [y[cl]] * 2, color=S.LOSS, lw=1.6, solid_capstyle='butt', zorder=3)
            if lo < y[cl] < hi: ax.text((H - dep[p] + H - dep[cl]) / 2, y[cl] - 0.45, f'{pp:.2f}', fontsize=5.5, ha='center', va='bottom', color=S.INK)
        for tp in tips:
            if not (lo < y[tp] < hi): continue
            r = tv[tp.name]; g = grp(r['nutrition']); yy = y[tp]
            ax.add_patch(Rectangle((-1, yy - 0.42), -2.6, 0.84, fc=S.C[g] if g != 'autotrophic' else S.light(S.C['autotrophic'], 0.55), ec='none', clip_on=False))
            if r['CAM'] == '1': ax.text(-5.9, yy, 'C', fontsize=5.5, va='center', ha='center', fontweight='bold')
            if r['aquatic'] == '1': ax.text(-5.9, yy, 'Aq', fontsize=5, va='center', ha='center', fontweight='bold')
            calls = r['ndh_calls']
            if tp.name == 'Byblis_gigantea': calls = ''.join('F' if g_ in 'CGIJK' else '?' for g_ in 'ABCDEFGHIJK')
            for j, ch in enumerate(calls):
                ax.add_patch(Rectangle((-10.5 - j * 2.3, yy - 0.42), -2.05, 0.84, fc=S.GS.get(ch, '#ffffff'), ec=S.INK3 if ch in '-?' else 'none', lw=0.15, hatch='////' if ch == '?' else None, clip_on=False))
            nm = tp.name.replace('_', ' '); nm = nm if len(nm) <= 30 else nm[:29] + '…'
            ax.text(-38, yy, nm, fontsize=6, va='center', ha='left', style='italic', color=S.INK if g != 'autotrophic' else S.INK2)
        ax.text(-10.5 - 1.0, lo - 0.7, 'A', fontsize=6, ha='center', va='bottom', style='italic'); ax.text(-10.5 - 10 * 2.3 - 1.0, lo - 0.7, 'K', fontsize=6, ha='center', va='bottom', style='italic'); ax.text(-10.5 - 5 * 2.3 - 1.0, lo - 0.7, 'ndh', fontsize=6, ha='center', va='bottom', style='italic')
        ax.text(-2.3, lo - 0.7, 'S', fontsize=6, ha='center', va='bottom', fontweight='bold')
        axb = fig.add_axes([0.005 + h * 0.5, 0.03, 0.49 * H / (H + 130), 0.008]); axb.set_xlim(H, 0); axb.set_yticks([])
        for sp in ('left', 'right', 'top'): axb.spines[sp].set_visible(False)
        axb.set_xticks([125, 100, 75, 50, 25, 0]); axb.patch.set_alpha(0); axb.set_xlabel('Ma', labelpad=0)
    fig.text(0.005, 0.995, f'Part {page + 1} of 2 (tips {page * 2 * per + 1}–{min(N, (page * 2 + 2) * per)} of {N}). ' + LEG, fontsize=6, va='top', color=S.INK2, wrap=True)
    S.save(fig, f'ed/ED{page + 1}_full_tree_part{page + 1}')
print('ok', N, len(lossnodes))
