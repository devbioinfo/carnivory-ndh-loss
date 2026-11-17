"""ED7 | Confounders: stratified Pagel tests and phylogenetic logistic regression model comparison."""
import csv, numpy as np, matplotlib.pyplot as plt
import figstyle as S
S.setup(); D = '/home/claude/work/pan/phaseC/'
st = list(csv.DictReader(open(D + 'conf311_GBM_stratified_pagel.csv'))); md = list(csv.DictReader(open(D + 'conf311_GBM_phyloglm_models.csv')))
fig = plt.figure(figsize=(S.W2, 75 * S.MM)); gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1], wspace=0.9, left=0.3, right=0.98, top=0.88, bottom=0.18)
ax = fig.add_subplot(gs[0]); S.panel(ax, 'a', x=-0.75)
for i, r in enumerate(st):
    y = len(st) - 1 - i; lp = -np.log10(float(r['p'])); carn = r['test'].startswith('carnivory')
    if i % 2 == 0: ax.axhspan(y - 0.45, y + 0.45, color='#f6f5f2', lw=0, zorder=0)
    col = S.C['carnivorous'] if carn else S.INK if 'CAM' in r['test'] else S.INK3
    ax.plot([0, lp], [y, y], color=col, lw=1); ax.plot(lp, y, 'o', ms=4, color=col, mec='white', mew=0.4)
    ax.text(lp + 0.7, y + 0.28, f"ratio {float(r['ratio']):.0f}×; n = {r['n_tips']} ({r['n_factor']} with factor)", va='center', fontsize=4.9, color=S.INK2)
ax.axvspan(0, -np.log10(0.05), color='#ecebe6', lw=0); ax.axvline(-np.log10(0.05), color=S.INK2, lw=0.6, ls=(0, (2, 1.5)))
ax.set_yticks(range(len(st))); ax.set_yticklabels([r['test'] for r in st][::-1], fontsize=5.3); ax.set_xlim(0, 13)
ax.set_xlabel(r'$-\log_{10}P$, Pagel dependent vs independent'); ax.set_title('Stratified tests', fontsize=5.8, loc='left', color=S.INK2)
ax = fig.add_subplot(gs[1]); S.panel(ax, 'b', x=-0.9)
lab = {'full': 'strategies + CAM + habitat', 'no_carnivory': 'without carnivory', 'habitat_only': 'CAM + habitat only', 'strategies_only': 'strategies only'}
for i, r in enumerate(md):
    y = len(md) - 1 - i; d = float(r['dAIC']); ax.barh(y, d, height=0.55, color='#8c9fbb' if d > 0 else S.INK, lw=0); ax.text(d + 1, y, f'{d:.1f}', va='center', fontsize=5.2)
ax.set_yticks(range(len(md))); ax.set_yticklabels([lab[r['model']] for r in md][::-1], fontsize=5.4); ax.set_xlim(0, 72)
ax.set_xlabel('ΔAIC vs full model'); ax.set_title('Phylogenetic logistic regression', fontsize=5.8, loc='left', color=S.INK2)
S.save(fig, 'ed/ED8_confounders'); print('ok')
