"""ED3 | Sensitivity of node ages to topology, clock model, calibration and data partition."""
import re, csv, numpy as np, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import figstyle as S
S.setup(); P = '/home/claude/work/pan/results/dating_final/results/'
rows = list(csv.DictReader(open(P + 'node_ages_wide.csv')))
cols = [('post_ML_r1', 'ML'), ('post_SR_r1', 'SR'), ('post_H4_r1', 'H4\n(rejected)'), ('post_H5_r1', 'H5'), ('post_H6_r1', 'H6'),
        ('sens_minima_r1', 'minima\nonly'), ('sens_softroot_r1', 'soft\nroot'), ('sens_pos12_r1', 'codon\n1+2'), ('sens_GBM_r1', 'GBM\n(preferred)')]
def parse(s):
    m = re.match(r'([\d.]+) \(([\d.]+)-([\d.]+)\)', s or ''); return tuple(map(float, m.groups())) if m else (np.nan,) * 3
rows = [r for r in rows if not r['node'].startswith(('root', 'crown mesangiosperms', 'crown Brocchinia reducta'))]
M = np.array([[parse(r[c])[0] for c, _ in cols] for r in rows]); ref = M[:, 0:1]
D = np.log2(M / ref)
fig, ax = plt.subplots(figsize=(S.W2, 170 * S.MM)); fig.subplots_adjust(left=0.27, right=0.9, top=0.88, bottom=0.04)
cm = LinearSegmentedColormap.from_list('d', ['#3e5373', '#f0efec', '#a8714a'])
im = ax.imshow(D, cmap=cm, vmin=-2, vmax=2, aspect='auto')
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if np.isfinite(M[i, j]): ax.text(j, i, f'{M[i, j]:.0f}', ha='center', va='center', fontsize=4.8, color=S.INK if abs(D[i, j]) < 1.3 else 'white')
ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[1] for c in cols], fontsize=5.4); ax.xaxis.tick_top()
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r['node'] for r in rows], fontsize=5.2)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.tick_params(length=0); ax.axvline(4.5, color='white', lw=2); ax.axvline(7.5, color='white', lw=2)
ax.text(2, -2.6, 'Topology (independent rates)', ha='center', fontsize=5.6, fontweight='bold', color=S.INK2)
ax.text(6, -2.6, 'Prior / data (ML, ILN)', ha='center', fontsize=5.6, fontweight='bold', color=S.INK2); ax.text(8, -2.6, 'Clock', ha='center', fontsize=5.6, fontweight='bold', color=S.INK2)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.outline.set_visible(False); cb.set_ticks([-2, -1, 0, 1, 2]); cb.set_ticklabels(['×¼', '×½', '1', '×2', '×4'])
cb.set_label('Age relative to ML, ILN (log$_2$ scale)', fontsize=5.6); cb.ax.tick_params(labelsize=5)
ax.text(-0.5, len(rows) + 0.4, 'Cell text: posterior mean age (Ma). Colour: ratio to the ML-topology, independent-rates analysis. Cross-validation fossils marked (X).', fontsize=5, color=S.INK2, va='top')
S.save(fig, 'ed/ED4_dating_sensitivity'); print('ok', len(rows))
