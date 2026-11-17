"""ED6 | Stochastic mapping: number of losses per map by regime, and posterior support of each independent loss."""
import csv, json, collections, numpy as np, matplotlib.pyplot as plt
import figstyle as S
S.setup(); D = '/home/claude/work/pan/phaseC/'
ev = list(csv.DictReader(open(D + 'f311_GBM_t3_loss_events.csv'))); sims = sorted({int(e['sim']) for e in ev})
G = {'autotrophic': 'autotrophic', 'carnivorous': 'carnivorous', 'hemiparasitic': 'parasitic', 'holoparasitic': 'parasitic', 'full_mycoheterotroph': 'mycoheterotrophic', 'partial_mycoheterotroph': 'mycoheterotrophic'}
cnt = collections.defaultdict(lambda: collections.Counter())
for e in ev: cnt[G[e['regime']]][int(e['sim'])] += 1
fig = plt.figure(figsize=(S.W2, 120 * S.MM)); gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=0.55, left=0.08, right=0.98, top=0.93, bottom=0.1)
ga = gs[0].subgridspec(4, 1, hspace=0.35)
for k, g in enumerate(['autotrophic', 'carnivorous', 'parasitic', 'mycoheterotrophic']):
    ax = fig.add_subplot(ga[k]); v = np.array([cnt[g][s_] for s_ in sims]); vals, n = np.unique(v, return_counts=True)
    ax.bar(vals, n / len(sims), width=0.8, color=S.C[g], lw=0); ax.set_xlim(-0.6, 16); ax.set_ylim(0, 0.9)
    ax.text(15.8, 0.8, f'{g.capitalize()}: median {np.median(v):.0f}', ha='right', va='top', fontsize=5.6)
    ax.yaxis.grid(True, color=S.GRID, lw=0.4); ax.set_axisbelow(True); ax.set_yticks([0, 0.4, 0.8])
    if k < 3: ax.tick_params(labelbottom=False)
    if k == 0: S.panel(ax, 'a', x=-0.2, y=1.1)
    if k == 1: ax.set_ylabel('Fraction of maps', y=-0.1)
ax.set_xlabel(f'ndh losses per stochastic map (n = {len(sims)} maps)')
ax = fig.add_subplot(gs[1]); S.panel(ax, 'b', x=-0.45)
L = json.load(open(D + 'loss_clades.json')); L.sort(key=lambda r: (['carnivorous', 'parasitic', 'mycoheterotrophic', 'CAM', 'aquatic', 'autotrophic'].index(r[0]), -r[2]))
for i, (nut, k, p, tl) in enumerate(L):
    y = len(L) - 1 - i; col = S.C.get(nut, S.INK3) if nut in ('carnivorous', 'parasitic', 'mycoheterotrophic') else S.INK3
    ax.barh(y, p, height=0.7, color=col, lw=0); t0 = tl[0].replace('_', ' ').split()
    ax.text(-0.02, y, (t0[0] if k > 1 else f'{t0[0][0]}. {t0[1]}') + (f' ({k})' if k > 1 else '') + ('' if nut in ('carnivorous', 'parasitic', 'mycoheterotrophic', 'autotrophic') else f' [{nut}]'), ha='right', va='center', fontsize=5.0, style='italic')
    ax.text(p + 0.01, y, f'{p:.2f}', va='center', fontsize=4.8, color=S.INK2)
ax.set_yticks([]); ax.spines['left'].set_visible(False); ax.set_xlim(0, 1.12); ax.axvline(0.5, color=S.INK3, lw=0.5, ls=(0, (2, 1.5)))
ax.set_xlabel('Posterior frequency of the loss event'); ax.set_ylim(-0.6, len(L) - 0.4)
ax.set_title(f'{len(L)} independent losses present in ≥50% of maps; clade size in brackets', fontsize=5.8, loc='left', color=S.INK2)
S.save(fig, 'ed/ED7_stochastic_maps'); print('ok')
