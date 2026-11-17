"""ED5 | Sensitivity of the rate estimates to model, loss threshold, taxon set and clock; Mk model comparison."""
import csv, re, numpy as np, matplotlib.pyplot as plt
import figstyle as S
S.setup(); D = '/home/claude/work/pan/phaseC/'
sets = [('f311_GBM_t3', '311 · GBM · ML · ≥3 (main)'), ('final_ML_GBM_t3', '295 · GBM · ML · ≥3'), ('final_ML_t3', '295 · ILN · ML · ≥3'),
        ('final_SR_t3', '295 · ILN · SR · ≥3'), ('final_ML_t1', '295 · ILN · ML · ≥1'), ('final_SR_t1', '295 · ILN · SR · ≥1')]
def p3(s): m = re.match(r'([\d.]+) \(([\d.]+)-([\d.]+)\)', s); return tuple(map(float, m.groups()))
fig = plt.figure(figsize=(S.W2, 135 * S.MM)); gs = fig.add_gridspec(2, 2, wspace=0.75, hspace=0.75, left=0.17, right=0.98, top=0.93, bottom=0.1, width_ratios=[1.2, 1], height_ratios=[1, 0.62])
ax = fig.add_subplot(gs[0, 0]); S.panel(ax, 'a', x=-0.4)
for j, (k, lab) in enumerate(sets):
    y = len(sets) - 1 - j
    try: R = {r['regime']: r for r in csv.DictReader(open(D + k + '_rates_summary.csv'))}
    except FileNotFoundError: continue
    if j % 2 == 0: ax.axhspan(y - 0.45, y + 0.45, color='#f6f5f2', lw=0, zorder=0)
    for dk, (reg, col, mk) in enumerate([('carnivorous', S.C['carnivorous'], 'o'), ('hemiparasitic', S.C['parasitic'], 's')]):
        m, lo, hi = p3(R[reg]['ratio_vs_autotrophic']); yy = y + 0.12 - dk * 0.24
        ax.plot([lo, hi], [yy, yy], color=col, lw=0.9); ax.plot(m, yy, mk, ms=3.8, color=col, mec='white', mew=0.4)
ax.axvline(1, color=S.INK2, lw=0.6); ax.set_xscale('log'); ax.set_xlim(0.8, 100)
ax.set_yticks(range(len(sets))); ax.set_yticklabels([s[1] for s in sets][::-1], fontsize=5.6)
ax.set_xlabel('Loss rate relative to autotrophs (median, 95% interval)')
ax.plot([], [], 'o-', color=S.C['carnivorous'], ms=3.5, lw=0.9, label='carnivorous'); ax.plot([], [], 's-', color=S.C['parasitic'], ms=3.5, lw=0.9, label='hemiparasitic')
ax.legend(loc='lower right', fontsize=5.4); ax.set_title('Taxa · clock · topology · loss threshold (genes)', fontsize=5.8, loc='left', color=S.INK2)
ax = fig.add_subplot(gs[0, 1]); S.panel(ax, 'b', x=-0.55)
for j, (k, lab) in enumerate(sets):
    y = len(sets) - 1 - j
    try: R = {r['model']: r for r in csv.DictReader(open(D + k + '_ndh_model_comparison.csv'))}
    except FileNotFoundError: continue
    if j % 2 == 0: ax.axhspan(y - 0.45, y + 0.45, color='#f6f5f2', lw=0, zorder=0)
    d = float(R['irreversible loss']['AIC']) - float(R['reversible (ARD)']['AIC'])
    ax.barh(y, d, height=0.5, color='#8c9fbb'); ax.text(d + 0.3, y, f'{d:.1f}', va='center', fontsize=5.2)
ax.axvline(2, color=S.INK3, lw=0.5, ls=(0, (2, 1.5))); ax.set_yticks(range(len(sets))); ax.set_yticklabels([s[1] for s in sets][::-1], fontsize=5.6)
ax.set_xlabel('ΔAIC, irreversible − reversible Mk'); ax.set_title('ndh character model (positive favours reversible)', fontsize=5.8, loc='left', color=S.INK2)
ax = fig.add_subplot(gs[1, :]); S.panel(ax, 'c', x=-0.115)
H = list(csv.DictReader(open(D + 'hrm_model_comparison.csv')))
lab = {'Mk irreversible': 'one rate, loss only (k = 1)', 'Mk reversible (ARD)': 'one loss + one regain rate (k = 2)',
       'HRM irreversible, 2 rates': 'two hidden rate classes, loss only (k = 3)', 'HRM reversible, 2 rates': 'two hidden classes, regain allowed (k = 5)'}
for i, r in enumerate(H):
    y = len(H) - 1 - i; d = float(r['dAIC']); best = d == 0
    if i % 2 == 0: ax.axhspan(y - 0.45, y + 0.45, color='#f6f5f2', lw=0, zorder=0)
    ax.barh(y, max(d, 0.3), height=0.55, color=S.INK if best else '#8c9fbb', lw=0)
    ax.text(max(d, 0.3) + 0.4, y, f"ΔAIC {d:.1f};  w = {float(r['w']):.3f}" + ('  (best)' if best else ''), va='center', fontsize=6.2)
ax.set_yticks(range(len(H))); ax.set_yticklabels([lab[r['model']] for r in H][::-1]); ax.set_xlim(0, 36)
ax.set_xlabel('ΔAIC relative to best model (311 taxa, GBM chronogram, ≥3-gene threshold)')
ax.set_title('Rate heterogeneity removes support for regain: estimated regain rates = 0; fast-class loss rate 0.031 Myr$^{-1}$ (carnivore rate 0.029)', fontsize=6.4, loc='left', color=S.INK2)
S.save(fig, 'ed/ED6_model_threshold'); print('ok')
