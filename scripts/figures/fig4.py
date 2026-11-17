"""Figure 4 | Statistical robustness of the heterotrophy-ndh association and of the timescale."""
import csv, numpy as np
import matplotlib.pyplot as plt
import figstyle as S
S.setup()
D = '/home/claude/work/pan/phaseC/'
fig = plt.figure(figsize=(S.W2, 150 * S.MM))
gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.62, left=0.15, right=0.98, top=0.95, bottom=0.07, width_ratios=[1.05, 1])

# ---- A: Pagel tests across analytical settings ----------------------------------------------
ax = fig.add_subplot(gs[0, 0]); S.panel(ax, 'a', x=-0.33)
sets = [('f311_GBM_t3', '311 taxa · GBM · ML · ≥3 genes (main)'), ('final_ML_GBM_t3', '295 · GBM · ML · ≥3'),
        ('final_ML_t3', '295 · ILN · ML · ≥3'), ('final_SR_t3', '295 · ILN · SR · ≥3'),
        ('final_ML_t1', '295 · ILN · ML · ≥1 gene'), ('final_SR_t1', '295 · ILN · SR · ≥1 gene')]
contr = [('carnivory vs autotrophy', 'Carnivory', S.C['carnivorous']), ('parasitism vs autotrophy', 'Parasitism', S.C['parasitic']),
         ('mycoheterotrophy vs autotrophy', 'Mycoheterotrophy', S.C['mycoheterotrophic'])]
mk = ['o', 's', '^']
for j, (key, lab) in enumerate(sets):
    rows = {r['contrast']: r for r in csv.DictReader(open(D + key + '_pagel.csv'))}
    y = len(sets) - 1 - j
    ax.axhspan(y - 0.45, y + 0.45, color=S.light(S.C['carnivorous'], 0.9) if j == 0 else ('#f6f5f2' if j % 2 == 0 else 'white'), lw=0, zorder=0)
    for k, (c, cl, col) in enumerate(contr):
        p = float(rows[c]['p_dep_ndh']); ax.plot(-np.log10(p), y + (1 - k) * 0.22, mk[k], color=col, ms=4.4, mec='white', mew=0.5, zorder=3)
ax.axvspan(0, -np.log10(0.05), color='#ecebe6', lw=0, zorder=0.5); ax.axvline(-np.log10(0.05), color=S.INK2, lw=0.6, ls=(0, (2, 1.5)))
ax.text(-np.log10(0.05) + 0.15, -0.6, 'P = 0.05', fontsize=5.4, color=S.INK2)
ax.text(0.1, len(sets) - 0.75, 'n.s.', fontsize=5, color=S.INK3)
ax.set_yticks(range(len(sets))); ax.set_yticklabels([l for _, l in sets][::-1], fontsize=6)
ax.set_xlabel(r'$-\log_{10}P$ (Pagel dependent vs independent model)', labelpad=1); ax.set_xlim(0, 11); ax.set_ylim(-0.7, len(sets) - 0.6)
ax.xaxis.grid(True, color='white', lw=0.5); ax.set_axisbelow(True)
for k, (c, cl, col) in enumerate(contr): ax.plot([], [], mk[k], color=col, ms=4, label=cl)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=3, handletextpad=0.2, columnspacing=1.0)
ax.set_title('Taxon set · clock · topology · loss definition', loc='left', fontsize=6.3, color=S.INK2)

# ---- B: leave-one-lineage-out -----------------------------------------------------------------
ax = fig.add_subplot(gs[0, 1]); S.panel(ax, 'b', x=-0.42)
jk = list(csv.DictReader(open(D + 'jk311_A.csv'))) + list(csv.DictReader(open(D + 'jk311_B.csv')))
names = {'none': 'None (all carnivores)', 'Lentibulariaceae': 'Lentibulariaceae', 'Utricularia': 'Utricularia', 'Droseraceae': 'Droseraceae',
         'Drosophyllum': 'Drosophyllum', 'Nepenthes': 'Nepenthes', 'Sarraceniaceae': 'Sarraceniaceae', 'Triantha': 'Triantha', 'Brocchinia': 'Brocchinia'}
ital = {'Utricularia', 'Drosophyllum', 'Nepenthes', 'Triantha', 'Brocchinia'}
for i, r in enumerate(jk):
    y = len(jk) - 1 - i; lp = -np.log10(float(r['p']))
    ax.axhspan(y - 0.45, y + 0.45, color='#f6f5f2' if i % 2 == 0 else 'white', lw=0, zorder=0)
    ax.plot([0, lp], [y, y], color=S.C['carnivorous'] if r['dropped'] != 'none' else S.INK, lw=1.0, solid_capstyle='butt')
    ax.plot(lp, y, 'o', ms=3.8, color=S.C['carnivorous'] if r['dropped'] != 'none' else S.INK, mec='white', mew=0.4)
    ax.text(lp + 0.25, y, f"{float(r['ratio']):.0f}×  (n = {r['n_carn']})", va='center', fontsize=5.7, color=S.INK2)
ax.set_yticks(range(len(jk))); ax.set_yticklabels([names[r['dropped']] for r in jk][::-1], fontsize=6)
for t, r in zip(ax.get_yticklabels(), jk[::-1]):
    if r['dropped'] in ital: t.set_style('italic')
ax.axvspan(0, -np.log10(0.05), color='#ecebe6', lw=0, zorder=0.5); ax.axvline(-np.log10(0.05), color=S.INK2, lw=0.6, ls=(0, (2, 1.5)))
rr = [float(r['ratio']) for r in jk]; ax.text(14.4, -0.75, f'all P < 4×10$^{{-6}}$; ratio range {min(rr):.0f}–{max(rr):.0f}×', ha='right', fontsize=5.3, color=S.INK2)
ax.set_ylim(-1.0, len(jk) - 0.5); ax.set_xlim(0, 14.5); ax.set_xlabel(r'$-\log_{10}P$, carnivory vs autotrophy')
ax.set_title('Lineage removed (loss-rate ratio, carnivores left)', loc='left', fontsize=6.3, color=S.INK2)
ax.xaxis.grid(True, color=S.GRID, lw=0.4); ax.set_axisbelow(True)

# ---- C: phylogenetic logistic regression + stratified tests ----------------------------------
ax = fig.add_subplot(gs[1, 0]); S.panel(ax, 'c', x=-0.33)
co = {r['term']: r for r in csv.DictReader(open(D + 'conf311_GBM_phyloglm_coefficients.csv'))}
terms = [('carn', 'Carnivory', S.C['carnivorous']), ('para', 'Parasitism', S.C['parasitic']), ('CAM', 'CAM photosynthesis', S.INK),
         ('aquatic', 'Aquatic', S.INK3), ('wet', 'Wet habitat', S.INK3)]
for i, (t, lab, col) in enumerate(terms):
    e, se, p = float(co[t]['estimate']), float(co[t]['SE']), float(co[t]['p']); y = len(terms) - 1 - i
    lo, hi = (e - 1.96 * se) / np.log(10), (e + 1.96 * se) / np.log(10)
    ax.plot([lo, hi], [y, y], color=col, lw=2.6, alpha=0.25, solid_capstyle='butt'); ax.plot([lo, hi], [y, y], color=col, lw=0.8); ax.plot(e / np.log(10), y, 'o', ms=4.4, color=col, mec='white', mew=0.5)
    ps = f'P = {p:.3f}' if p >= 0.001 else r'P = %.1f$\times10^{%d}$' % (p / 10 ** np.floor(np.log10(p)), np.floor(np.log10(p)))
    ax.text(5.3, y, ps, va='center', fontsize=6, color=S.INK2)
ax.axvline(0, color=S.INK2, lw=0.6)
ax.axhspan(2.5, 4.55, color=S.light(S.C['carnivorous'], 0.93), lw=0, zorder=0); ax.axhspan(-0.55, 2.5, color='#f6f5f2', lw=0, zorder=0)
ax.text(7.3, 4.45, 'nutritional strategy', ha='right', va='top', fontsize=5.0, color=S.INK3, style='italic')
ax.text(7.3, 2.4, 'physiology and habitat', ha='right', va='top', fontsize=5.0, color=S.INK3, style='italic')
ax.set_yticks(range(len(terms))); ax.set_yticklabels([t[1] for t in terms][::-1], fontsize=6.3)
ax.set_xticks([-2, 0, 2, 4]); ax.set_xticklabels(['$10^{-2}$', '1', '$10^{2}$', '$10^{4}$']); ax.set_xlim(-2.5, 7.4)
ax.set_xlabel('Odds ratio of ndh loss (95% CI)'); ax.set_ylim(-0.6, 5.3)
ax.xaxis.grid(True, color='white', lw=0.5); ax.set_axisbelow(True)
mods = {r['model']: float(r['dAIC']) for r in csv.DictReader(open(D + 'conf311_GBM_phyloglm_models.csv'))}
ax.set_title(r'logit $P(\mathrm{loss}_i)=\mathbf{x}_i^{\top}\boldsymbol{\beta}$, phylogenetic logistic regression', loc='left', fontsize=6.3, color=S.INK2)
ax.text(-2.4, 4.85, f"Dropping carnivory: ΔAIC = +{mods['no_carnivory']:.1f};  habitat + CAM only: +{mods['habitat_only']:.1f}", fontsize=5.8, color=S.INK)

# ---- D: dating sensitivity ---------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 1]); S.panel(ax, 'd', x=-0.42)
na = list(csv.DictReader(open('/home/claude/work/pan/results/dating_final/results/node_ages.csv')))
A = {(r['run'], r['node']): r for r in na}
nodes = [('crown carnivorous Caryophyllales', 'Crown carniv. Caryophyllales'), ('crown Droseraceae (X1)', 'Crown Droseraceae'),
         ('stem Drosophyllum', 'Stem Drosophyllum'), ('crown Nepenthes', 'Crown Nepenthes'), ('crown Sarraceniaceae', 'Crown Sarraceniaceae'),
         ('stem Roridula (C4)', 'Stem Roridula'), ('stem Cephalotus', 'Stem Cephalotus'), ('stem Byblis', 'Stem Byblis'),
         ('crown Lentibulariaceae', 'Crown Lentibulariaceae'), ('stem Triantha', 'Stem Triantha')]
runs = [('sens_GBM_r1', 'GBM, autocorrelated (preferred)', '#0b0b0b', 'o'), ('post_ML_r1', 'ILN, independent rates', '#3e5373', 's'),
        ('sens_softroot_r1', 'ILN, soft root', '#6d84a6', 'v'), ('sens_pos12_r1', 'ILN, codon pos. 1+2', '#98a9c3', 'D'),
        ('post_SR_r1', 'ILN, SR topology', '#b9c6da', '^')]
for i, (n, nl) in enumerate(nodes):
    y = len(nodes) - 1 - i
    ax.axhspan(y - 0.45, y + 0.45, color='#f6f5f2' if i % 2 == 0 else 'white', lw=0, zorder=0)
    for k, (run, rl, col, m) in enumerate(runs):
        r = A.get((run, n));
        if not r: continue
        yy = y + 0.3 - k * 0.15
        ax.plot([float(r['hpd95_low']), float(r['hpd95_high'])], [yy, yy], color=col, lw=1.1 if k == 0 else 0.7, zorder=2 if k else 3)
        ax.plot(float(r['mean_Ma']), yy, m, ms=3.2 if k == 0 else 2.8, color=col, mec='white', mew=0.3, zorder=4)
ax.set_yticks(range(len(nodes))); ax.set_yticklabels([n[1] for n in nodes][::-1], fontsize=6)
for t in ax.get_yticklabels():
    if any(g in t.get_text() for g in ('Drosophyllum', 'Nepenthes', 'Roridula', 'Cephalotus', 'Byblis', 'Triantha')): t.set_style('italic')
ax.set_xlim(150, 0); ax.set_xlabel('Age (Ma; mean, 95% HPD)'); ax.set_ylim(-0.6, len(nodes) - 0.4)
for b in (66, 33.9): ax.axvline(b, color=S.INK3, lw=0.5, ls=(0, (2, 1.5)), zorder=1)
ax.text(66, len(nodes) - 0.45, 'K–Pg', fontsize=5, color=S.INK3, ha='center', va='bottom'); ax.text(33.9, len(nodes) - 0.45, 'EOT', fontsize=5, color=S.INK3, ha='center', va='bottom')
for run, rl, col, m in runs: ax.plot([], [], m + '-', color=col, ms=3, lw=0.7, label=rl)
ax.legend(loc='upper center', bbox_to_anchor=(0.45, -0.13), ncol=3, fontsize=5.6, handlelength=1.4, columnspacing=0.8)
ax.text(0.99, 1.045, r'log BF$_{\mathrm{GBM\,vs\,ILN}}=48.8\pm0.8$', transform=ax.transAxes, ha='right', va='bottom', fontsize=6.2)
ax.set_axisbelow(True)
S.save(fig, 'Fig4_robustness')
print('ok')
