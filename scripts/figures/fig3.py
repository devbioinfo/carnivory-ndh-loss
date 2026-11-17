"""Figure 3 | Tempo and mode of plastid ndh loss across nutritional strategies (311 taxa, autocorrelated-rates
chronogram; ndh 'lost' = >=3 of 11 genes pseudogenised/absent)."""
import csv, collections, numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import figstyle as S
import icons as I
from matplotlib.patches import Rectangle
S.setup()
D = '/home/claude/work/pan/phaseC/'
rates = list(csv.DictReader(open(D + 'f311_GBM_t3_rates_per_map.csv')))
events = list(csv.DictReader(open(D + 'f311_GBM_t3_loss_events.csv')))
timing = list(csv.DictReader(open(D + 'f311_GBM_t3_timing_per_map.csv')))
nsim = len({r['sim'] for r in rates})
GROUP = {'autotrophic': 'autotrophic', 'carnivorous': 'carnivorous', 'hemiparasitic': 'parasitic', 'holoparasitic': 'parasitic',
         'full_mycoheterotroph': 'mycoheterotrophic', 'partial_mycoheterotroph': 'mycoheterotrophic'}

fig = plt.figure(figsize=(S.W2, 165 * S.MM))
gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.5, left=0.17, right=0.98, top=0.94, bottom=0.08)

# ---- A: loss rate per regime -------------------------------------------------------------
ax = fig.add_subplot(gs[0, 0]); S.panel(ax, 'a')
show = ['autotrophic', 'carnivorous', 'hemiparasitic']
by = collections.defaultdict(list); expo = collections.defaultdict(list)
for r in rates:
    expo[r['regime']].append(float(r['exposure_Myr']))
    if r['rate_per_Myr'] not in ('NA', ''): by[r['regime']].append(float(r['rate_per_Myr']))
med = {}
auto = {r['sim']: float(r['rate_per_Myr']) for r in rates if r['regime'] == 'autotrophic'}
for i, reg in enumerate(show):
    v = np.array(by[reg]); v = v[v > 0]
    lv = np.log10(v)
    parts = ax.violinplot([lv], positions=[i], widths=0.7, showextrema=False)
    col = S.C['parasitic'] if reg == 'hemiparasitic' else S.C[reg]
    for b in parts['bodies']: b.set_facecolor(S.light(col, 0.45)); b.set_edgecolor(col); b.set_linewidth(0.8); b.set_alpha(0.75)
    q = np.percentile(lv, [2.5, 50, 97.5]); qq = np.percentile(lv, [25, 75])
    ax.plot([i, i], [q[0], q[2]], color=S.INK, lw=0.7); ax.plot([i, i], qq, color=S.INK, lw=2.2, solid_capstyle='butt')
    ax.plot(i, q[1], 'o', ms=4, color='white', mec=S.INK, mew=0.8, zorder=5)
    med[reg] = q[1]
    if reg != 'autotrophic':
        ratio = np.array([float(r['rate_per_Myr']) / auto[r['sim']] for r in rates if r['regime'] == reg and r['rate_per_Myr'] not in ('NA', '') and auto[r['sim']] > 0])
        rq = np.percentile(ratio, [2.5, 50, 97.5])
        ax.annotate('', xy=(i - 0.28, q[1]), xytext=(i - 0.28, med['autotrophic']), arrowprops=dict(arrowstyle='-|>', lw=0.7, color=col, mutation_scale=6, shrinkA=0, shrinkB=2))
        ax.plot([0.12, i - 0.28], [med['autotrophic']] * 2, color=S.INK3, lw=0.5, ls=(0, (2, 1.5)))
        ax.text(i - 0.33, (q[1] + med['autotrophic']) / 2, f'{rq[1]:.0f}×\n({rq[0]:.0f}–{rq[2]:.0f})', ha='right', va='center', fontsize=6.2, color=S.INK, fontweight='bold', linespacing=1.1)
ax.set_xticks(range(len(show))); ax.set_xticklabels([S.LABEL[r] for r in show])
for i, (k_, c_) in enumerate([('autotroph', S.C['autotrophic']), ('sundew', S.C['carnivorous']), ('hemiparasite', S.C['parasitic'])]): I.icon_in(ax, k_, [i / 3 + 0.2, 0.02, 0.09, 0.13], c_)
ax.text(0.02, 0.805, 'dot: median; thick bar: IQR; thin bar: 95% interval\nbold numbers: fold change vs autotrophs (95% interval)', transform=ax.transAxes, ha='left', va='top', fontsize=5.6, color=S.INK3, linespacing=1.25, bbox=dict(fc='white', ec='none', pad=0.6))
ax.set_ylabel('ndh loss rate (events per lineage-Myr)')
ticks = [-4, -3, -2, -1]; ax.set_yticks(ticks); ax.set_yticklabels(['$10^{%d}$' % t for t in ticks]); ax.set_ylim(-4.2, -0.15)
ax.yaxis.grid(True, color=S.GRID, lw=0.4); ax.set_axisbelow(True)
ax.text(0.02, 0.97, r'$\hat{\lambda}_r=N_r\,/\,T_r$', transform=ax.transAxes, va='top', fontsize=7.5)
ax.text(0.02, 0.895, r'$N_r$: losses in regime $r$;  $T_r$: time intact in $r$', transform=ax.transAxes, va='top', fontsize=5.8, color=S.INK2)
he = np.median(expo['holoparasitic']); me = np.median(expo['full_mycoheterotroph']) + np.median(expo['partial_mycoheterotroph'])
ax.text(0.98, 0.165, f'Not estimated: holoparasites\n(ndh lost before holoparasitism);\nmycoheterotrophs ({me:.1f} Myr exposure)',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=5.3, color=S.INK2, linespacing=1.2)
ax.set_title(f'Posterior loss rates from {nsim} paired stochastic maps', loc='left', color=S.INK2, fontsize=6.5)

# ---- B: cumulative independent losses through time --------------------------------------------
ax = fig.add_subplot(gs[0, 1]); S.panel(ax, 'b')
grid = np.linspace(0, 140, 281)
cnt = {g: np.zeros((nsim, grid.size)) for g in ['autotrophic', 'carnivorous', 'parasitic', 'mycoheterotrophic']}
sims = sorted({int(r['sim']) for r in rates}); si = {s: k for k, s in enumerate(sims)}
for e in events:
    g = GROUP[e['regime']]; a = float(e['age_Ma'])
    cnt[g][si[int(e['sim'])]] += (grid <= a)
for g in ['autotrophic', 'parasitic', 'carnivorous', 'mycoheterotrophic']:
    m = cnt[g]; lo, md, hi = np.percentile(m, [2.5, 50, 97.5], axis=0)
    ax.fill_between(grid, lo, hi, color=S.C[g], alpha=0.22, lw=0, step='mid', zorder=2)
    ax.step(grid, md, where='mid', color='white', lw=2.4, zorder=3); ax.step(grid, md, where='mid', color=S.C[g], lw=1.4, zorder=4)
    ax.text(-1.5, md[0],  {'autotrophic': 'Autotrophic', 'parasitic': 'Parasitic', 'carnivorous': 'Carnivorous',
            'mycoheterotrophic': 'Mycoheterotrophic'}[g], color=S.INK, fontsize=6.3, va='center', ha='left', clip_on=False)
ax.set_xlim(140, 0); ax.set_xlabel('Age (Ma)'); ax.set_ylabel('Cumulative ndh losses (median, 95% interval)')
ax.set_ylim(-0.5, 15.5); S.geobands(ax, 140, ypos=1.005, alpha_f=0.88, short_below=15)
ax.yaxis.grid(True, color='white', lw=0.5); ax.set_axisbelow(True)
for b, lab in [(66, 'K–Pg'), (33.9, 'EOT')]: ax.text(b - 1, 15.2, lab, fontsize=5.3, color=S.INK2, va='top')
ax.set_title('Independent losses accumulate in each regime', loc='left', color=S.INK2, fontsize=6.8, pad=16)

# ---- C: timing of carnivory origin vs first ndh loss --------------------------------------
ax = fig.add_subplot(gs[1, 0]); S.panel(ax, 'c')
lin = collections.defaultdict(lambda: {'c': [], 'l': [], 'nl': 0, 'n': 0})
for r in timing:
    d = lin[r['lineage']]; d['n'] += 1
    if r['carnivory_origin_Ma'] not in ('NA', ''): d['c'].append(float(r['carnivory_origin_Ma']))
    if r['first_ndh_loss_Ma'] not in ('NA', ''): d['l'].append(float(r['first_ndh_loss_Ma']))
    else: d['nl'] += 1
    if r['carnivory_origin_Ma'] not in ('NA', '') and r['first_ndh_loss_Ma'] not in ('NA', ''):
        d.setdefault('after', []).append(float(r['first_ndh_loss_Ma']) < float(r['carnivory_origin_Ma']))
order = ['Lentibulariaceae', 'Caryophyllales_carnivores', 'Sarraceniaceae', 'Cephalotus', 'Roridula', 'Byblis', 'Triantha', 'Brocchinia']
name = {'Caryophyllales_carnivores': 'Carnivorous Caryophyllales', 'Brocchinia': 'Brocchinia (carnivorous)'}
for k, L in enumerate(order):
    d = lin[L]; y = len(order) - 1 - k
    for vals, col, dy, mk in [(d['c'], S.C['carnivorous'], 0.15, 'o'), (d['l'], S.INK, -0.15, 's')]:
        if len(vals) < 0.2 * d['n']: continue
        q = np.percentile(vals, [2.5, 50, 97.5])
        ax.plot([q[0], q[2]], [y + dy] * 2, color=col, lw=3.2, alpha=0.22, solid_capstyle='butt', zorder=2)
        ax.plot([q[0], q[2]], [y + dy] * 2, color=col, lw=0.7, solid_capstyle='butt', zorder=3); ax.plot(q[1], y + dy, mk, ms=3.8, color=col, mec='white', mew=0.5, zorder=4)
        d.setdefault('med', {})[mk] = q[1]
    lab = name.get(L, L); it = L in ('Cephalotus', 'Roridula', 'Byblis', 'Triantha')
    import matplotlib.transforms as mtr
    ax.text(-0.25, y, lab, ha='right', va='center', fontsize=6.3, style='italic' if it else 'normal', transform=mtr.blended_transform_factory(ax.transAxes, ax.transData))
    after = d.get('after', [])
    if after and np.mean(after) >= 0.95:
        ax.annotate('', xy=(d['med']['s'], y - 0.15), xytext=(d['med']['o'], y + 0.15), arrowprops=dict(arrowstyle='-|>', lw=0.8, color=S.INK, mutation_scale=6, shrinkA=2, shrinkB=2), zorder=5)
    txt = 'no loss' if d['nl'] / d['n'] > 0.8 else f"P(after) = {np.mean(after):.2f}" if after else ''
    pa = np.mean(after) if after else None
    txt = 'no loss' if d['nl'] / d['n'] > 0.8 else (f'{pa:.2f} after' if pa >= 0.95 else f'{pa:.2f} unresolved') if pa is not None else ''
    ax.text(-2, y, txt, ha='left', va='center', fontsize=5.6,
            color=S.INK if (pa is not None and pa >= 0.95) else S.INK2, fontweight='bold' if (pa is not None and pa >= 0.95) else 'normal')
    ikind = {'Lentibulariaceae': ['bladderwort'], 'Caryophyllales_carnivores': ['sundew', 'nepenthes'], 'Sarraceniaceae': ['trumpet'], 'Cephalotus': ['cephalotus'],
             'Roridula': ['roridula'], 'Byblis': ['byblis'], 'Triantha': ['triantha'], 'Brocchinia': ['tank']}[L]
    for m_, kk in enumerate(ikind):
        I.icon_in(ax, kk, [-0.075 - m_ * 0.07, (y + 0.7) / (len(order) + 0.4) - 0.05, 0.065, 0.1], S.C['carnivorous'])
ax.set_xlim(140, 0); ax.set_ylim(-0.7, len(order) - 0.3); ax.set_yticks([]); ax.spines['left'].set_visible(False)
S.geobands(ax, 140, labels=False, alpha_f=0.88)
ax.set_xlabel('Age (Ma)')
ax.text(-2, len(order) - 0.45, 'P(loss after carnivory)', fontsize=5.4, color=S.INK2, ha='left')
ax.plot([], [], 'o-', color=S.C['carnivorous'], ms=3.5, lw=1, label='Origin of carnivory'); ax.plot([], [], 's-', color=S.INK, ms=3.5, lw=1, label='First ndh loss')
ax.legend(loc='upper left', bbox_to_anchor=(0.0, 1.13), ncol=2, handlelength=1.6)

# ---- D: Pagel dependent model ------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 1]); S.panel(ax, 'd', x=-0.04)
ax.set_xlim(0, 10); ax.set_ylim(0, 7.6); ax.axis('off')
Q = {}; rows = list(csv.reader(open(D + 'pagelQ_carn_depx.csv'))); hdr = rows[0][1:]
for r in rows[1:]:
    for j, v in enumerate(r[1:]): Q[(r[0], hdr[j])] = float(v)
pos = {'intact|autotrophic': (2.0, 5.4), 'intact|carnivorous': (8.0, 5.4), 'lost|autotrophic': (2.0, 1.75), 'lost|carnivorous': (8.0, 1.75)}
ax.add_patch(FancyBboxPatch((0.55, 0.85), 2.9, 5.95, boxstyle='round,pad=0,rounding_size=0.3', fc=S.light(S.C['autotrophic'], 0.85), ec='none', zorder=0))
ax.add_patch(FancyBboxPatch((6.55, 0.85), 2.9, 5.95, boxstyle='round,pad=0,rounding_size=0.3', fc=S.light(S.C['carnivorous'], 0.85), ec='none', zorder=0))
ax.text(2.0, 6.33, 'Autotrophic', ha='center', fontsize=6.2, fontweight='bold', color=S.INK2); ax.text(8.0, 6.33, 'Carnivorous', ha='center', fontsize=6.2, fontweight='bold', color=S.DARK['carnivorous'])
for s_, (x, y) in pos.items():
    lost = s_.startswith('lost'); col = S.C['carnivorous'] if 'carn' in s_ else S.INK3
    ax.add_patch(FancyBboxPatch((x - 1.15, y - 0.6), 2.3, 1.2, boxstyle='round,pad=0.02,rounding_size=0.18', fc='white', ec=col, lw=1.0, ls='--' if lost else '-', zorder=2))
    ax.text(x, y + 0.27, 'ndh lost' if lost else 'ndh intact', ha='center', va='center', fontsize=6.2, style='normal', zorder=3)
    for j in range(11):   # 11-gene glyph in the gene-state colours
        st = ('P' if j % 3 else '-') if lost else 'F'
        ax.add_patch(Rectangle((x - 0.88 + j * 0.16, y - 0.38), 0.13, 0.3, fc=S.GS[st], ec=S.INK3 if st == '-' else 'none', lw=0.2, zorder=3))
def arrow(p1, p2, rate, label, lab_xy, col):
    w = 0.5 + 1.1 * np.log10(rate / 2e-4) if rate > 0 else 0.5
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle='-|>', mutation_scale=7 + w, lw=max(w, 0.5), color=col, zorder=1.5))
    ax.text(*lab_xy, label, fontsize=6.0, ha='center', va='center', bbox=dict(fc='white', ec='none', pad=0.8), zorder=4)
qa = Q[('intact|autotrophic', 'lost|autotrophic')]; qc = Q[('intact|carnivorous', 'lost|carnivorous')]
qg = Q[('intact|autotrophic', 'intact|carnivorous')]; ql = Q[('lost|autotrophic', 'lost|carnivorous')]
arrow((2.0, 4.75), (2.0, 2.4), qa, f'$q_{{IL}}^{{A}}$ = {qa*1e3:.2f}', (2.0, 3.58), S.INK2)
arrow((8.0, 4.75), (8.0, 2.4), qc, f'$q_{{IL}}^{{C}}$ = {qc*1e3:.1f}', (8.0, 3.58), S.DARK['carnivorous'])
arrow((3.2, 5.55), (6.8, 5.55), qg, f'{qg*1e3:.2f}', (5.0, 5.55), S.INK3)
arrow((3.2, 1.6), (6.8, 1.6), ql, f'{ql*1e3:.2f}', (5.0, 1.6), S.INK3)
ax.text(5.0, 3.58, f'$q_{{IL}}^{{C}}\\,/\\,q_{{IL}}^{{A}}$ = {qc/qa:.0f}', ha='center', va='center', fontsize=7.0, fontweight='bold')
ax.text(5.0, 3.05, 'faster ndh loss\nin carnivorous lineages', ha='center', va='top', fontsize=5.4, color=S.INK2)
ax.text(5.0, 5.9, 'gain of carnivory', ha='center', fontsize=5.2, color=S.INK3)
ax.text(0.2, 7.75, r'$\mathbf{Q}_{\mathrm{dep}}$: ndh loss rate depends on diet' + '\n' + r'$q\times10^{-3}$ Myr$^{-1}$; arrow width $\propto\log q$', fontsize=6.3, va='top')
ax.add_patch(FancyBboxPatch((0.25, 0.0), 9.5, 0.62, boxstyle='round,pad=0,rounding_size=0.15', fc='#f6f5f2', ec='none'))
ax.text(5.0, 0.31, r'LR $=2(\ell_{\mathrm{dep}}-\ell_{\mathrm{ind}})=43.2$  ·  $\chi^2_{2}$, $P=4.1\times10^{-10}$  ·  $\Delta$AIC $=39.2$', fontsize=6.0, ha='center', va='center')
ax.text(5.0, -0.28, 'A statistical association between states; it does not by itself identify a molecular mechanism.', fontsize=4.9, ha='center', va='top', color=S.INK3, style='italic')
S.save(fig, 'Fig3_tempo_mode')
print('ok', qa, qc, qg)
