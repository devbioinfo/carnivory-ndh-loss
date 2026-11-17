"""ED4 | Nuclear survey, gene by gene."""
import csv, collections, numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import figstyle as S
S.setup(); P = '/home/claude/work/pan/nuclear/'
M = list(csv.DictReader(open(P + 'nuclear_matrix.tsv'), delimiter='\t'))
genes = [l.split('\t')[:2] for l in open(P + 'genes.tsv').read().strip().split('\n')]
sp_info = {r['species']: (r['nutrition'], r['family'], r['evidence'], r['controls_present_of_4']) for r in M}
call = {(r['species'], r['gene']): r['call'] for r in M}
ordr = ['Drosera_rotundifolia', 'Drosera_capensis', 'Dionaea_muscipula', 'Cephalotus_follicularis', 'Sarracenia_purpurea', 'Genlisea_aurea', 'Pinguicula_primuliflora',
        'Drosophyllum_lusitanicum', 'Nepenthes_gracilis', 'Triphyophyllum_peltatum', 'Utricularia_aurea', 'Utricularia_gibba', 'Striga_asiatica', 'Cuscuta_australis']
ordr += [s for s in sp_info if s not in ordr]
CC = {'present': S.GS['F'], 'partial': S.GS['P'], 'absent': S.GS['-']}
fig, ax = plt.subplots(figsize=(S.W2, 150 * S.MM)); fig.subplots_adjust(left=0.16, right=0.99, top=0.83, bottom=0.08)
mods = []; 
for i, (m, g) in enumerate(genes):
    y = len(genes) - 1 - i
    for j, sp in enumerate(ordr):
        c = call.get((sp, g), 'NA')
        ax.add_patch(Rectangle((j, y), 0.92, 0.88, fc=CC.get(c, 'white'), ec=S.INK3 if c in ('absent', 'NA') else 'none', lw=0.2, hatch='////' if c == 'NA' else None))
    ax.text(-0.3, y + 0.44, g.replace('_', ' / '), ha='right', va='center', fontsize=4.8, style='italic')
    if not mods or mods[-1][0] != m: mods.append([m, y, y])
    else: mods[-1][1] = y
mlab = {'NDH_subA': 'Subcomplex A', 'NDH_subED': 'Subcomplex E/D', 'NDH_subB': 'Subcomplex B', 'NDH_lumen': 'Lumen', 'NDH_PSI_link': 'PSI linkers', 'NDH_assembly': 'Assembly',
        'NDH_expression': 'Plastid RNA', 'PGR5_route': 'PGR5 route', 'control': 'Controls'}
for m, y0, y1 in mods:
    ax.plot([-5.3, -5.3], [y0 + 0.05, y1 + 0.85], color=S.INK2, lw=0.8, clip_on=False); ax.text(-5.6, (y0 + y1 + 0.9) / 2, mlab[m], ha='right', va='center', fontsize=5, fontweight='bold', clip_on=False)
    if m != 'NDH_subA': ax.axhline(y1 + 0.94, color='white', lw=1.5)
for j, sp in enumerate(ordr):
    nut, fam, evd, ctl = sp_info[sp]; g = 'carnivorous' if nut == 'carnivorous' else 'parasitic' if 'parasit' in nut else 'autotrophic'
    ax.add_patch(Rectangle((j, len(genes) + 0.2), 0.92, 0.6, fc=S.C[g] if g != 'autotrophic' else S.light(S.C['autotrophic'], 0.55), ec='none'))
    ax.text(j + 0.46, len(genes) + 1.1, sp.replace('_', ' '), rotation=60, ha='left', va='bottom', fontsize=4.8, style='italic')
    ax.text(j + 0.46, -0.6, f'{ctl}/4', ha='center', va='top', fontsize=4.3, color=S.INK2)
ax.text(-0.3, -0.6, 'controls present', ha='right', va='top', fontsize=4.5, color=S.INK2)
ax.set_xlim(-0.2, len(ordr)); ax.set_ylim(-1.5, len(genes) + 1); ax.axis('off')
for k, (c, l) in enumerate([(CC['present'], 'present (reciprocal best hit)'), (CC['partial'], 'partial / fragment'), (CC['absent'], 'absent (confirmed in genome)'), ('white', 'not assessed')]):
    ax.add_patch(Rectangle((k * 5.5, -3.2), 0.7, 0.7, fc=c, ec=S.INK3, lw=0.2, hatch='////' if l == 'not assessed' else None, clip_on=False)); ax.text(k * 5.5 + 0.9, -2.85, l, fontsize=5, va='center')
S.save(fig, 'ed/ED5_nuclear_matrix'); print('ok', len(ordr), set(r['call'] for r in M))
