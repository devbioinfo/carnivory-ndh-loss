"""ED10 | Read-level support for Triantha ndh calls and ORF completeness of Byblis transcript assemblies."""
import csv, numpy as np, matplotlib.pyplot as plt
import figstyle as S
S.setup(); P = '/home/claude/work/pan/results/phaseA/'
tr = list(csv.DictReader(open(P + 'triantha_ndh_readsupport.tsv'), delimiter='\t'))
calls = dict(zip('ABCDEFGHIJK', next(r['ndh_calls'] for r in csv.DictReader(open('/home/claude/work/pan/phaseC/phaseC_traits_cov311.csv')) if r['tip'] == 'Triantha_occidentalis')))
fig = plt.figure(figsize=(S.W2, 70 * S.MM)); gs = fig.add_gridspec(1, 2, wspace=0.35, left=0.08, right=0.98, top=0.8, bottom=0.2)
ax = fig.add_subplot(gs[0]); S.panel(ax, 'a', x=-0.12)
names = [r['gene'].replace('_locus_ndhE-ndhI', ' locus') for r in tr]
for i, r in enumerate(tr):
    g = r['gene'][3] if r['gene'].startswith('ndh') else '?'
    st = 'deleted' if 'locus' in r['gene'] else calls.get(g, 'F')
    col = S.GS['F'] if st == 'F' else S.GS['P'] if st == 'P' else S.light(S.GS['P'],0.5)
    ax.bar(i, float(r['mean_depth']), color=col, width=0.7, ec=S.INK2 if st != 'F' else 'none', lw=0.4)
    ax.plot([i - 0.3, i + 0.3], [float(r['min_depth'])] * 2, color=S.INK, lw=0.8)
    ax.text(i, float(r['mean_depth']) + 10, f"{r['positions_<80%_agree']}", ha='center', fontsize=5, color=S.INK2)
ax.axhline(160, color=S.INK3, lw=0.5, ls=(0, (2, 1.5))); ax.text(10.4, 165, '160×', fontsize=5, color=S.INK3, ha='right')
ax.set_xticks(range(len(tr))); ax.set_xticklabels(names, rotation=45, ha='right', fontsize=5.4, style='italic'); ax.set_ylabel('Read depth (bar mean; tick minimum)')
ax.set_title('Triantha occidentalis: depth at each ndh locus; number above bar =\nsites with <80% read agreement. Dark: functional; grey outlined:\npseudogene; ndhG locus: deletion spanned by reads', fontsize=5.4, loc='left', color=S.INK2)
ax = fig.add_subplot(gs[1]); S.panel(ax, 'b', x=-0.12)
by = list(csv.DictReader(open(P + 'byblis_ndh_orfcalls.tsv'), delimiter='\t'))
for i, r in enumerate(by):
    f = float(r['len_frac_of_ref']); ok = r['call'] == 'intact'
    ax.bar(i, f, color=S.GS['F'] if ok else 'white', ec=S.INK2, hatch=None if ok else '////', lw=0.4, width=0.7)
    ax.text(i, f + 0.02, f'{f:.2f}', ha='center', fontsize=5)
ax.axhline(1, color=S.INK3, lw=0.5, ls=(0, (2, 1.5))); ax.set_ylim(0, 1.25)
ax.set_xticks(range(len(by))); ax.set_xticklabels([r['gene'] for r in by], style='italic', fontsize=5.6)
ax.set_ylabel('Assembled ORF length / reference CDS')
ax.set_title('Byblis gigantea transcripts: intact ORFs (dark) vs truncated or\nframe-broken (hatched); ndhA, ndhB, ndhE unassembled. Truncation is\nnot separable from low-coverage assembly fragmentation', fontsize=5.4, loc='left', color=S.INK2)
S.save(fig, 'ed/ED10_read_support'); print('ok')
