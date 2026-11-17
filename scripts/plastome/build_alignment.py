#!/usr/bin/env python3
"""Concatenated codon alignment for the backbone tree.

One intact annotated copy of each protein-coding gene per plastome (ORF test as in pangene.py),
genes present in >= 50% of plastomes, proteins aligned with MAFFT (--auto), back-translated,
codon columns with gaps in > 50% of the sequences of that gene removed, concatenated with a
NEXUS partition file (one partition per gene).
"""
import csv, glob, os, statistics, subprocess, sys
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pangene as P
from genbank_utils import read_genbank, extract

TAXA = os.environ.get('TAXA', 'taxa_selected.csv'); PHY = os.environ.get('PHY', 'results/phylogeny')
meta = {r['accession']: r for r in csv.DictReader(open(TAXA))}
recs = {os.path.basename(p)[:-3]: next(read_genbank(p)) for p in sorted(glob.glob('genbank/*.gb')) if os.path.basename(p)[:-3] in meta}
assert len(recs) == len(meta), set(meta) - set(recs)
name = {a: m['species'].replace(' x ', '_x_').replace(' ', '_').replace('.', '') for a, m in meta.items()}
lens = defaultdict(list)
for acc, rec in recs.items():
    if meta[acc]['nutrition'] != 'autotrophic': continue
    for f in rec['features']:
        if f['type'] == 'CDS':
            g = P.gname(f['qualifiers'])
            if g in P.GENES and 'pseudo' not in f['qualifiers']:
                try: nt = extract(rec['seq'], f['location'])
                except Exception: continue
                if len(nt) % 3 == 0: lens[g].append(len(nt))
reflen = {g: statistics.median(v) for g, v in lens.items()}

seqs = defaultdict(dict)
for acc, rec in recs.items():
    for f in rec['features']:
        if f['type'] != 'CDS' or 'pseudo' in f['qualifiers']: continue
        g = P.gname(f['qualifiers'])
        if g not in P.GENES or acc in seqs[g]: continue
        try: nt = extract(rec['seq'], f['location']).upper()
        except Exception: continue
        frac = 0.7 if g in P.VARLEN else 0.8
        ok, _ = P.orf_ok(nt, reflen.get(g, len(nt)), frac)
        if not ok:
            nt2 = P.extend_start(rec['seq'], f['location'], nt)
            ok, _ = P.orf_ok(nt2, reflen.get(g, len(nt2)), frac); nt = nt2
        if ok: seqs[g][acc] = nt[:-3]

N = len(recs)
os.makedirs(f'{PHY}/genes', exist_ok=True)
concat = {a: [] for a in recs}; parts = []; pos = 1
for g in P.GENES:
    if len(seqs[g]) < 0.5 * N: continue
    prot = {a: ''.join(P.CODE.get(s[i:i + 3], 'X') for i in range(0, len(s), 3)) for a, s in seqs[g].items()}
    fa = f'{PHY}/genes/{g}.faa'
    open(fa, 'w').write(''.join(f'>{a}\n{p}\n' for a, p in prot.items()))
    aln = subprocess.run(['mafft', '--auto', '--quiet', '--thread', '2', fa], capture_output=True, text=True).stdout
    al = {}; cur = None
    for line in aln.splitlines():
        if line.startswith('>'): cur = line[1:].strip(); al[cur] = []
        else: al[cur].append(line.strip())
    al = {a: ''.join(v).upper() for a, v in al.items()}
    L = len(next(iter(al.values())))
    cod = {}
    for a, pa in al.items():
        s = seqs[g][a]; k = 0; out = []
        for ch in pa:
            if ch == '-': out.append('---')
            else: out.append(s[3 * k:3 * k + 3]); k += 1
        cod[a] = out
    keep = [j for j in range(L) if sum(al[a][j] == '-' for a in al) <= 0.5 * len(al)]
    for a in recs:
        concat[a].append(''.join(cod[a][j] for j in keep) if a in cod else '---' * len(keep))
    parts.append((g, pos, pos + 3 * len(keep) - 1, len(seqs[g]))); pos += 3 * len(keep)
    print(g, len(seqs[g]), 'taxa', len(keep), 'codons', flush=True)

with open(f'{PHY}/concat.fasta', 'w') as fh:
    for a in recs: fh.write(f'>{name[a]}\n{"".join(concat[a])}\n')
with open(f'{PHY}/partitions.nex', 'w') as fh:
    fh.write('#nexus\nbegin sets;\n')
    for g, s, e, n in parts: fh.write(f'  charset {g} = {s}-{e};\n')
    fh.write('end;\n')
with open(f'{PHY}/genes_used.tsv', 'w') as fh:
    fh.write('gene\tstart\tend\tn_taxa\n'); fh.writelines(f'{g}\t{s}\t{e}\t{n}\n' for g, s, e, n in parts)
print('genes', len(parts), 'sites', pos - 1)
