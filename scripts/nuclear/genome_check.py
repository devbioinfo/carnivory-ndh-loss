#!/usr/bin/env python3
"""Genome-level check of the nuclear survey (catches genes missed by annotation).
tblastn Arabidopsis queries against a genome assembly; for each query, the best-scoring locus (HSPs on one
contig within 20 kb, same strand) is translated per HSP and sent back to the Arabidopsis proteome (blastp);
call 'present' if the reciprocal best hit is the query (or its isoform/family as in rbh.py) and HSPs cover
>=50% of the query, 'fragment' if reciprocal but <50% (possible pseudogene or assembly gap), 'absent' otherwise.
usage: genome_check.py <species> <genome.fna>"""
import sys, os, subprocess, collections, re
sp, fna = sys.argv[1:3]; TH = '2'
out = f'gchk/{sp}'; os.makedirs('gchk', exist_ok=True)
FAMILY_LEVEL = {'PsaD1', 'PsbO1', 'RbcS1A', 'PGRL1A', 'PGRL1B', 'PetC'}
if not os.path.exists(out + '.ndb') and not os.path.exists(out + '.00.nsq') and not os.path.exists(out + '.nsq'):
    subprocess.run(['makeblastdb', '-in', fna, '-dbtype', 'nucl', '-out', out, '-max_file_sz', '4GB'], check=True, capture_output=True)
tb = out + '_tbn.tsv'
if not os.path.exists(tb):
    subprocess.run(['tblastn', '-query', 'ath.fa', '-db', out, '-evalue', '1e-8', '-outfmt',
                    '6 qseqid sseqid qstart qend sstart send bitscore qlen sseq', '-max_target_seqs', '10',
                    '-num_threads', TH, '-out', tb], check=True)
iso = collections.defaultdict(set)
for l in open('rbh/ath_self.tsv'):
    a, b, pid, ln, ql, sl = l.split('\t')
    if float(pid) >= 98 and int(ln) >= 0.8 * min(int(ql), int(sl)): iso[a].add(b)
    if a.split('|')[0] in FAMILY_LEVEL and float(pid) >= 60 and int(ln) >= 0.7 * int(ql): iso[a].add(b)
H = collections.defaultdict(list)
for l in open(tb):
    q, s, qs, qe, ss, se, bs, ql, sq = l.rstrip('\n').split('\t')
    H[q].append((s, int(qs), int(qe), int(ss), int(se), float(bs), int(ql), sq))
queries = [l[1:].strip() for l in open('ath.fa') if l.startswith('>')]
res = []
for q in queries:
    hs = sorted(H.get(q, []), key=lambda x: -x[5])
    if not hs: res.append((q, 'absent', 0, '')); continue
    top = hs[0]; strand = top[4] > top[3]
    loc = [h for h in hs if h[0] == top[0] and (h[4] > h[3]) == strand and abs(min(h[3], h[4]) - min(top[3], top[4])) < 20000]
    cov = set()
    for h in loc: cov.update(range(h[1], h[2] + 1))
    covf = len(cov) / top[6]
    pep = ''.join(h[7].replace('-', '') for h in sorted(loc, key=lambda x: x[1]))
    pep = re.sub(r'[^A-Z]', 'X', pep.replace('*', 'X'))
    open('gchk/tmp.fa', 'w').write(f'>x\n{pep}\n')
    o = subprocess.run(['blastp', '-query', 'gchk/tmp.fa', '-db', 'athdb', '-outfmt', '6 sseqid bitscore', '-max_target_seqs', '1', '-evalue', '1e-5'],
                       capture_output=True, text=True).stdout.split('\n')
    back = o[0].split('\t')[0] if o and o[0] else ''
    stops = sum(h[7].count('*') for h in loc)
    if back in iso[q]: call = 'present' if covf >= 0.5 else 'fragment'
    else: call = 'absent'
    res.append((q, call, round(covf, 2), f'{top[0]}:{min(top[3],top[4])} stops={stops}'))
with open(f'gchk/{sp}_calls.tsv', 'w') as f:
    f.write('gene\tcall\tquery_cov\tlocus\n')
    for r in res: f.write('\t'.join(map(str, r)) + '\n')
print(sp, collections.Counter(r[1] for r in res))
