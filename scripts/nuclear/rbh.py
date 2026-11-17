#!/usr/bin/env python3
"""Nuclear NDH / PGR5 survey by reciprocal best hits (RBH) against the Arabidopsis proteome.
For each species proteome: blastp Arabidopsis queries -> species (e<=1e-10); for the top 5 hits per query,
blastp back to the full Arabidopsis proteome; a hit is an ortholog call if its best Arabidopsis hit is the query
protein (any isoform of the same locus). Call: 'present' if the RBH covers >=50% of the query,
'partial' if <50%, 'absent' if no RBH (to be checked against the genome assembly before being called a loss).
usage: rbh.py   (reads ath.fa, ath_proteome.faa, panel.tsv, prot/*.faa; writes rbh_calls.tsv)"""
import subprocess, os, csv, collections
TH = 2
FAMILY_LEVEL = {'PsaD1', 'PsbO1', 'RbcS1A', 'PGRL1A', 'PGRL1B', 'PetC'}
def run(cmd): subprocess.run(cmd, check=True, capture_output=True)
# Arabidopsis isoforms: map every RefSeq protein to the query locus whose sequence it equals (for the query set),
# and treat proteins with identical sequence or same gene product as one locus
q = {}; n = None
for l in open('ath.fa'):
    if l.startswith('>'): n = l[1:].strip(); q[n] = ''
    else: q[n] += l.strip()
seq = {}; n = None
for l in open('ath_proteome.faa'):
    if l.startswith('>'): n = l[1:].split()[0]; seq[n] = ''
    else: seq[n] += l.strip()
# isoform groups: blast queries against Arabidopsis; proteins with >=98% identity over >=80% of query = same locus
os.makedirs('rbh', exist_ok=True)
run(['blastp', '-query', 'ath.fa', '-db', 'athdb', '-outfmt', '6 qseqid sseqid pident length qlen slen', '-evalue', '1e-20',
     '-max_target_seqs', '50', '-num_threads', str(TH), '-out', 'rbh/ath_self.tsv'])
iso = collections.defaultdict(set)
for l in open('rbh/ath_self.tsv'):
    a, b, pid, ln, ql, sl = l.split('\t');
    if float(pid) >= 98 and int(ln) >= 0.8 * min(int(ql), int(sl)): iso[a].add(b)
    # gene-family level for controls and the Brassicaceae-specific PGRL1A/B duplication: any close paralog counts
    if a.split('|')[0] in FAMILY_LEVEL and float(pid) >= 60 and int(ln) >= 0.7 * int(ql): iso[a].add(b)
rows = []
for sp, acc, nut, fam in csv.reader(open('panel.tsv'), delimiter='\t'):
    faa = f'prot/{sp}.faa'; db = f'rbh/{sp}'
    if not os.path.exists(db + '.pin') and not os.path.exists(db + '.00.pin'):
        run(['makeblastdb', '-in', faa, '-dbtype', 'prot', '-out', db])
    fw = f'rbh/{sp}_fw.tsv'
    if not os.path.exists(fw):
        run(['blastp', '-query', 'ath.fa', '-db', db, '-outfmt', '6 qseqid sseqid pident length qlen slen evalue bitscore qstart qend',
             '-evalue', '1e-10', '-max_target_seqs', '5', '-num_threads', str(TH), '-out', fw])
    hits = collections.defaultdict(list)
    for l in open(fw):
        p = l.rstrip('\n').split('\t'); hits[p[0]].append(p)
    cand = sorted({h[1] for v in hits.values() for h in v})
    # reverse search
    sseq = {}; n = None; want = set(cand)
    for l in open(faa):
        if l.startswith('>'): n = l[1:].split()[0]; keep = n in want
        elif keep: sseq[n] = sseq.get(n, '') + l.strip()
    rv = f'rbh/{sp}_rv.tsv'
    if not os.path.exists(rv):
        open('rbh/tmp.fa', 'w').write(''.join(f'>{k}\n{v}\n' for k, v in sseq.items()))
        run(['blastp', '-query', 'rbh/tmp.fa', '-db', 'athdb', '-outfmt', '6 qseqid sseqid bitscore', '-evalue', '1e-5',
             '-max_target_seqs', '1', '-num_threads', str(TH), '-out', rv])
    best = {}
    for l in open(rv):
        a, b, bs = l.rstrip('\n').split('\t')
        if a not in best or float(bs) > best[a][1]: best[a] = (b, float(bs))
    for g in q:
        call, hit, cov, pid = 'absent', '', 0, 0
        for h in sorted(hits.get(g, []), key=lambda x: -float(x[7])):
            if h[1] in best and best[h[1]][0] in iso[g]:
                # query coverage from all HSPs of this subject
                cov_pos = set()
                for x in hits[g]:
                    if x[1] == h[1]: cov_pos.update(range(int(x[8]), int(x[9]) + 1))
                cov = len(cov_pos) / int(h[4]); pid = float(h[2]); hit = h[1]
                call = 'present' if cov >= 0.5 else 'partial'; break
        rows.append([sp, nut, fam, g, call, hit, round(cov, 2), pid])
with open('rbh_calls.tsv', 'w') as f:
    f.write('species\tnutrition\tfamily\tgene\tcall\thit\tquery_cov\tpident\n')
    for r in rows: f.write('\t'.join(map(str, r)) + '\n')
print('done')
