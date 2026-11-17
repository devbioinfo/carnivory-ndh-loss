#!/usr/bin/env python3
"""Pan-plastome gene status for 79 protein-coding genes and 4 rRNAs in every plastome.

Status per gene:  F functional | P pseudogene/fragment | - absent
Evidence combines (1) the GenBank annotation (CDS present, not /pseudo, intact ORF of
>= 80% of the reference length; >= 70% for the length-variable accD, ycf1, ycf2, rpoC2;
ACG/GTG starts allowed) and (2) annotation-independent tblastn/blastn searches with reference
proteins from five well-annotated autotrophs (one per order).

  annotated + ORF intact                              -> F
  annotated, ORF fails, but homology intact (>=90%,
     no stop codons) for an intron-containing gene    -> F  (extraction issue noted)
  annotated + /pseudo or ORF fails                    -> P
  not annotated, homology >= 80% and no stop codons   -> F  (unannotated, flagged)
  not annotated, homology 20-80% or stop codons       -> P
  otherwise                                           -> -
Homology hits must have >= 50% identity for ndh genes (paralogue exclusion) and >= 35% otherwise.
"""
import csv, glob, os, re, subprocess, sys, tempfile
from collections import defaultdict
sys.path.insert(0, '/home/claude/work/repo/scripts')
from genbank_utils import read_genbank, extract, revcomp, _segments

GENES = ('accD atpA atpB atpE atpF atpH atpI ccsA cemA clpP infA matK '
         'ndhA ndhB ndhC ndhD ndhE ndhF ndhG ndhH ndhI ndhJ ndhK '
         'petA petB petD petG petL petN psaA psaB psaC psaI psaJ '
         'psbA psbB psbC psbD psbE psbF psbH psbI psbJ psbK psbL psbM psbN psbT psbZ rbcL '
         'rpl2 rpl14 rpl16 rpl20 rpl22 rpl23 rpl32 rpl33 rpl36 rpoA rpoB rpoC1 rpoC2 '
         'rps2 rps3 rps4 rps7 rps8 rps11 rps12 rps14 rps15 rps16 rps18 rps19 ycf1 ycf2 ycf3 ycf4').split()
RRNA = ['rrn16', 'rrn23', 'rrn4.5', 'rrn5']
INTRON = {'atpF', 'clpP', 'ndhA', 'ndhB', 'petB', 'petD', 'rpl2', 'rpl16', 'rpoC1', 'rps12', 'rps16', 'ycf3'}
VARLEN = {'accD', 'ycf1', 'ycf2', 'rpoC2'}
SYN = {'pafi': 'ycf3', 'pafii': 'ycf4', 'pbf1': 'psbN', 'clpp1': 'clpP', 'ycf9': 'psbZ', 'lhba': 'psbZ',
       'ycf5': 'ccsA', 'ycf10': 'cemA', 'rrn23s': 'rrn23', 'rrn16s': 'rrn16', 'rrn5s': 'rrn5', 'rrn4.5s': 'rrn4.5',
       'rrn4_5': 'rrn4.5', 'rrn4.5s': 'rrn4.5'}
LOWER = {g.lower(): g for g in GENES + RRNA}
REFS = {'Amborella trichopoda': 'NC_005086.1', 'Nepenthes khasiana': None, 'Paulownia': None, 'Actinidia': None, 'Ananas comosus': None}
STOPS = {'TAA', 'TAG', 'TGA'}
B = 'TCAG'; AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
tr = lambda s: ''.join(CODE.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))


def norm(name):
    k = name.strip().replace(' ', '').lower()
    if k.startswith('rrn'):
        k = re.sub(r'(rrn[0-9.]+)s$', r'\1', k)
    k = SYN.get(k, k).lower()
    return LOWER.get(k)


def gname(q):
    for key in ('gene', 'gene_synonym', 'locus_tag'):
        for v in q.get(key, []):
            g = norm(v)
            if g: return g
    for v in q.get('product', []):                 # e.g. 'NdhF', 'Ycf2', 'NADH dehydrogenase subunit F'
        for tok in re.findall(r'[A-Za-z]{3}[A-Z0-9][0-9]?', v):
            g = norm(tok)
            if g: return g
        m = re.search(r'NADH[- ]\w*\s*dehydrogenase subunit ([A-K])\b', v, re.I)
        if m: return 'ndh' + m.group(1).upper()
    return None


def extend_start(seq, loc, nt):
    a, b, strand = _segments(loc)[0]
    best = nt
    for k in range(1, 400):
        if strand == 1:
            s0 = a - 3 * k
            if s0 < 0: break
            cod = seq[s0:s0 + 3].upper()
        else:
            e0 = b + 3 * k
            if e0 > len(seq): break
            cod = revcomp(seq[e0 - 3:e0]).upper()
        if cod in STOPS: break
        if cod in ('ATG', 'ACG', 'GTG'):
            best = (seq[s0:a] + nt) if strand == 1 else (revcomp(seq[b:e0]) + nt)
    return best


def orf_ok(nt, ref_len, frac):
    nt = nt.upper()
    if len(nt) % 3 or len(nt) < 90: return False, 'frame'
    cod = [nt[i:i + 3] for i in range(0, len(nt), 3)]
    if cod[-1] not in STOPS: return False, 'no terminal stop'
    if any(c in STOPS for c in cod[:-1]): return False, 'internal stop'
    if len(nt) < frac * ref_len: return False, f'truncated ({len(nt) / ref_len:.2f})'
    return True, f'intact ({len(nt) / ref_len:.2f})'


def gene_orf(rec, g, reflen):
    for f in rec['features']:
        if f['type'] == 'gene' and gname(f['qualifiers']) == g:
            try: nt = extract(rec['seq'], f['location'])
            except Exception: continue
            frac = 0.7 if g in VARLEN else 0.8
            ok, r = orf_ok(nt, reflen.get(g, len(nt)), frac)
            if not ok:
                nt2 = extend_start(rec['seq'], f['location'], nt)
                ok, r = orf_ok(nt2, reflen.get(g, len(nt2)), frac)
            if ok: return True
    return False


def annotation_calls(rec, reflen):
    st = {g: '-' for g in GENES + RRNA}; why = {g: 'not annotated' for g in GENES + RRNA}
    for f in rec['features']:
        if f['type'] not in ('gene', 'CDS', 'rRNA'): continue
        g = gname(f['qualifiers'])
        if not g: continue
        q = f['qualifiers']
        if g in RRNA:
            if f['type'] == 'rRNA' or f['type'] == 'gene':
                if 'pseudo' in q or 'pseudogene' in q:
                    if st[g] != 'F': st[g], why[g] = 'P', 'flagged /pseudo'
                else:
                    st[g], why[g] = 'F', 'annotated'
            continue
        if 'pseudo' in q or 'pseudogene' in q:
            if st[g] != 'F': st[g], why[g] = 'P', 'flagged /pseudo'
        elif st[g] == '-':
            st[g], why[g] = 'P', 'gene feature only'
    for f in rec['features']:
        if f['type'] != 'CDS': continue
        g = gname(f['qualifiers'])
        if g not in GENES or st[g] == 'F': continue
        q = f['qualifiers']
        if 'pseudo' in q or 'pseudogene' in q: continue
        try:
            nt = extract(rec['seq'], f['location'])
        except Exception:
            why[g] = 'location parse error'; continue
        frac = 0.7 if g in VARLEN else 0.8
        ok, r = orf_ok(nt, reflen.get(g, len(nt)), frac)
        if not ok and ('truncated' in r or 'no terminal' in r or 'internal' in r):
            nt2 = extend_start(rec['seq'], f['location'], nt)
            if len(nt2) > len(nt):
                ok2, r2 = orf_ok(nt2, reflen.get(g, len(nt2)), frac)
                if ok2: ok, r = True, r2 + ', start extended'
        if ok: st[g] = 'F'
        else: st[g] = 'P'
        why[g] = r
    return st, why


def main():
    meta = {r['accession']: r for r in csv.DictReader(open(os.environ.get('TAXA', 'taxa_selected.csv')))}
    gbs = [p for p in sorted(glob.glob('genbank/*.gb')) if os.path.basename(p)[:-3] in meta]
    recs = {}
    for p in gbs:
        rec = next(read_genbank(p)); recs[os.path.basename(p)[:-3]] = rec
    # reference lengths: median intact-looking annotated CDS length among autotrophs
    lens = defaultdict(list)
    for acc, rec in recs.items():
        if meta[acc]['nutrition'] != 'autotrophic': continue
        for f in rec['features']:
            if f['type'] == 'CDS':
                g = gname(f['qualifiers'])
                if g in GENES and 'pseudo' not in f['qualifiers']:
                    try: nt = extract(rec['seq'], f['location'])
                    except Exception: continue
                    if len(nt) % 3 == 0: lens[g].append(len(nt))
    import statistics
    reflen = {g: statistics.median(v) for g, v in lens.items()}
    # reference proteins / rRNAs from one autotroph per order
    refacc = {}
    for order, sp in (('Amborellales', 'Amborella trichopoda'), ('Caryophyllales', 'Nepenthes khasiana'),
                      ('Lamiales', 'Paulownia'), ('Ericales', 'Actinidia'), ('Poales', 'Ananas comosus'), ('Oxalidales', 'Oxalis')):
        for acc, m in meta.items():
            if m['species'].startswith(sp): refacc[order] = acc; break
    prot, rna = {}, {}
    for order, acc in refacc.items():
        rec = recs[acc]
        for f in rec['features']:
            g = gname(f['qualifiers'])
            if not g or 'pseudo' in f['qualifiers']: continue
            try: nt = extract(rec['seq'], f['location']).upper()
            except Exception: continue
            if f['type'] == 'CDS' and g in GENES:
                p = tr(nt).rstrip('*')
                if '*' not in p and len(p) * 3 >= 0.7 * reflen.get(g, 0):
                    prot.setdefault(f'{g}|{order}', p)
            if f['type'] == 'rRNA' and g in RRNA:
                rna.setdefault(f'{g}|{order}', nt)
    print('reference proteins', len(prot), 'rRNAs', len(rna), 'from', refacc, file=sys.stderr)
    tmp = tempfile.mkdtemp()
    qp = os.path.join(tmp, 'q.faa'); open(qp, 'w').write(''.join(f'>{k}\n{v}\n' for k, v in prot.items()))
    qn = os.path.join(tmp, 'q.fna'); open(qn, 'w').write(''.join(f'>{k}\n{v}\n' for k, v in rna.items()))
    plen = {k: len(v) for k, v in prot.items()}; nlen = {k: len(v) for k, v in rna.items()}
    rows, evid = [], []
    for acc, rec in recs.items():
        st, why = annotation_calls(rec, reflen)
        gfa = os.path.join(tmp, 'g.fna'); open(gfa, 'w').write('>g\n' + rec['seq'] + '\n')
        out = subprocess.run(['tblastn', '-query', qp, '-subject', gfa, '-evalue', '1e-5', '-db_gencode', '11', '-seg', 'no',
                              '-outfmt', '6 qseqid qstart qend pident sseq sframe sstart send bitscore'], capture_output=True, text=True).stdout
        hsps = defaultdict(list)
        for line in out.splitlines():
            qid, qs, qe, pid, sseq, sf, ss, se, bs = line.split('\t')
            g = qid.split('|')[0]
            if float(pid) < (50 if g.startswith('ndh') else 35): continue
            hsps[qid].append((float(bs), int(qs), int(qe), sseq, sf, min(int(ss), int(se)), max(int(ss), int(se))))
        hom = {}
        for qid in prot:
            g = qid.split('|')[0]
            hs = sorted(hsps[qid], reverse=True)
            if not hs:
                c, bad = 0.0, False
            else:
                # keep only HSPs at the best-scoring locus (same strand, within the expected gene span)
                b0 = hs[0]; sgn = lambda f: f.startswith('-'); strand = sgn(b0[4])
                span = plen[qid] * 3 + (3000 if g in INTRON else 300)
                if g == 'rps12': span = 10 ** 9            # trans-spliced
                loc = [h for h in hs if sgn(h[4]) == strand and h[5] < b0[6] + span and h[6] > b0[5] - span]
                cov = set()
                for h in loc: cov.update(range(h[1], h[2] + 1))
                c = len(cov) / plen[qid]
                stops = any('*' in h[3].rstrip('*') for h in loc)
                fs = g not in INTRON and len({h[4] for h in loc}) > 1
                bad = stops or fs
            if g not in hom or c > hom[g][0] or (abs(c - hom[g][0]) < 1e-9 and not bad): hom[g] = (c, bad)
        out = subprocess.run(['blastn', '-query', qn, '-subject', gfa, '-evalue', '1e-10', '-outfmt', '6 qseqid qstart qend pident'],
                             capture_output=True, text=True).stdout
        ncov = defaultdict(set)
        for line in out.splitlines():
            qid, qs, qe, pid = line.split('\t')
            if float(pid) >= 70: ncov[qid].update(range(int(qs), int(qe) + 1))
        for qid in rna:
            g = qid.split('|')[0]; c = len(ncov[qid]) / nlen[qid]
            if g not in hom or c > hom[g][0]: hom[g] = (c, False)
        final = {}
        for g in GENES + RRNA:
            c, s = hom.get(g, (0.0, False))
            a = st[g]; note = why[g]
            if a == 'F':
                f_ = 'F'
            elif a == 'P':
                if g in INTRON and c >= 0.9 and not s and note not in ('flagged /pseudo',):
                    f_, note = 'F', note + '; homology intact (intron gene, extraction issue)'
                elif note == 'gene feature only' and c >= 0.8 and not s and gene_orf(rec, g, reflen):
                    f_, note = 'F', f'gene feature only, ORF and homology intact {c:.2f}'
                else:
                    f_ = 'P'
            else:
                if c >= 0.8 and not s: f_, note = 'F', f'unannotated, homology {c:.2f}'
                elif c >= 0.2: f_, note = 'P', f'unannotated fragment, homology {c:.2f}' + (', stop codons' if s else '')
                else: f_ = '-'
            final[g] = f_
            evid.append([acc, meta[acc]['species'], g, a, f'{c:.2f}', int(s), f_, note])
        m = meta[acc]
        nF = sum(final[g] == 'F' for g in GENES)
        rows.append([m['order'], m['family'], m['species'], acc, m['nutrition'], len(rec['seq'])] + [final[g] for g in GENES + RRNA]
                    + [nF, sum(final[g] == 'F' for g in GENES if g.startswith('ndh'))])
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    os.makedirs('results', exist_ok=True)
    with open('results/pan_gene_matrix.csv', 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['order', 'family', 'species', 'accession', 'nutrition', 'length'] + GENES + RRNA + ['n_protein_functional', 'n_ndh_functional']); w.writerows(rows)
    with open('results/pan_gene_evidence.csv', 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['accession', 'species', 'gene', 'annotation_call', 'homology_coverage', 'stop_or_frameshift', 'final', 'note']); w.writerows(evid)
    print('wrote', len(rows), 'plastomes', file=sys.stderr)


if __name__ == '__main__':
    main()
