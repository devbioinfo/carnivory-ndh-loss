#!/usr/bin/env python3
"""Prepare MCMCTree (approximate likelihood) inputs for the five topologies.

Inputs : concat.fasta (295 taxa, codon alignment), taxa_phylo295.csv,
         IQ-TREE trees ../mac_run_v7/out/{ml295_r1,c_SR,c_H4,c_H5,c_H6}.treefile (override with TREEDIR=...)
Outputs: aln_both.phy (2 loci: codon positions 1+2, position 3), aln_pos12.phy (1 locus)
         trees/<topo>_<calset>.tre  (rooted on Amborella, calibrations in MCMCTree notation, time unit 100 Myr)
         runs/<run>/mcmctree.ctl     for every Hessian, prior and posterior run
         calibration_nodes.tsv       which clade each calibration landed on, per topology (check this)
Calibrations (approved table, ICS v2024/12 stage tops):
  C1 root      main: hard uniform 136.0-139.35      sens: soft 125.9-247.3
  C2 crown mesangiosperms  min 123.0 / 121.4
  C3 crown core Ericales (Sarraceniaceae, Roridulaceae, Actinidiaceae, Clethraceae, Ericaceae)  min 85.7 / 89.8
  C4 stem Roridula (Roridula + its sister)   min 35.0 / 33.9
  C5 stem Aldrovanda (Aldrovanda + its sister) min 33.9 / 33.9
Minimum bounds use MCMCTree's truncated-Cauchy L(tL, p=0.1, c=1, pL=0.025) (soft minimum, heavy right tail).
"""
import os, sys, statistics
from Bio import Phylo

TREEDIR = os.environ.get('TREEDIR', '../mac_run_v7/out')
# ML = Ericaceae + (Roridula, Actinidiaceae) (identical to c_H3 with 295 taxa); SR = Sarraceniaceae + (Roridula, Actinidiaceae)
TOPOS = {'ML': 'ml295_r1', 'SR': 'c_SR', 'H4': 'c_H4', 'H5': 'c_H5', 'H6': 'c_H6'}
CAL = {
    'main': {'C1': "B(1.360,1.3935,1e-300,1e-300)", 'C2': 1.230, 'C3': 0.857, 'C4': 0.350, 'C5': 0.339},
    'softroot': {'C1': "B(1.259,2.473,0.025,0.025)", 'C2': 1.230, 'C3': 0.857, 'C4': 0.350, 'C5': 0.339},
    'sensmin': {'C1': "B(1.360,1.3935,1e-300,1e-300)", 'C2': 1.214, 'C3': 0.898, 'C4': 0.339, 'C5': 0.339},
}
def L(t): return f"L({t:.4f},0.1,1,0.025)"

# ---------- alignment ----------
seqs, order, cur = {}, [], None
for line in open('concat.fasta'):
    line = line.strip()
    if line.startswith('>'): cur = line[1:]; seqs[cur] = []; order.append(cur)
    elif cur: seqs[cur].append(line)
seqs = {k: ''.join(v).upper() for k, v in seqs.items()}
# dating subset (taxa_dating.txt, 124 taxa chosen by select_dating_taxa.py; full 295-taxon run was ~3 days per chain)
KEEP = [l.strip() for l in open('taxa_dating.txt') if l.strip()] if os.path.exists('taxa_dating.txt') else order
missing = [k for k in KEEP if k not in seqs]
if missing: print('taxa_dating.txt names not in alignment:', missing); sys.exit(1)
order = [k for k in order if k in set(KEEP)]
print('dating subset:', len(order), 'taxa')
n = len(next(iter(seqs.values())))
keepcols = [i for i in range(n) if any(seqs[k][i] not in '-?N' for k in order)]
def phy(cols_fn):
    out = []
    for k in order:
        s = seqs[k]; out.append(''.join(s[i] for i in keepcols if cols_fn(i)))
    L_ = len(out[0])
    return f"{len(order)} {L_}\n" + ''.join(f"{k}  {s}\n" for k, s in zip(order, out)) + '\n'
p12 = phy(lambda i: i % 3 != 2); p3 = phy(lambda i: i % 3 == 2)
open('aln_both.phy', 'w').write(p12 + p3)
open('aln_pos12.phy', 'w').write(p12)
print('alignments written: pos12', p12.split('\n')[0], '| pos3', p3.split('\n')[0])

# ---------- trees ----------
os.makedirs('trees', exist_ok=True)
log = open('calibration_nodes.tsv', 'w'); log.write('topology\tcalibration\tn_tips\tfamilies_or_sister\n')
fam = {r.split(',')[2].replace(' x ', '_x_').replace(' ', '_').replace('.', ''): r.split(',')[1]
       for r in open('taxa_phylo295.csv').read().splitlines()[1:]}
rates = {}
def tipnames(c): return [t.name for t in c.get_terminals()]
def find(t, prefix): return next(x for x in t.get_terminals() if x.name.startswith(prefix))
for topo, fn in TOPOS.items():
    path = f'{TREEDIR}/{fn}.treefile'
    if not os.path.exists(path): print('MISSING', path); sys.exit(1)
    t = Phylo.read(path, 'newick')
    t.root_with_outgroup('Amborella_trichopoda')
    ks = set(order)
    for tip in [x.name for x in t.get_terminals() if x.name not in ks]: t.prune(tip)
    assert len(t.get_terminals()) == len(order)
    # mean rate prior: median root-to-tip (subs/site) / 1.38 (x100 Myr)
    d = t.depths(); rtt = [d[x] for x in t.get_terminals()]
    rates[topo] = statistics.median(rtt) / 1.38
    for c in t.find_clades():
        c.confidence = None
        if not c.is_terminal(): c.name = None
    root = t.root
    meso = [c for c in root.clades if not (c.is_terminal() and c.name == 'Amborella_trichopoda')][0]
    core = t.common_ancestor([find(t, g).name for g in ('Sarracenia', 'Roridula', 'Actinidia', 'Clethra', 'Rhododendron')])
    def stem(prefix):
        leaf = find(t, prefix); p = t.get_path(leaf)
        return p[-2] if len(p) > 1 else root
    rori, aldro = stem('Roridula'), stem('Aldrovanda')
    nodes = {'C1': root, 'C2': meso, 'C3': core, 'C4': rori, 'C5': aldro}
    assert len({id(v) for v in nodes.values()}) == 5, f'{topo}: two calibrations on one node'
    for k, c in nodes.items():
        tips = tipnames(c)
        if k in ('C4', 'C5'):
            info = 'sister: ' + ','.join(sorted({fam.get(x, '?') for x in tips if not x.startswith(('Roridula', 'Aldrovanda'))}))
        else:
            info = ','.join(sorted({fam.get(x, '?') for x in tips})) if k == 'C3' else ''
        log.write(f'{topo}\t{k}\t{len(tips)}\t{info}\n')
        if k == 'C3':
            outside = {fam.get(x) for x in tips} - {'Sarraceniaceae', 'Roridulaceae', 'Actinidiaceae', 'Clethraceae', 'Ericaceae'}
            if outside: print(f'WARNING {topo}: crown core Ericales also contains {outside}')
    for cs, cal in CAL.items():
        lab = {id(nodes['C1']): cal['C1']}
        for k in ('C2', 'C3', 'C4', 'C5'): lab[id(nodes[k])] = L(cal[k])
        def nw(c):
            if c.is_terminal(): return c.name
            s = '(' + ','.join(nw(x) for x in c.clades) + ')'
            if id(c) in lab: s += "'" + lab[id(c)] + "'"
            return s
        assert all(len(c.clades) == 2 for c in t.get_nonterminals()), f'{topo}: tree not binary'
        open(f'trees/{topo}_{cs}.tre', 'w').write(f"{len(t.get_terminals())} 1\n{nw(root)};\n")
log.close()
print('mean rate (subs/site/100 Myr) per topology:', {k: round(v, 4) for k, v in rates.items()})

# ---------- control files ----------
rate = statistics.median(rates.values())
RG = f"2 {2 / rate:.2f} 1"   # gamma(2, 2/mean): mean = observed rate
def ctl(run, aln, tree, usedata, clock, seed, ndata):
    os.makedirs(f'runs/{run}', exist_ok=True)
    open(f'runs/{run}/mcmctree.ctl', 'w').write(f"""          seed = {seed}
       seqfile = ../../{aln}
      treefile = ../../trees/{tree}.tre
      mcmcfile = mcmc.txt
       outfile = out.txt
         ndata = {ndata}
       seqtype = 0
       usedata = {usedata}
         clock = {clock}
         model = 4
         alpha = 0.5
         ncatG = 5
     cleandata = 0
       BDparas = 1 1 0.1 C
   kappa_gamma = 6 2
   alpha_gamma = 1 1
   rgene_gamma = {RG}
  sigma2_gamma = 1 10 1
         print = 1
        burnin = 20000
      sampfreq = {200 if usedata != 0 else 20}
       nsample = 10000
""")
jobs_h, jobs_m = [], []
for topo in TOPOS:
    ctl(f'hessian_{topo}', 'aln_both.phy', f'{topo}_main', 3, 2, 1, 2); jobs_h.append(f'hessian_{topo}')
    ctl(f'prior_{topo}', 'aln_both.phy', f'{topo}_main', 0, 2, 11, 2); jobs_m.append(f'prior_{topo}|-')
    for r, sd in (('r1', 101), ('r2', 202)):
        ctl(f'post_{topo}_{r}', 'aln_both.phy', f'{topo}_main', '2 in.BV', 2, sd, 2); jobs_m.append(f'post_{topo}_{r}|hessian_{topo}')
for r, sd in (('r1', 101), ('r2', 202)):
    ctl(f'sens_GBM_{r}', 'aln_both.phy', 'ML_main', '2 in.BV', 3, sd, 2); jobs_m.append(f'sens_GBM_{r}|hessian_ML')
    ctl(f'sens_softroot_{r}', 'aln_both.phy', 'ML_softroot', '2 in.BV', 2, sd, 2); jobs_m.append(f'sens_softroot_{r}|hessian_ML')
    ctl(f'sens_minima_{r}', 'aln_both.phy', 'ML_sensmin', '2 in.BV', 2, sd, 2); jobs_m.append(f'sens_minima_{r}|hessian_ML')
    ctl(f'sens_pos12_{r}', 'aln_pos12.phy', 'ML_main', '2 in.BV', 2, sd, 1); jobs_m.append(f'sens_pos12_{r}|hessian_ML:pos12')
ctl('prior_ML_softroot', 'aln_both.phy', 'ML_softroot', 0, 2, 11, 2); jobs_m.append('prior_ML_softroot|-')
ctl('prior_ML_sensmin', 'aln_both.phy', 'ML_sensmin', 0, 2, 11, 2); jobs_m.append('prior_ML_sensmin|-')
open('jobs_hessian.txt', 'w').write('\n'.join(jobs_h) + '\n')
open('jobs_mcmc.txt', 'w').write('\n'.join(jobs_m) + '\n')
print(f'rgene_gamma = {RG}; {len(jobs_h)} Hessian jobs, {len(jobs_m)} MCMC jobs. Check calibration_nodes.tsv.')
