#!/usr/bin/env python3
"""Summarise MCMCTree runs: focal node ages (mean, 95% HPD), convergence of replicate runs, prior checks.
Writes results/node_ages.csv, results/convergence.csv, results/node_ages_wide.csv."""
import os, re, csv, glob, io, math
from Bio import Phylo
import numpy as np
os.makedirs('results', exist_ok=True)
fam = {}; genus = {}
for r in list(csv.DictReader(open('taxa_phylo295.csv'))):
    n = r['species'].replace(' x ', '_x_').replace(' ', '_').replace('.', '')
    fam[n] = r['family']; genus[n] = n.split('_')[0]
def by(f=None, g=None, sp=None):
    return [n for n in fam if (f and fam[n] in f) or (g and genus[n] in g) or (sp and n in sp)]
# name: (tips, which, validation age from Fleischmann et al. 2018 / calibration, note)
FOCAL = {
 'root (C1)': (list(fam), 'crown', '136-139.35 prior', ''),
 'crown mesangiosperms (C2)': ([n for n in fam if not n.startswith('Amborella')], 'crown', '>=123.0 prior', ''),
 'crown core Ericales (C3)': (by(f={'Sarraceniaceae','Roridulaceae','Actinidiaceae','Clethraceae','Ericaceae'}), 'crown', '>=85.7 prior', ''),
 'stem Roridula (C4)': (by(g={'Roridula'}), 'stem', '>=35.0 prior; V 38.1', ''),
 'stem Aldrovanda (C5)': (by(g={'Aldrovanda'}), 'stem', '>=33.9 prior', ''),
 'crown carnivorous Caryophyllales': (by(f={'Nepenthaceae','Droseraceae','Drosophyllaceae','Dioncophyllaceae'}), 'crown', 'V 84.8 (Droseraceae/Nepenthaceae stem)', ''),
 'stem Droseraceae': (by(f={'Droseraceae'}), 'stem', '', ''),
 'crown Droseraceae (X1)': (by(f={'Droseraceae'}), 'crown', 'X1 >=33.9', ''),
 'crown Drosera': (by(g={'Drosera'}), 'crown', 'V 53.4', ''),
 'Dionaea+Aldrovanda': (by(g={'Dionaea','Aldrovanda'}), 'crown', 'V 48.0', ''),
 'stem Nepenthes': (by(g={'Nepenthes'}), 'stem', 'X2 (Miocene pollen)', ''),
 'crown Nepenthes': (by(g={'Nepenthes'}), 'crown', 'V ~8.7', ''),
 'stem Drosophyllum': (by(g={'Drosophyllum'}), 'stem', 'V 70.4', ''),
 'stem Triphyophyllum': (by(g={'Triphyophyllum'}), 'stem', 'V 6.9', ''),
 'stem Sarraceniaceae': (by(f={'Sarraceniaceae'}), 'stem', '', ''),
 'crown Sarraceniaceae': (by(f={'Sarraceniaceae'}), 'crown', 'V 48.6', ''),
 'stem Cephalotus': (by(g={'Cephalotus'}), 'stem', 'V 32.2', ''),
 'Elaeocarpus+Crinodendron (X4)': (by(g={'Elaeocarpus','Crinodendron'}), 'crown', 'X4 (horizon uncertain)', ''),
 'stem Byblis': (by(g={'Byblis'}), 'stem', '', ''),
 'stem Lentibulariaceae': (by(f={'Lentibulariaceae'}), 'stem', '', ''),
 'crown Lentibulariaceae': (by(f={'Lentibulariaceae'}), 'crown', '', ''),
 'crown Pinguicula': (by(g={'Pinguicula'}), 'crown', '', ''),
 'crown Utricularia': (by(g={'Utricularia'}), 'crown', 'V 31.0', ''),
 'crown Genlisea': (by(g={'Genlisea'}), 'crown', '', ''),
 'Catalpa+Chilopsis (X3)': (by(g={'Catalpa','Chilopsis'}), 'crown', 'X3 >=33.9', ''),
 'stem Brocchinia reducta+hechtioides': (by(sp={'Brocchinia_reducta','Brocchinia_hechtioides'}), 'stem', '', ''),
 'crown Brocchinia reducta+hechtioides': (by(sp={'Brocchinia_reducta','Brocchinia_hechtioides'}), 'crown', 'V 1.9 (B. reducta)', ''),
 'stem Catopsis': (by(g={'Catopsis'}), 'stem', '', ''),
 'crown Catopsis': (by(g={'Catopsis'}), 'crown', '', ''),
 'stem Triantha': (by(g={'Triantha'}), 'stem', '', ''),
 'crown Tofieldiaceae': (by(f={'Tofieldiaceae'}), 'crown', '', ''),
}
def read_figtree(p):
    s = open(p).read(); i = s.index('(', s.index('UTREE')); s = s[i:s.index(';', i) + 1]
    s = re.sub(r'\)\s*\[&95%(?:HPD)?=\{\s*([\d.eE+-]+)\s*,\s*([\d.eE+-]+)\s*\}\]', lambda m: f')HPD_{m.group(1)}_{m.group(2)}', s)
    s = re.sub(r'\[[^\]]*\]', '', s)
    return Phylo.read(io.StringIO(s), 'newick')
def node_age(t, c):
    tip = c.get_terminals()[0]; return t.distance(c, tip)
rows, wide = [], {}
for d in sorted(glob.glob('runs/*/FigTree.tre')):
    run = d.split('/')[1]; t = read_figtree(d)
    names = {x.name for x in t.get_terminals()}
    for k, (tips, which, val, _) in FOCAL.items():
        tips = [x for x in tips if x in names]
        if not tips: continue
        if len(tips) == 1: crown = next(x for x in t.get_terminals() if x.name == tips[0])
        else: crown = t.common_ancestor(tips)
        if which == 'crown' and len(tips) == 1: continue
        node = crown
        if which == 'stem':
            path = t.get_path(crown); node = path[-2] if len(path) > 1 else t.root
        age = node_age(t, node) * 100
        lo = hi = ''
        if node.name and node.name.startswith('HPD_'):
            _, a, b = node.name.split('_'); lo, hi = float(a) * 100, float(b) * 100
        mono = len(crown.get_terminals()) == len(tips)
        rows.append([run, k, round(age, 2), round(lo, 2) if lo != '' else '', round(hi, 2) if hi != '' else '', val, '' if mono else 'group not monophyletic in this tree'])
        wide.setdefault(k, {})[run] = f'{age:.1f} ({lo:.1f}-{hi:.1f})' if lo != '' else f'{age:.1f}'
with open('results/node_ages.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['run', 'node', 'mean_Ma', 'hpd95_low', 'hpd95_high', 'calibration_or_validation', 'note']); w.writerows(rows)
runs = sorted({r[0] for r in rows})
with open('results/node_ages_wide.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['node'] + runs)
    for k in FOCAL:
        if k in wide: w.writerow([k] + [wide[k].get(r, '') for r in runs])
# convergence of replicate pairs
def load(p):
    with open(p) as fh:
        h = fh.readline().split('\t'); idx = [i for i, x in enumerate(h) if x.startswith('t_n')]
        a = np.array([[float(l.split('\t')[i]) for i in idx] for l in fh if l.strip()])
    return a
def ess(x):
    x = x - x.mean(); n = len(x); v = x.var()
    if v == 0: return float(n)
    f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n] / (v * n)
    s = 0.0
    for k in range(1, n):
        if ac[k] < 0.05: break
        s += ac[k]
    return n / (1 + 2 * s)
conv = []
for r1 in sorted(glob.glob('runs/*_r1/mcmc.txt')):
    r2 = r1.replace('_r1/', '_r2/')
    if not os.path.exists(r2): continue
    a, b = load(r1), load(r2)
    ma, mb = a.mean(0), b.mean(0)
    rel = np.abs(ma - mb) / ((ma + mb) / 2)
    ea = min(ess(a[:, j]) for j in range(a.shape[1])); eb = min(ess(b[:, j]) for j in range(b.shape[1]))
    esa = np.array([ess(a[:, j]) for j in range(a.shape[1])]); esb = np.array([ess(b[:, j]) for j in range(b.shape[1])])
    esc = esa + esb   # ESS of the two chains combined
    big = (ma + mb) / 2 >= 0.05   # nodes older than 5 Myr (relative differences of very young nodes are not informative)
    conv.append([r1.split('/')[1][:-3], a.shape[0], b.shape[0], round(float(np.corrcoef(ma, mb)[0, 1]), 5),
                 round(float(rel[big].max()) * 100, 2), round(float(np.median(rel)) * 100, 3), round(ea), round(eb),
                 round(float(np.median(esc))), round(float((esc >= 200).mean()) * 100, 1),
                 'OK' if rel[big].max() < 0.05 and np.median(esc) >= 200 and (esc >= 200).mean() >= 0.9 else 'CHECK'])
with open('results/convergence.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['analysis', 'samples_r1', 'samples_r2', 'corr_node_means', 'max_rel_diff_%_nodes_>5Myr', 'median_rel_diff_%', 'min_ESS_r1', 'min_ESS_r2', 'median_ESS_combined', 'pct_nodes_ESS>=200_combined', 'verdict']); w.writerows(conv)
for c in conv: print(*c)
print('wrote results/node_ages.csv, node_ages_wide.csv, convergence.csv')
