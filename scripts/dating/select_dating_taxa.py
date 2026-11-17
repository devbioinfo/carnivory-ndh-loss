#!/usr/bin/env python3
"""Select the dating subset (~130 taxa) from the 295-taxon ML tree.
Rules: keep Amborella; keep every carnivore outside Lentibulariaceae and every lineage-defining taxon
(Tofieldiaceae, Byblidaceae, Sarraceniaceae, Roridulaceae, Cephalotaceae, Droseraceae, Nepenthaceae, Drosophyllaceae,
Dioncophyllaceae, Ancistrocladaceae); Lentibulariaceae: up to Q_LENT taxa per genus chosen to maximise phylogenetic spread
(so each genus crown is preserved); every other family x nutrition group: up to Q taxa maximising spread (the two most
distant tips first, so family crowns are preserved); plus calibration/cross-validation taxa."""
import csv,collections,sys
from Bio import Phylo
Q, Q_LENT = int(sys.argv[1]), int(sys.argv[2])
r=list(csv.DictReader(open('/home/claude/work/pan/taxa_phylo295.csv')))
nm=lambda s:s.replace(' x ','_x_').replace(' ','_').replace('.','')
info={nm(x['species']):x for x in r}
t=Phylo.read('/home/claude/work/pan/results/step1_v7/out/ml295_r1.treefile','newick'); t.root_with_outgroup('Amborella_trichopoda')
def pdmax(names,k):
    names=list(names)
    if len(names)<=k: return names
    best=max(((a,b) for i,a in enumerate(names) for b in names[i+1:]),key=lambda p:t.distance(p[0],p[1])); sel=list(best)
    while len(sel)<k: sel.append(max((n for n in names if n not in sel),key=lambda n:min(t.distance(n,s) for s in sel)))
    return sel
keep={'Amborella_trichopoda'}; groups=collections.defaultdict(list)
for n,x in info.items(): groups[(x['family'],x['nutrition'])].append(n)
ALL={'Tofieldiaceae','Byblidaceae','Sarraceniaceae','Roridulaceae','Cephalotaceae','Droseraceae','Nepenthaceae','Drosophyllaceae','Dioncophyllaceae','Ancistrocladaceae'}
for (f,nut),ns in groups.items():
    if f in ALL or ('carniv' in nut and f!='Lentibulariaceae'): keep|=set(ns)
    elif f=='Lentibulariaceae':
        g=collections.defaultdict(list)
        for n in ns: g[n.split('_')[0]].append(n)
        for gn in g.values(): keep|=set(pdmax(gn,Q_LENT))
    else: keep|=set(pdmax(ns,Q))
for p in ['Catalpa','Chilopsis','Elaeocarpus','Crinodendron','Clethra','Rhododendron','Actinidia','Brunellia']:
    c=sorted(x for x in info if x.startswith(p))
    if c and not keep&set(c): keep.add(c[0])
print(len(keep), dict(collections.Counter(info[n]['nutrition'] for n in keep)))
open('taxa_dating.txt','w').write('\n'.join(sorted(keep))+'\n')
