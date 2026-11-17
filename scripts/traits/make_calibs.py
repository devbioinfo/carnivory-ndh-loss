#!/usr/bin/env python3
"""Node ages from an MCMCTree FigTree.tre (124-taxon dating subset) -> secondary calibrations for the full tree.
Each subset node is identified by two tips that span it (one from each daughter clade), so the MRCA can be found
in the 295-taxon tree. Ages in Ma (mean, 95% HPD).  usage: make_calibs.py FigTree.tre out.csv"""
import sys,re,io,csv
from Bio import Phylo
def read_figtree(p):
    s=open(p).read(); i=s.index('(',s.index('UTREE')); s=s[i:s.index(';',i)+1]
    s=re.sub(r'\)\s*\[&95%(?:HPD)?=\{\s*([\d.eE+-]+)\s*,\s*([\d.eE+-]+)\s*\}\]',lambda m:f')HPD_{m.group(1)}_{m.group(2)}',s)
    s=re.sub(r'\[[^\]]*\]','',s); return Phylo.read(io.StringIO(s),'newick')
t=read_figtree(sys.argv[1]); w=csv.writer(open(sys.argv[2],'w',newline=''))
w.writerow(['tipA','tipB','n_tips','age_mean','hpd_low','hpd_high'])
for c in t.get_nonterminals():
    a=c.clades[0].get_terminals()[0].name; b=c.clades[1].get_terminals()[0].name
    age=t.distance(c,c.get_terminals()[0])*100
    lo,hi=(float(x)*100 for x in c.name.split('_')[1:]) if c.name and c.name.startswith('HPD_') else (age,age)
    w.writerow([a,b,len(c.get_terminals()),round(age,4),round(lo,4),round(hi,4)])
print('nodes:',len(t.get_nonterminals()))
