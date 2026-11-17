#!/usr/bin/env python3
"""Compare two unrooted trees on their shared tips: normalised RF distance and conflicting splits
that touch carnivores or key relatives. Usage: compare_trees.py treeA treeB out.tsv"""
import sys
from Bio import Phylo

def splits(tree, taxa):
    allt = frozenset(taxa); out = {}
    for cl in tree.get_nonterminals():
        s = frozenset(t.name for t in cl.get_terminals()) & allt
        if 1 < len(s) < len(allt) - 1:
            key = min(s, allt - s, key=lambda x: (len(x), sorted(x)))
            sup = cl.confidence if cl.confidence is not None else (cl.name or '')
            out[frozenset(key)] = sup
    return out

a, b = (Phylo.read(p, 'newick') for p in sys.argv[1:3])
shared = {t.name for t in a.get_terminals()} & {t.name for t in b.get_terminals()}
sa, sb = splits(a, shared), splits(b, shared)
only_a, only_b = set(sa) - set(sb), set(sb) - set(sa)
rf = len(only_a) + len(only_b); n = len(shared)
print(f'shared tips {n}; splits A {len(sa)} B {len(sb)}; RF {rf}; nRF {rf / (2 * (n - 3)):.4f}')
with open(sys.argv[3], 'w') as f:
    f.write('tree\tsize\tsupport\ttaxa\n')
    for tag, sset, sup in (('A_only', only_a, sa), ('B_only', only_b, sb)):
        for s in sorted(sset, key=len):
            f.write(f'{tag}\t{len(s)}\t{sup[s]}\t{",".join(sorted(s))[:400]}\n')
for tag, sset, sup in (('A_only', only_a, sa), ('B_only', only_b, sb)):
    for s in sorted(sset, key=len)[:25]:
        print(tag, len(s), sup[s], ','.join(sorted(s))[:160])
