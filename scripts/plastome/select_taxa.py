"""Select one plastome per species and subsample non-carnivorous relatives per family.

Nutritional mode is coded at genus level (species level for Brocchinia/Catopsis/Paepalanthus);
codes must be checked against the literature before analysis.
"""
import csv, random, re, sys
from collections import defaultdict, OrderedDict

random.seed(42)
EXCLUDE = {'CM045551.1', 'CM120794.1', 'ON005622.1',   # unannotated assembly records (no features)
           'PX309297.1', 'PP706316.1'}   # Chorizanthe novoana, Frankenia triandra: misplaced outside their family in the backbone tree (probable misidentification)
rows = [r for r in csv.DictReader(open('census_all.csv')) if r['accession'] not in EXCLUDE]
prev = {r['species']: r['genbank_accession'].split()[0] for r in csv.DictReader(open('/home/claude/work/repo/data/plastome_accessions.csv'))}

CARN = {'Nepenthes', 'Drosera', 'Dionaea', 'Aldrovanda', 'Drosophyllum', 'Utricularia', 'Genlisea', 'Pinguicula',
        'Sarracenia', 'Heliamphora', 'Darlingtonia', 'Roridula', 'Cephalotus', 'Byblis', 'Philcoxia'}
FACULT = {'Triphyophyllum'}
LOST = {'Habropetalum', 'Dioncophyllum', 'Ancistrocladus'}
CARN_SP = {'Brocchinia reducta', 'Brocchinia hechtioides', 'Catopsis berteroniana', 'Paepalanthus bromelioides'}
HOLO = {'Orobanche', 'Phelipanche', 'Aphyllon', 'Conopholis', 'Epifagus', 'Boschniakia', 'Kopsiopsis', 'Xylanche',
        'Cistanche', 'Christisonia', 'Aeginetia', 'Lathraea', 'Hyobanche', 'Mannagettaea', 'Phacellanthus', 'Gleadovia',
        'Platypholis', 'Diphelypaea', 'Harveya', 'Orobanchaceae'}
ORO_AUTO = {'Lindenbergia', 'Rehmannia', 'Triaenophora'}
MYCO = {'Monotropa', 'Monotropsis', 'Hypopitys', 'Monotropastrum', 'Allotropa', 'Pityopus', 'Pterospora', 'Sarcodes', 'Hemitomes', 'Cheilotheca'}
MIXO = {'Pyrola', 'Orthilia', 'Chimaphila', 'Moneses'}

QUOTA = {'Polygonaceae': 14, 'Plumbaginaceae': 8, 'Tamaricaceae': 6, 'Frankeniaceae': 9, 'Cactaceae': 12,
         'Martyniaceae': 9, 'Bignoniaceae': 8, 'Acanthaceae': 8, 'Lamiaceae': 8, 'Plantaginaceae': 10, 'Paulowniaceae': 3,
         'Actinidiaceae': 6, 'Clethraceae': 9, 'Elaeocarpaceae': 8, 'Brunelliaceae': 9, 'Oxalidaceae': 8, 'Rapateaceae': 9}


def binomial(org):
    org = re.sub(r'\s+(subsp\.|var\.|f\.)\s+.*$', '', org.strip())
    w = org.split()
    if len(w) >= 5 and w[2] == 'x':
        return ' '.join(w[:5])
    return ' '.join(w[:2])


def mode(sp, fam):
    g = sp.split()[0]
    if sp in CARN_SP or g in CARN: return 'carnivorous'
    if g in FACULT: return 'carnivorous (facultative)'
    if g in LOST: return 'carnivory lost'
    if fam == 'Orobanchaceae':
        if g in ORO_AUTO: return 'autotrophic'
        return 'holoparasitic' if g in HOLO else 'hemiparasitic'
    if g in MYCO: return 'mycoheterotrophic'
    if g in MIXO: return 'partially mycoheterotrophic'
    return 'autotrophic'


bysp = OrderedDict()
for r in rows:
    sp = binomial(r['organism'])
    if any(t in sp for t in (' sp.', ' cf.', ' aff.')) or len(sp.split()) < 2 or sp.split()[1][0].isupper():
        continue
    fam = r['group'] if r['group'] not in ('Droseraceae',) else 'Droseraceae'
    if r['query'] in ('Nepenthes',): fam = 'Nepenthaceae'
    bysp.setdefault((r['order'], fam, sp), []).append(r)


def pick(recs, sp):
    accs = [r['accession'] for r in recs]
    if sp in prev:
        for a in accs:
            if a.split('.')[0] == prev[sp]: return a
    nc = [a for a in accs if a.startswith('NC_')]
    return nc[0] if nc else sorted(accs)[-1]


cands = defaultdict(list)
for (order, fam, sp), recs in bysp.items():
    cands[(order, fam)].append((sp, pick(recs, sp), max(int(r['length']) for r in recs)))

sel = []
for (order, fam), lst in cands.items():
    random.seed(f'{fam}-42')
    keep = []
    special = [x for x in lst if mode(x[0], fam) != 'autotrophic']
    auto = [x for x in lst if mode(x[0], fam) == 'autotrophic']
    if fam == 'Orobanchaceae':
        # all holoparasites and autotrophs, hemiparasites round-robin by genus up to 20
        keep = [x for x in lst if mode(x[0], fam) in ('holoparasitic', 'autotrophic')]
        hemi = [x for x in lst if mode(x[0], fam) == 'hemiparasitic']
        bygen = defaultdict(list)
        for x in sorted(hemi): bygen[x[0].split()[0]].append(x)
        for v in bygen.values(): random.shuffle(v)
        while len([k for k in keep if mode(k[0], fam) == 'hemiparasitic']) < 20 and any(bygen.values()):
            for g in sorted(bygen):
                if bygen[g] and len([k for k in keep if mode(k[0], fam) == 'hemiparasitic']) < 20:
                    keep.append(bygen[g].pop())
    elif fam == 'Ericaceae':
        keep = special[:]
        bygen = defaultdict(list)
        for x in sorted(auto): bygen[x[0].split()[0]].append(x)
        for v in bygen.values(): random.shuffle(v)
        while len(keep) - len(special) < 14 and any(bygen.values()):
            for g in sorted(bygen):
                if bygen[g] and len(keep) - len(special) < 14: keep.append(bygen[g].pop())
    elif fam == 'Bromeliaceae':
        # all Brocchinia, Catopsis, Navia, Lindmania, Pitcairnia; 1 per other genus up to 12
        keep = [x for x in lst if x[0].split()[0] in ('Brocchinia', 'Catopsis', 'Navia', 'Lindmania', 'Pitcairnia')]
        bygen = defaultdict(list)
        for x in sorted(lst):
            if x not in keep: bygen[x[0].split()[0]].append(x)
        gens = sorted(bygen); random.shuffle(gens)
        for g in gens[:12]: keep.append(random.choice(bygen[g]))
    elif fam in QUOTA:
        bygen = defaultdict(list)
        for x in sorted(lst): bygen[x[0].split()[0]].append(x)
        for v in bygen.values(): random.shuffle(v)
        q = QUOTA[fam]
        while len(keep) < q and any(bygen.values()):
            for g in sorted(bygen):
                if bygen[g] and len(keep) < q: keep.append(bygen[g].pop())
    else:
        keep = lst
    for sp, acc, ln in keep:
        sel.append([order, fam, sp, acc, ln, mode(sp, fam)])

sel.sort(key=lambda r: (r[0], r[1], r[2]))
with open('taxa_selected.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['order', 'family', 'species', 'accession', 'length', 'nutrition']); w.writerows(sel)
from collections import Counter
print(len(sel), 'plastomes')
print(Counter(r[5] for r in sel))
print(Counter(r[0] for r in sel))
