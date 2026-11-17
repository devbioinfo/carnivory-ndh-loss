#!/usr/bin/env python3
"""Trait table v2: implements the external scientific review (traits_scientific_review.xlsx)
and adds model-ready columns plus book-derived carnivory descriptors (Ellison & Adamec 2018).

Input : results/traits.csv (v1, left untouched)
Output: results/traits_v2.csv, results/trait_decisions_v2.csv, results/trait_references_v2.csv
Every changed cell is logged with the review decision and its status (ACCEPT / PROVISIONAL / CHECK),
so provisional codes can be switched off before the final analysis.
"""
import csv, pandas as pd

t = pd.read_csv('results/traits.csv')
genus = t.species.str.split().str[0]
log = []

def setv(mask, col, val, decision, status, source):
    for i in t.index[mask]:
        old = t.at[i, col] if col in t.columns else ''
        if str(old) != str(val):
            log.append(dict(species=t.at[i, 'species'], column=col, old=old, new=val,
                            decision=decision, status=status, source=source))
        t.at[i, col] = val

sp = lambda *names: t.species.isin(names)

# ---------------------------------------------------------------- nutrition (model state)
# nutrition_v2: autotrophic | carnivorous | hemiparasitic | holoparasitic | partial_mycoheterotroph | full_mycoheterotroph
NUT = {'autotrophic': 'autotrophic', 'autotrophic (non-parasitic Orobanchaceae)': 'autotrophic',
       'autotrophic (myrmecophyte)': 'autotrophic', 'carnivorous': 'carnivorous',
       'carnivorous (facultative, juvenile leaves)': 'carnivorous',
       'carnivorous (indirect digestion via mutualist bugs)': 'carnivorous',
       'carnivory lost': 'autotrophic', 'hemiparasitic': 'hemiparasitic', 'holoparasitic': 'holoparasitic',
       'mycoheterotrophic': 'full_mycoheterotroph', 'partially mycoheterotrophic': 'partial_mycoheterotroph'}
t['nutrition_v2'] = t.nutrition.map(NUT)
assert t.nutrition_v2.notna().all()
t['carnivory_terminal'] = (t.nutrition_v2 == 'carnivorous').astype(int)
t['review_status'] = 'unchanged'

# secondary loss: terminal 0, ancestry recorded separately, reconstructed on the tree
m = sp('Ancistrocladus tectorius', 'Dioncophyllum thollonii', 'Habropetalum dawei')
setv(m, 'nutrition_v2', 'autotrophic', 'RECODE FOR MODEL: terminal state non-carnivorous', 'ACCEPT', 'Heubl2006')
t['ancestral_carnivory_hypothesis'] = ''
t.loc[m, 'ancestral_carnivory_hypothesis'] = 'secondary loss inferred (Heubl2006); reconstruct on tree, topology-dependent'
t.loc[sp('Brocchinia tatei'), 'ancestral_carnivory_hypothesis'] = 'possible derivation from carnivorous ancestor (EllisonAdamec2018 s3.7); low confidence'
t.loc[m, 'review_status'] = 'ACCEPT'

# Roridula, Triphyophyllum: carnivorous in main model; subtype kept as metadata; sensitivity flags
setv(sp('Roridula gorgonias'), 'nutrition', 'carnivorous (indirect, bug-mediated)', 'RETAIN MAIN + SENSITIVITY', 'ACCEPT', 'Roridula_PMC2729641')
setv(sp('Triphyophyllum peltatum'), 'nutrition', 'carnivorous (facultative, transient)', 'RETAIN MAIN + SENSITIVITY', 'ACCEPT', 'Triphyophyllum_PMC10731798')
t['sens_drop_roridula'] = sp('Roridula gorgonias').astype(int)
t['sens_drop_triphyophyllum'] = sp('Triphyophyllum peltatum').astype(int)
t.loc[sp('Roridula gorgonias', 'Triphyophyllum peltatum'), 'review_status'] = 'ACCEPT'

# Brocchinia
setv(sp('Brocchinia acuminata'), 'nutrition', 'autotrophic (myrmecophyte)', 'ACCEPT: myrmecotrophy kept separate from carnivory', 'ACCEPT', 'Givnish1997')
setv(sp('Brocchinia acuminata'), 'nutrition_conf', 'high', 'ACCEPT', 'ACCEPT', 'Givnish1997')
setv(sp('Brocchinia tatei'), 'nutrition', 'autotrophic (tank epiphyte)', 'ACCEPT PROVISIONALLY: non-carnivorous tank epiphyte', 'PROVISIONAL', 'Givnish1997')
setv(sp('Brocchinia tatei'), 'nutrition_conf', 'medium', 'ACCEPT PROVISIONALLY', 'PROVISIONAL', 'Givnish1997')
setv(sp('Brocchinia tatei'), 'habitat', 'epiphyte', 'tank epiphyte (Givnish1997)', 'PROVISIONAL', 'Givnish1997')
t.loc[sp('Brocchinia acuminata'), 'review_status'] = 'ACCEPT'
t.loc[sp('Brocchinia tatei'), 'review_status'] = 'PROVISIONAL'

# Lathraea: accept holoparasitic
setv(sp('Lathraea japonica', 'Lathraea squamaria'), 'nutrition_conf', 'high', 'ACCEPT', 'ACCEPT', 'Lathraea_PMC12724003')
t.loc[sp('Lathraea japonica', 'Lathraea squamaria'), 'review_status'] = 'ACCEPT'

# Pyrola: split evidence
t['nutrition_evidence_level'] = 'genus/clade'
t.loc[t.nutrition.isin(['autotrophic', 'autotrophic (non-parasitic Orobanchaceae)']), 'nutrition_evidence_level'] = 'default'
setv(sp('Pyrola rotundifolia'), 'nutrition_conf', 'high', 'SPLIT EVIDENCE: direct isotopic support', 'ACCEPT', 'Tedersoo2007')
t.loc[sp('Pyrola rotundifolia'), 'nutrition_evidence_level'] = 'species (isotopes)'
setv(sp('Pyrola atropurpurea', 'Pyrola decorata'), 'nutrition_conf', 'check', 'SPLIT EVIDENCE: genus-inferred, verify species isotopes', 'PROVISIONAL', 'Tedersoo2007')
t.loc[sp('Pyrola atropurpurea', 'Pyrola decorata'), 'review_status'] = 'PROVISIONAL'
t.loc[sp('Pyrola rotundifolia'), 'review_status'] = 'ACCEPT'
t['sens_pyrola_autotrophic'] = sp('Pyrola atropurpurea', 'Pyrola decorata').astype(int)

# ---------------------------------------------------------------- habitat (not collapsed)
setv(sp('Utricularia cucullata', 'Utricularia olivacea', 'Utricularia tenuicaulis'), 'habitat_conf', 'high', 'ACCEPT aquatic', 'ACCEPT', 'review_WFO/NCparks/Tenuicaulis2022')
t.loc[sp('Utricularia cucullata', 'Utricularia olivacea', 'Utricularia tenuicaulis'), 'review_status'] = 'ACCEPT'
setv(sp('Utricularia nephrophylla'), 'habitat_conf', 'high', 'ACCEPT terrestrial/lithophytic', 'ACCEPT', 'POWO')
t.loc[sp('Utricularia nephrophylla'), 'review_status'] = 'ACCEPT'
setv(sp('Utricularia longifolia'), 'habitat_conf', 'medium', 'PROVISIONAL terrestrial/lithophytic (mixed habit)', 'PROVISIONAL', 'POWO')
t.loc[sp('Utricularia longifolia'), 'review_status'] = 'PROVISIONAL'
setv(sp('Deinostema violacea'), 'habitat_conf', 'check', 'PROVISIONAL wetland; seek specimen/flora account', 'CHECK', 'review')
t.loc[sp('Deinostema violacea'), 'review_status'] = 'CHECK'

# primary habitat + full habitat set + moisture/salinity, per protocol (no default collapsing)
HAB = {  # v1 string -> (primary, all, moisture)
    'terrestrial': ('terrestrial', 'terrestrial', 'mesic'),
    'terrestrial (wet)': ('terrestrial', 'terrestrial', 'wet'),
    'terrestrial (dry)': ('terrestrial', 'terrestrial', 'dry'),
    'terrestrial (wetland)': ('wetland', 'wetland;terrestrial', 'wet'),
    'terrestrial/climber': ('terrestrial', 'terrestrial', 'wet'),
    'arid': ('terrestrial', 'terrestrial', 'arid'),
    'arid (host roots)': ('terrestrial', 'terrestrial', 'arid'),
    'arid/halophytic': ('terrestrial', 'terrestrial', 'arid;saline'),
    'halophytic': ('terrestrial', 'terrestrial', 'saline'),
    'mangrove': ('wetland', 'wetland', 'saline;inundated'),
    'riparian': ('terrestrial', 'terrestrial', 'wet'),
    'aquatic': ('aquatic', 'aquatic', 'submerged'),
    'aquatic/wetland': ('wetland', 'wetland;aquatic', 'wet'),
    'aquatic/wetland (P. amphibia amphibious)': ('wetland', 'wetland;aquatic', 'amphibious'),
    'wetland': ('wetland', 'wetland', 'wet'),
    'lithophyte': ('lithophytic', 'lithophytic', 'wet'),
    'lithophyte/terrestrial': ('lithophytic', 'lithophytic;terrestrial', 'wet'),
    'terrestrial/lithophyte': ('terrestrial', 'terrestrial;lithophytic', 'wet'),
    'epiphyte': ('epiphytic', 'epiphytic', 'mesic'),
    'aquatic epiphyte (bromeliad tanks)': ('aquatic', 'aquatic;epiphytic', 'submerged'),
    'rheophyte': ('rheophytic', 'rheophytic;aquatic', 'submerged'),
}
assert set(t.habitat) <= set(HAB), set(t.habitat) - set(HAB)
t['habitat_primary'] = t.habitat.map(lambda h: HAB[h][0])
t['habitat_all'] = t.habitat.map(lambda h: HAB[h][1])
t['moisture'] = t.habitat.map(lambda h: HAB[h][2])
t.loc[sp('Utricularia longifolia', 'Utricularia nephrophylla'), ['habitat_primary', 'habitat_all']] = ['terrestrial', 'terrestrial;lithophytic']
# binary used by the NDH hypothesis: foliage submerged/floating in water
t['submerged_foliage'] = t.moisture.str.contains('submerged').astype(int)

# ---------------------------------------------------------------- photosynthetic pathway
# CAM_primary: species-level or clade-constitutive evidence only; CAM_sensitivity: adds genus-derived calls
t['CAM_evidence'] = 'default C3'
t.loc[t.photosynthesis == 'non-photosynthetic', 'CAM_evidence'] = 'non-photosynthetic'
cact = t.family.eq('Cactaceae') & t.photosynthesis.eq('CAM')
t.loc[cact, 'CAM_evidence'] = 'clade-constitutive (Cactoideae/Opuntioideae); approve for primary'
t.loc[sp('Leuenbergeria portulacifolia'), 'CAM_evidence'] = 'C3 with weak CAM (Edwards2008)'
brom = t.family.eq('Bromeliaceae') & t.photosynthesis.str.contains('CAM')
t.loc[brom, 'CAM_evidence'] = 'genus-derived (Crayn2015); verify species d13C in Crayn2015 appendix'
t['CAM_primary'] = 'C3'
t.loc[t.photosynthesis == 'non-photosynthetic', 'CAM_primary'] = 'NA'
t.loc[cact, 'CAM_primary'] = 'CAM'
t.loc[brom, 'CAM_primary'] = 'unknown'
t['CAM_sensitivity'] = t.CAM_primary
t.loc[brom & t.photosynthesis.eq('CAM'), 'CAM_sensitivity'] = 'CAM'
for s in ['Tillandsia rubroviolacea']:
    setv(sp(s), 'photosynthesis', 'unknown (CAM vs C3)', 'LEAVE UNKNOWN', 'ACCEPT', 'Crayn2015')
    t.loc[sp(s), ['CAM_primary', 'CAM_sensitivity', 'review_status']] = ['unknown', 'unknown', 'ACCEPT']
setv(sp('Pseudalcantarea viridiflora'), 'photosynthesis_conf', 'medium', 'PROVISIONAL C3; verify species measurement/synonymy', 'PROVISIONAL', 'Crayn2015')
t.loc[sp('Pseudalcantarea viridiflora'), 'review_status'] = 'PROVISIONAL'
t.loc[sp('Pseudalcantarea viridiflora'), 'CAM_evidence'] = 'C3 provisional (Tillandsioideae CAM confined to Tillandsia; Crayn2015)'

# ---------------------------------------------------------------- carnivory descriptors (metadata; EllisonAdamec2018)
TRAP = {  # genus -> (trap_type, motility, continuity, trap_photosynthesis_coupling)
    'Nepenthes': ('pitfall', 'passive', 'constitutive', 'partly decoupled (lamina photosynthesises; pitcher low Amax)'),
    'Sarracenia': ('pitfall', 'passive', 'constitutive', 'coupled (pitcher is the main leaf)'),
    'Darlingtonia': ('pitfall', 'passive', 'constitutive', 'coupled (pitcher is the main leaf)'),
    'Heliamphora': ('pitfall', 'passive', 'constitutive', 'coupled (pitcher is the main leaf)'),
    'Cephalotus': ('pitfall', 'passive', 'constitutive', 'decoupled (separate flat leaves and pitchers)'),
    'Brocchinia': ('pitfall (tank)', 'passive', 'constitutive', 'coupled (rosette leaves form the tank)'),
    'Drosera': ('flypaper', 'active (slow)', 'constitutive', 'coupled'),
    'Drosophyllum': ('flypaper', 'passive', 'constitutive', 'coupled'),
    'Pinguicula': ('flypaper', 'active (slow)', 'constitutive (some seasonal heterophylly)', 'coupled'),
    'Roridula': ('flypaper, indirect (bug-mediated)', 'passive', 'constitutive', 'coupled'),
    'Triphyophyllum': ('flypaper', 'passive', 'facultative, transient', 'decoupled (separate trap leaves, limited life stage)'),
    'Dionaea': ('snap', 'active (rapid)', 'constitutive', 'partly decoupled (photosynthetic petiole)'),
    'Aldrovanda': ('snap', 'active (rapid)', 'constitutive', 'coupled'),
    'Utricularia': ('suction', 'active (rapid)', 'constitutive', 'decoupled (bladders separate from photosynthetic shoots)'),
    'Genlisea': ('eel trap', 'passive', 'constitutive', 'decoupled (subterranean rhizophylls)'),
}
for c in ['trap_type', 'trap_motility', 'carnivory_continuity', 'trap_photosynthesis_coupling']:
    t[c] = ''
car = t.carnivory_terminal.eq(1)
for i in t.index[car]:
    t.loc[i, ['trap_type', 'trap_motility', 'carnivory_continuity', 'trap_photosynthesis_coupling']] = TRAP[genus[i]]
# ordinal degree of carnivory (provisional, for exploratory models only):
# 0 none; 1 indirect or transient; 2 constitutive, passive or slow-moving; 3 constitutive, rapid motile trap
def degree(i):
    if not car[i]: return 0
    tt, mo, co, _ = TRAP[genus[i]]
    if 'indirect' in tt or 'transient' in co: return 1
    return 3 if 'rapid' in mo else 2
t['carnivory_degree'] = [degree(i) for i in t.index]

BROC = {'Brocchinia reducta': 'impoundment carnivory', 'Brocchinia hechtioides': 'impoundment carnivory',
        'Brocchinia acuminata': 'myrmecotrophy (ant-fed)', 'Brocchinia tatei': 'tank epiphyte (litter/canopy debris)',
        'Brocchinia micrantha': 'tank, litter impoundment', 'Brocchinia paniculata': 'non-impounding, root uptake'}
t['brocchinia_nutrient_mode'] = t.species.map(BROC).fillna('')

# ---------------------------------------------------------------- write
t['overall_conf'] = t[['nutrition_conf', 'habitat_conf', 'photosynthesis_conf']].apply(
    lambda r: 'check' if 'check' in r.values else ('medium' if 'medium' in r.values else 'high'), axis=1)
front = ['order', 'family', 'species', 'accession', 'nutrition_v2', 'carnivory_terminal', 'habitat_primary',
         'habitat_all', 'moisture', 'submerged_foliage', 'CAM_primary', 'CAM_sensitivity', 'review_status', 'overall_conf']
rest = [c for c in t.columns if c not in front]
t[front + rest].to_csv('results/traits_v2.csv', index=False)
pd.DataFrame(log).to_csv('results/trait_decisions_v2.csv', index=False)

refs = list(csv.DictReader(open('results/trait_references.csv')))
extra = {
    'Givnish1997': 'Givnish TJ et al. (1997) Molecular evolution, adaptive radiation and geographic speciation in Brocchinia ... In: Givnish & Sytsma (eds) Molecular Evolution and Adaptive Radiation. CUP. https://sytsma.botany.wisc.edu/pdf/Givnish%26al1997.pdf',
    'Roridula_PMC2729641': 'Roridula prey-derived nutrient uptake (as cited in external review). https://pmc.ncbi.nlm.nih.gov/articles/PMC2729641/',
    'Triphyophyllum_PMC10731798': 'Triphyophyllum transient carnivory (as cited in external review). https://pmc.ncbi.nlm.nih.gov/articles/PMC10731798/',
    'Lathraea_PMC12724003': 'Lathraea holoparasitism (as cited in external review). https://pmc.ncbi.nlm.nih.gov/articles/PMC12724003/',
    'POWO': 'Plants of the World Online, Royal Botanic Gardens Kew. https://powo.science.kew.org/',
    'review_WFO/NCparks/Tenuicaulis2022': 'World Flora Online wfo-0000415747; NC Parks flora (U. olivacea); Mitochondrial DNA B 2022 doi:10.1080/23802359.2022.2080597',
    'Petersen2025': 'Petersen G (2025) Repeated loss of plastid NDH during evolution of land plants. Annals of Botany. doi:10.1093/aob/mcaf198',
}
keys = {r[list(r)[0]] for r in refs}
with open('results/trait_references_v2.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(list(refs[0]))
    for r in refs: w.writerow(list(r.values()))
    for k, v in extra.items():
        if k not in keys: w.writerow([k, v])

print(t.nutrition_v2.value_counts().to_string(), '\n')
print(t.habitat_primary.value_counts().to_string(), '\n')
print(t.CAM_primary.value_counts().to_string(), '\n')
print(t.review_status.value_counts().to_string(), '\n')
print(t[car].groupby(['carnivory_degree', 'trap_photosynthesis_coupling']).size().to_string())
print(len(log), 'logged changes')
