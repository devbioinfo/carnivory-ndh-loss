# Carnivorous plants lose the chloroplast NDH complex repeatedly, irreversibly and at parasite-like rates

Code, data and results for:

> Konhar, R. & Biswal, D. K. Carnivorous plants lose the chloroplast NDH complex repeatedly, irreversibly and at parasite-like rates. (manuscript in preparation for *Nature Ecology & Evolution*)

Corresponding author: Devendra Kumar Biswal (devbioinfo@gmail.com), North-Eastern Hill University, Shillong, India.

## What is here

| Folder | Contents |
| --- | --- |
| `data/plastomes/` | 311 taxa with GenBank accessions and plastome lengths; gene-by-gene calls (F functional, P pseudogene, − absent) for 79 plastid genes; evidence table behind each call |
| `data/traits/` | Nutritional strategy, CAM, aquatic, wet-habitat and epiphyte coding, 11 *ndh* calls per taxon, and literature sources for each trait |
| `data/nuclear/` | 36 query genes (32 NDH/cyclic-flow genes, 4 controls), *Arabidopsis* query proteins, and the 22-genome panel |
| `data/verification/` | Read-level support for *Triantha occidentalis*, *Byblis gigantea* transcript ORF calls, and the *Roridula gorgonias* junction-read test |
| `results/phylogeny/` | Concatenated alignments (295 and 311 taxa, gzipped), partition files, ML trees, constrained alternative topologies, AU/ELW tests, concordance factors, the 311-taxon placement tree |
| `results/dating/` | MCMCTree node ages for every run, convergence summary, calibration nodes, and the FigTree output of each chain |
| `results/traits/` | Chronogram used for trait analyses, stochastic-map loss events and rates, Pagel tests, phylogenetic logistic regression, jackknife, stratified tests, hidden-rates models, the 27 loss clades |
| `results/nuclear/` | Reciprocal-best-hit calls, genome tblastn confirmations and the combined nuclear matrix |
| `scripts/` | All analysis code, in pipeline order (see below) |
| `figures/` | Main Figures 1–6 and Extended Data Figures 1–10 (vector PDF, 183 mm wide) |

High-resolution TIFF versions of the figures are archived with the Zenodo release.

## Pipeline

1. **Plastome gene calls** — `scripts/plastome/`: `fetch.py` (GenBank download), `select_taxa.py`, `pangene.py` (79-gene F/P/− calls from annotation plus tblastn), `traits_v2.py`, `build_alignment.py`.
2. **Phylogeny** — `scripts/phylogeny/iqtree_commands_v7*.txt`: IQ-TREE 3.1.4, three codon-position partitions under GTR+F+Γ4, 1,000 ultrafast bootstraps and 1,000 SH-aLRT replicates; constrained topologies; AU test with 10,000 RELL replicates; gCF/sCF.
3. **Dating** — `scripts/dating/`: `select_dating_taxa.py`, `prep_dating.py` (124-taxon subset, five calibrations, control files), `run_dating.sh` (Hessians in baseml under GTR+Γ4, then MCMCTree), `recover.sh`, `summarize_dating.py`; clock comparison with `clockbf.sh` and `bf_summary.py` (power posteriors, 8 Gauss–Legendre points). PAML 4.10.10.
4. **Trait models** — `scripts/traits/` (R 4, ape, phytools, phylolm): `make_calibs.py` then `phaseC.R` (chronogram, Mk models, 200 stochastic maps, loss rates, timing), `pagel.R`, `pagelQ.R`, `jackknife.R`, `confounder.R`, `hrm.R`, `pgls_len.R`.
5. **Nuclear survey** — `scripts/nuclear/`: `fetch_ath.py`, `rbh.py` (reciprocal best hits, BLAST+), `genome_check.py` and `gbatch.sh` (tblastn confirmation of absences).
6. **Figures** — `scripts/figures/`: `fig1.py`–`fig5.py`, `ed1.py`–`ed10.py` (`ed8.py` produces main Fig. 6), shared style in `figstyle.py`; `finalize.py` fits each figure to 183 × 247 mm and writes 600-dpi TIFFs.

Scripts were run from a project folder and contain absolute paths from the original environment; set them to your own paths before re-running.

## Software

IQ-TREE 3.1.4; PAML 4.10.10 (MCMCTree, baseml); R 4 with ape, phytools and phylolm; BLAST+; Python 3 with Biopython 1.88, NumPy, SciPy, pandas and matplotlib.

## Licence

Code: MIT licence (`LICENSE`). Data, results and figures: Creative Commons Attribution 4.0 International (CC BY 4.0). Plastome sequences are from GenBank and remain subject to their original terms.

## Citation

See `CITATION.cff`. The Zenodo DOI will be added here once the archive is minted.
