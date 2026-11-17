mac_dating_v1 - Phase B: MCMCTree dating of the five topologies (295 taxa)

Put this folder NEXT TO mac_run_v7 (both in ~/Downloads). Start only after mac_run_v7 has finished
(it reads ../mac_run_v7/out/ml295_r1, c_SR, c_H4, c_H5, c_H6 .treefile; run mac_run_v7b first).

One-time setup (in Terminal):
  cd ~/Downloads/mac_dating_v1
  conda activate iqtree
  conda install -y -c conda-forge biopython numpy
  git clone --depth 1 https://github.com/abacus-gene/paml.git paml_src
  (cd paml_src/src && make baseml mcmctree)      # needs Xcode command-line tools; if 'make' is missing: xcode-select --install
  ./paml_src/src/mcmctree --version 2>&1 | head -2  # expect PAML 4.10.10 (the version the scripts were tested with)

Run:
  JOBS=12 nohup caffeinate -i bash run_dating.sh > run_dating.txt 2>&1 &

What it does:
  1. prep_dating.py   alignment split into codon positions 1+2 and 3 (two loci), rooted trees with calibrations C1-C5,
                      control files. CHECK calibration_nodes.tsv: every calibration must sit on the intended clade
                      (C3 = only Sarraceniaceae, Roridulaceae, Actinidiaceae, Clethraceae, Ericaceae; C4 sister = Actinidiaceae).
  2. Hessians         GTR+G4 for each topology and locus (baseml, approximate likelihood).
  3. MCMC (25 runs)   5 topologies x (prior-only + 2 posterior runs), independent-rates clock;
                      on the ML topology: autocorrelated clock, soft root, sensitivity minima, positions 1+2 only (2 runs each),
                      plus prior-only runs for the soft-root and sensitivity-minima calibration sets.
                      Each run: burn-in 20,000, then 10,000 samples every 20 iterations, on a 124-taxon subset (taxa_dating.txt).
  4. summarize_dating.py   results/node_ages.csv, node_ages_wide.csv (focal carnivore nodes, calibrated and
                      cross-validation nodes, with Fleischmann et al. 2018 ages for comparison), results/convergence.csv
                      (r1 vs r2: node-mean correlation, max difference, minimum ESS).

Expected time: Hessians about 2-5 h, MCMC about 4-8 h (overnight). Safe to re-run: finished steps are skipped.
Check progress:  cat run_dating.txt ; ls runs/*/.done | wc -l   (30 when everything is done)
Upload at the end:
  zip -r dating_out.zip results runs/*/FigTree.tre runs/*/out.txt runs/*/mcmctree.log logs calibration_nodes.tsv versions_dating.txt run_dating.txt
