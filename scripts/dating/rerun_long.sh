#!/usr/bin/env bash
# Rerun the 18 posterior chains 10x longer (sample every 200 iterations: 2.02 M iterations per chain).
# The first set (every 20) gave consistent ages between replicate runs but too few effective samples (min ESS 8-37).
# Prior-only runs and Hessians are kept. Run inside mac_dating_v1:
#   JOBS=12 nohup caffeinate -i bash rerun_long.sh > rerun_long.txt 2>&1 &
cd "$(dirname "$0")"
mkdir -p first_short_runs
for d in runs/post_* runs/sens_*; do
  r=$(basename "$d")
  [ -f "$d/FigTree.tre" ] && mkdir -p "first_short_runs/$r" && cp "$d/FigTree.tre" "first_short_runs/$r/" 2>/dev/null
  rm -f "$d/.done" "$d/mcmc.txt" "$d/FigTree.tre" "$d/out.txt" "$d/mcmctree.log" "$d"/SeedUsed "$d"/rates.txt "$d"/r1.txt
  perl -pi -e 's/^\s*sampfreq\s*=.*/      sampfreq = 200/' "$d/mcmctree.ctl"
done
cp results/convergence.csv first_short_runs/convergence_short.csv 2>/dev/null
cp results/node_ages.csv first_short_runs/node_ages_short.csv 2>/dev/null
echo "chains reset: $(ls -d runs/post_* runs/sens_* | wc -l)"
bash run_dating.sh
