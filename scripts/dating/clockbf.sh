#!/usr/bin/env bash
# Clock-model comparison (independent vs autocorrelated rates) on the ML topology by Bayes factors.
# Marginal likelihoods by Gauss-Legendre quadrature over 8 power-posterior steps (BayesFactorBeta),
# as in dos Reis et al. 2018 (Syst. Biol.) / mcmc3r. 16 chains, reusing the ML Hessian (approximate likelihood).
# Run inside ~/Downloads/mac_dating_v1:
#   JOBS=16 nohup nice -n 10 caffeinate -i bash clockbf.sh > clockbf.txt 2>&1 & disown
# Safe to re-run after an interruption: finished chains are skipped, unfinished ones resume from their checkpoint.
cd "$(dirname "$0")"
JOBS=${JOBS:-16}; NS=${NS:-5000}; SF=${SF:-100}
[ -x paml_src/src/mcmctree ] && export PATH="$PWD/paml_src/src:$PATH"
MCMCTREE=$(command -v mcmctree) || { echo "mcmctree not found"; exit 1; }; export MCMCTREE
[ -s runs/hessian_ML/in.BV ] || { echo "runs/hessian_ML/in.BV missing"; exit 1; }
# Gauss-Legendre nodes on [-1,1] (n = 8); beta = (x + 1) / 2
X=(-0.9602898565 -0.7966664774 -0.5255324099 -0.1834346425 0.1834346425 0.5255324099 0.7966664774 0.9602898565)
mkdir -p runs_bf; : > runs_bf/jobs.txt
for clk in 2 3; do name=$([ $clk = 2 ] && echo ILN || echo GBM)
  for i in 0 1 2 3 4 5 6 7; do
    d=runs_bf/${name}_b$i; b=$(python3 -c "print(f'{(${X[$i]}+1)/2:.10f}')")
    mkdir -p "$d"; echo "$d" >> runs_bf/jobs.txt
    [ -e "$d/.done" ] && continue
    cp runs/hessian_ML/in.BV "$d/in.BV"
    perl -pe "s/^\s*seed\s*=.*/          seed = $((7000 + 100*clk + i))/; s/^\s*clock\s*=.*/         clock = $clk/;
              s/^\s*nsample\s*=.*/       nsample = $NS/; s/^\s*sampfreq\s*=.*/      sampfreq = $SF/; s/^\s*checkpoint\s*=.*//;
              s/^\s*BayesFactorBeta\s*=.*//" runs/post_ML_r1/mcmctree.ctl > "$d/mcmctree.ctl"
    echo "BayesFactorBeta = $b" >> "$d/mcmctree.ctl"
    if [ -s "$d/mcmctree.ckpt" ]; then
      n=1; while [ -e "$d/mcmc_part$n.txt" ]; do n=$((n+1)); done
      [ -s "$d/mcmc.txt" ] && awk -F'\t' 'NR==1{nf=NF; print; next} NF==nf' "$d/mcmc.txt" > "$d/mcmc_part$n.txt"
      have=0; for f in "$d"/mcmc_part*.txt; do [ -e "$f" ] && have=$((have + $(wc -l < "$f") - 1)); done
      left=$((NS - have)); [ $left -lt 100 ] && left=100
      perl -pi -e "s/^\s*nsample\s*=.*/       nsample = $left/" "$d/mcmctree.ctl"
      echo "checkpoint = 2 0.05 mcmctree.ckpt" >> "$d/mcmctree.ctl"
    else echo "checkpoint = 1 0.05 mcmctree.ckpt" >> "$d/mcmctree.ctl"; fi
  done
done
runbf () { local d=$1; [ -e "$d/.done" ] && return 0
  echo "[$(date '+%H:%M')] start $d"
  ( cd "$d" && "$MCMCTREE" mcmctree.ctl > mcmctree.log 2>&1 ) && touch "$d/.done" && echo "[$(date '+%H:%M')] done  $d" || echo "[$(date '+%H:%M')] FAILED $d"; }
export -f runbf
xargs -P "$JOBS" -I{} bash -c 'runbf "$@"' _ {} < runs_bf/jobs.txt
python3 bf_summary.py | tee runs_bf/clock_bayes_factor.txt
echo "All done $(date '+%H:%M')."
