#!/usr/bin/env bash
# Phase B: MCMCTree dating of the five topologies (approximate likelihood, GTR+G Hessians).
# Run inside this folder after mac_run_v7 has finished (needs ../mac_run_v7/out/*.treefile):
#   conda activate iqtree
#   JOBS=12 caffeinate -i bash run_dating.sh
# Safe to re-run: finished steps are skipped (runs/<run>/.done).
set -uo pipefail
JOBS=${JOBS:-12}
[ -x paml_src/src/mcmctree ] && export PATH="$PWD/paml_src/src:$PATH"   # PAML built from source (README step 2)
MCMCTREE=$(command -v mcmctree) || { echo "mcmctree not found: conda install -c bioconda paml"; exit 1; }
BASEML=$(command -v baseml)     || { echo "baseml not found"; exit 1; }
export MCMCTREE BASEML
mkdir -p logs
{ echo "mcmctree: $MCMCTREE"; "$MCMCTREE" 2>&1 </dev/null | head -3; python3 -c "import Bio;print('biopython',Bio.__version__)"; } > versions_dating.txt 2>&1

if [ ! -e .prep.done ]; then
  python3 prep_dating.py | tee logs/prep.log || exit 1
  touch .prep.done
fi

hessian () { local run=$1; local d="runs/$run"
  [ -e "$d/.done" ] && { echo "[skip] $run"; return 0; }
  echo "[$(date '+%H:%M')] start $run"
  ( cd "$d" || exit 1
    # 1) let mcmctree write the per-locus baseml control files (baseml hidden so it does not run them with HKY)
    PATH=/usr/bin:/bin "$MCMCTREE" mcmctree.ctl > mcmctree_usedata3.log 2>&1
    ls tmp0*.ctl >/dev/null 2>&1 || { echo "no tmp*.ctl written"; exit 1; }
    # 2) GTR+G4 Hessian for each locus
    # each locus in its own subfolder so the two baseml runs can go in parallel (both write rst2)
    for f in tmp0*.ctl; do
      b=${f%.ctl}; mkdir -p "$b"; cp "$f" "$b.txt" "$b.trees" "$b/"
      perl -pi -e 's/^\s*model\s*=.*/model = 7/; s/^\s*ncatG\s*=.*/ncatG = 4/' "$b/$f"
      grep -q "getSE" "$b/$f" || echo "getSE = 2" >> "$b/$f"
      # baseml writes the Hessian to rst2 BEFORE inverting it for standard errors; with zero-length branches
      # (identical sequences at codon positions 1+2) that inversion fails and baseml exits non-zero.
      # MCMCTree needs only rst2, so accept it whenever it holds a complete Hessian.
      ( cd "$b" && { "$BASEML" "$f" > baseml.log 2>&1; grep -q "Hessian" rst2 2>/dev/null && cp rst2 "../$b.rst2"; } ) &
    done
    wait
    for f in tmp0*.ctl; do [ -s "${f%.ctl}.rst2" ] || { echo "baseml failed for $f"; exit 1; }; done
    cat tmp0*.rst2 > in.BV
    cp tmp0001.rst2 in_pos12.BV ) > "logs/$run.log" 2>&1 \
  && { touch "$d/.done"; echo "[$(date '+%H:%M')] done  $run"; } \
  || echo "[$(date '+%H:%M')] FAILED $run (logs/$run.log)"; }

mcmc () { local spec=$1; local run=${spec%%|*} dep=${spec#*|}; local d="runs/$run"
  [ -e "$d/.done" ] && { echo "[skip] $run"; return 0; }
  if [ "$dep" != "-" ]; then
    local h=${dep%%:*} bv=in.BV; [ "$dep" != "$h" ] && bv=in_pos12.BV
    [ -e "runs/$h/.done" ] || { echo "[wait] $run needs $h (not done)"; return 1; }
    cp "runs/$h/$bv" "$d/in.BV"
  fi
  echo "[$(date '+%H:%M')] start $run"
  ( cd "$d" && "$MCMCTREE" mcmctree.ctl > mcmctree.log 2>&1 ) && [ -s "$d/FigTree.tre" ] \
    && { touch "$d/.done"; echo "[$(date '+%H:%M')] done  $run"; } \
    || echo "[$(date '+%H:%M')] FAILED $run (runs/$run/mcmctree.log)"; }
export -f hessian mcmc

echo "== Hessians ($(wc -l < jobs_hessian.txt) jobs, $JOBS in parallel)"
xargs -P "$JOBS" -I{} bash -c 'hessian "$@"' _ {} < jobs_hessian.txt
echo "== MCMC runs ($(wc -l < jobs_mcmc.txt) jobs, $JOBS in parallel)"
xargs -P "$JOBS" -I{} bash -c 'mcmc "$@"' _ {} < jobs_mcmc.txt
echo "== Summary"
python3 summarize_dating.py | tee logs/summary.log
echo "All done $(date '+%H:%M'). Upload: zip -r dating_out.zip results runs/*/FigTree.tre runs/*/out.txt runs/*/mcmctree.log logs calibration_nodes.tsv versions_dating.txt"
