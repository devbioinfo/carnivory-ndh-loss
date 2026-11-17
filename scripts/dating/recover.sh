#!/usr/bin/env bash
# Resume the long dating chains after any interruption (power cut, restart) WITHOUT losing samples already collected.
# Safe to run again after every interruption. Do NOT use rerun_long.sh any more (it deletes collected samples).
#   cd ~/Downloads/mac_dating_v1
#   JOBS=12 nohup nice -n 10 caffeinate -i bash recover.sh > recover.txt 2>&1 & disown
# For every unfinished chain:
#   - samples already written (mcmc.txt) are kept as mcmc_partN.txt (a cut-off last line is dropped);
#   - the chain continues for the samples still missing (target 10,000 per chain), with checkpointing switched on
#     (state saved every 5% to mcmctree.ckpt; after a later interruption the chain resumes from that state);
#   - when it finishes, all parts are joined (mcmc_all.txt) and summarised by MCMCTree (print = -1) into FigTree.tre.
cd "$(dirname "$0")"
JOBS=${JOBS:-12}; TARGET=${TARGET:-10000}
[ -x paml_src/src/mcmctree ] && export PATH="$PWD/paml_src/src:$PATH"
MCMCTREE=$(command -v mcmctree) || { echo "mcmctree not found"; exit 1; }; export MCMCTREE TARGET

prepare () { local d=$1; [ -d "$d" ] || return 0
  [ -e "$d/.done" ] && return 0
  if [ -f "$d/mcmc.txt" ] && [ "$(wc -l < "$d/mcmc.txt")" -gt 1 ]; then
    local n=1; while [ -e "$d/mcmc_part$n.txt" ]; do n=$((n+1)); done
    awk -F'\t' 'NR==1{nf=NF; print; next} NF==nf' "$d/mcmc.txt" > "$d/mcmc_part$n.txt" && rm -f "$d/mcmc.txt"
  fi
  local have=0 f
  for f in "$d"/mcmc_part*.txt; do [ -e "$f" ] && have=$((have + $(wc -l < "$f") - 1)); done
  local need=$((TARGET - have)); [ $need -lt 0 ] && need=0
  perl -pi -e "s/^\s*nsample\s*=.*/       nsample = $need/" "$d/mcmctree.ctl"
  perl -ni -e 'print unless /^\s*checkpoint\s*=/' "$d/mcmctree.ctl"
  if [ -s "$d/mcmctree.ckpt" ]; then echo "    checkpoint = 2 0.05 mcmctree.ckpt" >> "$d/mcmctree.ctl"
  else
    echo "    checkpoint = 1 0.05 mcmctree.ckpt" >> "$d/mcmctree.ctl"
    # fresh chain: a different seed from the earlier part, so the parts are not copies of each other
    local s; s=$(awk '/^[ \t]*seed[ \t]*=/{print $3}' "$d/mcmctree.ctl"); n=0; for f in "$d"/mcmc_part*.txt; do [ -e "$f" ] && n=$((n+1)); done
    perl -pi -e "s/^\s*seed\s*=.*/          seed = $((s + 1000 * n))/" "$d/mcmctree.ctl"
  fi
  echo "$(basename "$d"): $have samples kept, $need to go"
}

combine () { local d=$1
  ( cd "$d" || exit 1
    ls mcmc_part*.txt >/dev/null 2>&1 || exit 0
    { head -1 mcmc_part1.txt; for f in mcmc_part*.txt mcmc.txt; do [ -s "$f" ] && tail -n +2 "$f"; done; } > mcmc_all.txt
    perl -pe 's/^\s*print\s*=.*/         print = -1/; s/^\s*mcmcfile\s*=.*/      mcmcfile = mcmc_all.txt/; s/^\s*outfile\s*=.*/       outfile = out_combined.txt/; s/^\s*checkpoint\s*=.*//' mcmctree.ctl > summary.ctl
    "$MCMCTREE" summary.ctl > summary.log 2>&1
    cp mcmc_all.txt mcmc.txt )
}

runone () { local spec=$1; local run=${spec%%|*}; local d="runs/$run"
  [ -e "$d/.done" ] && return 0
  local need; need=$(awk '/^[ \t]*nsample[ \t]*=/{print $3}' "$d/mcmctree.ctl")
  echo "[$(date '+%H:%M')] start $run ($need samples)"
  if [ "$need" -gt 0 ]; then ( cd "$d" && "$MCMCTREE" mcmctree.ctl > mcmctree.log 2>&1 ) || { echo "[$(date '+%H:%M')] FAILED $run"; return 1; }; fi
  combine "$d"
  [ -s "$d/FigTree.tre" ] && touch "$d/.done" && echo "[$(date '+%H:%M')] done  $run" || echo "[$(date '+%H:%M')] FAILED $run (no FigTree.tre)"
}
export -f runone combine

for d in runs/post_* runs/sens_*; do prepare "$d"; done
grep -v "^prior_" jobs_mcmc.txt | xargs -P "$JOBS" -I{} bash -c 'runone "$@"' _ {}
echo "== Summary"; python3 summarize_dating.py | tee logs/summary_long.log
echo "All done $(date '+%H:%M')."
