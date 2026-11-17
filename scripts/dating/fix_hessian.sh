#!/usr/bin/env bash
# Recover the Hessians already computed: baseml wrote rst2 for codon positions 1+2 and then stopped only at the
# (unneeded) standard-error step. This copies those rst2 files into place and marks the Hessian steps done.
cd "$(dirname "$0")"
for d in runs/hessian_*; do
  r=$(basename "$d")
  if [ -s "$d/tmp0001/rst2" ] && grep -q "Hessian" "$d/tmp0001/rst2" && [ -s "$d/tmp0002.rst2" ]; then
    n1=$(head -2 "$d/tmp0001/rst2" | tr -d ' \n')
    cp "$d/tmp0001/rst2" "$d/tmp0001.rst2"
    cat "$d/tmp0001.rst2" "$d/tmp0002.rst2" > "$d/in.BV"
    cp "$d/tmp0001.rst2" "$d/in_pos12.BV"
    touch "$d/.done"; echo "$r: recovered (taxa line: $n1; in.BV $(wc -l < "$d/in.BV") lines)"
  else
    echo "$r: NOT recoverable - check $d/tmp0001/baseml.log"
  fi
done
