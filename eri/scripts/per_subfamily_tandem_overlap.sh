#!/bin/bash
set -euo pipefail
cd /data/W/toki/Genomes/Mammalia/Eulipotyphla/Erniacidae/run_20260820_221537

echo "=== per-subfamily tandem overlap ==="
for sf in e1-1 e1-2 e1-3consensus_32seqs e1-4 e2-1 e2-2 e2-3 e2-4; do
  total=$(awk -F'\t' -v s="$sf" 'index($4, s"|")==1' sine_hits.sorted.bed | wc -l)
  tand=$(awk -F'\t' -v s="$sf" 'index($4, s"|")==1' sine_hits.in_tandem.bed | wc -l)
  pct=$(awk -v a="$tand" -v b="$total" 'BEGIN{if(b>0) printf "%.1f", 100*a/b; else print 0}')
  echo "$sf: total=$total tandem=$tand pct=${pct}%"
done
