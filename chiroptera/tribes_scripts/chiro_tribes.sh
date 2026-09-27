#!/usr/bin/env bash
# Tribes (99% vsearch clusters, top 10 aligned with 50L/70R flanks) for the bat runs.
# Usage: chiro_tribes.sh CODE:RUN_ROOT ...   (PAR genomes at once, THREADS each)
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
export TMPDIR=$HOME/tmp THREADS=${THREADS:-12}
PAR=${PAR:-4}
LOG=~/chiro/tribes/tribes.log; mkdir -p ~/chiro/tribes
one(){
  c=${1%%:*}; r=${1#*:}; out=~/chiro/tribes/$c
  if ls $out/${c}_tribe*_*seqs.aln.fa >/dev/null 2>&1 && [[ -s $out/tribes_metadata.txt ]]; then echo "[$(date '+%F %T')] SKIP $c (done)" >> $LOG; return; fi
  echo "[$(date '+%F %T')] START $c $r" >> $LOG
  bash ~/chiro/extract_tribes.sh "$r" "$c" "$out" > $out.log 2>&1; rc=$?
  n=$(( $(wc -l < $out/${c}_tribes_summary.tsv 2>/dev/null || echo 1) - 1 ))
  echo "[$(date '+%F %T')] END $c exit=$rc top_tribes=$n" >> $LOG
}
export -f one; export LOG
printf '%s\n' "$@" | xargs -P "$PAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date '+%F %T')] ALL_DONE" >> $LOG
