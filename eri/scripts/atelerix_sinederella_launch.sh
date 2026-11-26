#!/bin/bash
# SINEderella on Atelerix albiventris GWHEQWD00000000 with the ten Tal eri consensuses
# (ERI-1, ERI-2 literature + e1-1..e2-4 E. europaeus candidates). PMC12026660 check, 2026-11-26.
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella; export TMPDIR=$HOME/tmp
D=/home/toki/sine_runs/Ata/ten10; cd $D || exit 1; L=$D/pipeline.log
rm -f $D/DONE $D/FAILED
{
echo "== start $(date '+%F %T') code $(git -C /home/toki/SINEderella_main_20261126 log -1 --format=%h) bank $D/bank10.fa"
THREADS=64 SPECIES_CODE=aal /usr/bin/time -v -o $D/run.time \
  /home/toki/SINEderella_main_20261126/SINEderella /home/toki/tmp/atelerix/GWHEQWD00000000.genome.fasta $D/bank10.fa
rc=$?
echo "== SINEderella rc=$rc $(date '+%F %T')"
R=$(ls -dt $D/run_* 2>/dev/null | head -1); echo "$R" > $D/RUN
if [ $rc -eq 0 ] && [ -s $R/results/report.html ]; then touch $D/DONE; else touch $D/FAILED; fi
} >> $L 2>&1
