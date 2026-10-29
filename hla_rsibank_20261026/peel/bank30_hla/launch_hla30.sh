#!/bin/bash
# hla run 2026-10-30: bank30_hla (rsi bank adjusted from the hla peel, see build_bank.sh), SINEderella branch end-extension (bf2480a)
# families30_hla.tsv. Usage: launch_hla30.sh CODE GENOME
set -uo pipefail
SP=$1; G=$2
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella; export TMPDIR=$HOME/tmp
C=/home/toki/SINEderella_endext; D=$(pwd -P); cd "$D" || exit 1; L=$D/pipeline.log
B=/home/toki/rhin/hla_bank_20261030/bank30_hla.fa; FT=/home/toki/rhin/hla_bank_20261030/families30_hla.tsv
rm -f DONE FAILED
export TRNA_REF=/home/toki/refs/smallrna/hg19-tRNAs.fa RFAM_CM=/home/toki/refs/smallrna/Rfam.cm RNATOOLS_BIN=/home/toki/miniforge3/envs/rnatools/bin
RAW=https://raw.githubusercontent.com/Toki-bio/Tal/main/$(basename "$D")/alignments/
{
echo "== start $(date '+%F %T') $SP code $(git -C $C log -1 --format=%h) bank $B genome $G"
THREADS=32 SPECIES_CODE=$SP PUBLISH=1 RAW_ALN_BASE=$RAW /usr/bin/time -v -o $D/run.time $C/SINEderella --families $FT "$G" $B
rc=$?; echo "== SINEderella rc=$rc $(date '+%F %T')"
R=$(ls -dt $D/run_* 2>/dev/null | head -1); echo "$R" > $D/RUN
[ $rc -eq 0 ] && [ -s "$R/results/report.html" ] || { touch $D/FAILED; exit 1; }
echo "== flankscan 1-7 $(date '+%F %T')"
/usr/bin/time -v -o $D/fs.time bash $C/flankscan/fs_all.sh $R $R/flankscan 32 > $D/fs.log 2>&1 || { echo "fs_all failed"; touch $D/FAILED; exit 1; }
echo "== report rebuild through publish_run $(date '+%F %T')"
cd $R || exit 1
SKIP_ALIGN=1 SKIP_STEP4=1 THREADS=32 SINEDERELLA_BIN=$C RAW_ALN_BASE=$RAW bash $C/publish_run.sh "$(pwd -P)" $SP > $D/report_rebuild.log 2>&1 < /dev/null || { echo "report rebuild failed"; touch $D/FAILED; exit 1; }
echo "== done $(date '+%F %T')"
touch $D/DONE
} >> $L 2>&1
