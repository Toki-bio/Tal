#!/bin/bash
# pack_final.sh <code> <run dir>: publish set (layout as rsi_main_20261022) + family view tables -> <dir>/publish_tal.tgz
set -euo pipefail
code=$1; D=$2
cd $D || exit 1
[ -f DONE ] || { echo "$code not DONE"; exit 1; }
R=$(cat RUN); P=$D/publish_tal
rm -rf "$P"; mkdir -p "$P/flankscan" "$P/alignments" "$P/family_view"
cp $R/results/report.html $R/manifest.txt $R/consensuses.clean.fa "$P/"
for f in loci_verdict.tsv flank_verdict.tsv composite_units.tsv; do [ -f $R/step2/step2_output/$f ] && cp $R/step2/step2_output/$f "$P/"; done
for f in assignment_stats.tsv summary.by_subfam.tsv flank_twins.tsv array_flag.tsv consensus_blocks.tsv; do [ -e $R/results/$f ] && cp -L $R/results/$f "$P/"; done
cp /home/toki/rhin/final_bank_20261026/rsi_families_nor1.tsv "$P/families.tsv"; [ -f $R/families/families.summary.tsv ] && cp $R/families/families.summary.tsv "$P/"
FV=$R/family_view; [ -d $R/family_view_seeds ] && FV=$R/family_view_seeds   # regrouped with SINEderella main 9d484fa (seeds)
for f in family_view.summary.tsv family_view.tsv partial_copies.tsv partial_stretches.tsv; do [ -f $FV/$f ] && cp $FV/$f "$P/family_view/"; done
[ -f $R/genome.clean_step1/searches/composite_units.bank.tsv ] && cp $R/genome.clean_step1/searches/composite_units.bank.tsv "$P/family_view/pieces_from_composite_units.tsv"
for d in satellites length_variants consensus_audit continuation; do [ -e $R/results/$d ] && cp -rL $R/results/$d "$P/"; done
cp $R/results/alignments/* "$P/alignments/"
cp $R/genome.clean_step1/subfam_input/input.clw.al "$P/alignments/${code}_all_copies_subfam30k.aln.fa"   # step 1 SubFam: 600 chunk consensuses of 30,000 copies + the bank
for f in accepted.fa candidates.fa candidates.tsv chains.tsv ends.tsv family_summary.tsv hierarchy.tsv peaks.tsv reassign.tsv reassign_moves.tsv; do [ -f $R/flankscan/$f ] && cp $R/flankscan/$f "$P/flankscan/"; done
[ -d $R/flankscan/cand ] && cp -r $R/flankscan/cand "$P/flankscan/"
grep "^==" $D/pipeline.log > "$P/pipeline_stages.txt"
tar czf publish_tal.tgz -C $D publish_tal && echo "$code $(du -sh $P | cut -f1) $(ls $P/alignments | wc -l) alignments"
