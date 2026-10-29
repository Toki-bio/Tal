#!/bin/bash
# hla bank 2026-10-30: rsi bank21_nor1 adjusted from the peel of the hla step 1 SubFam plate (owner, rounds 1-2) and two groups
# from the colleague's notes on the hla subfamily plates (Tal hla_rsibank_20261026/peel/README.md).
#   replaced: r8_83seqs -> r8_hla; r6_210seqs -> r5_r6_short_hla; r5_r6_P26 -> r5_r6_long_hla; r10_r8_P18 -> consensus_386seqs
#   added: r9_ext, consensus_28seqs, consensus_47seqs, consensus_10seqs_noterm, consensus_12seqs, consensus_10seqs_b,
#          consensus_29seqs (owner); r1_r3_P1_hla020, r5_r5_P48_hla007 (plate rows named in the colleague's notes)
set -euo pipefail
cd /home/toki/rhin/hla_bank_20261030 || exit 1
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
B=/home/toki/rhin/final_bank_20261026
seqkit grep -n -v -p r8_83seqs -p r6_210seqs -p r5_r6_P26 -p r10_r8_P18 $B/bank21_nor1.fa > keep17.fa
cat keep17.fa peel_round1_owner_consensuses.fa peel_round2_owner_consensuses.fa colleague_rows.fa \
  | seqkit seq -u -w 0 | awk '/^>/{print $1; next}{print}' > bank30_hla.fa
awk -F'\t' '$1!="r8_83seqs" && $1!="r6_210seqs" && $1!="r5_r6_P26" && $1!="r10_r8_P18"' $B/rsi_families_nor1.tsv > families30_hla.tsv
grep '^>' peel_round1_owner_consensuses.fa peel_round2_owner_consensuses.fa colleague_rows.fa | sed 's/.*>//; s/ .*//' \
  | awk '{print $1"\tRhin"}' >> families30_hla.tsv
n=$(grep -c '^>' bank30_hla.fa); f=$(wc -l < families30_hla.tsv)
[ "$n" -eq 30 ] && [ "$f" -eq 30 ] || { echo "expected 30 consensuses and 30 family rows, got $n and $f"; exit 1; }
[ -z "$(grep '^>' bank30_hla.fa | sort | uniq -d)" ] || { echo "duplicate names"; exit 1; }
[ -z "$(grep -v '^>' bank30_hla.fa | grep '[^ACGT]')" ] || { echo "non-ACGT bases"; exit 1; }
diff <(grep '^>' bank30_hla.fa | tr -d '>' | sort) <(cut -f1 families30_hla.tsv | sort) && echo "bank30_hla.fa: $n consensuses, families match"
