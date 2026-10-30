#!/bin/bash
# old (hla_rsibank_20261026) -> new (hla_hlabank_20261030) per-copy cross-tab: copies matched by >= 50 % reciprocal overlap, same strand
set -euo pipefail
mkdir -p /home/toki/tmp/hla_cmp; cd /home/toki/tmp/hla_cmp || exit 1
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
O=/home/toki/rhin/hla_rsibank_20261026/run_20261026_203020/results
N=/home/toki/rhin/hla_hlabank_20261030/run_20261029_220442/results
tobed() {  # assignment_full.tsv (every locus; Status != assigned -> "unassigned") -> chrom start end label . strand
  awk -F'\t' 'NR>1{print $1"\t"($5=="assigned"?$2:"unassigned")}' "$1/assignment_full.tsv" \
  | awk -F'\t' '{n=$1; st=substr(n,length(n)-1,1); sub(/\([+-]\)$/,"",n); i=match(n,/:[0-9]+-[0-9]+$/);
      c=substr(n,1,i-1); split(substr(n,i+1),a,"-"); print c"\t"a[1]"\t"a[2]"\t"$2"\t.\t"st}' | sort -k1,1 -k2,2n
}
tobed $O > old.bed; tobed $N > new.bed
wc -l old.bed new.bed
bedtools intersect -a old.bed -b new.bed -s -f 0.5 -r -wa -wb | awk -F'\t' '!s[$1":"$2":"$3]++{print $4"\t"$10}' > pairs.tsv
bedtools intersect -a old.bed -b new.bed -s -f 0.5 -r -v | awk -F'\t' '{print $4"\tNO_NEW"}' >> pairs.tsv
bedtools intersect -a new.bed -b old.bed -s -f 0.5 -r -v | awk -F'\t' '{print "NO_OLD\t"$4}' >> pairs.tsv
sort pairs.tsv | uniq -c | sort -k1,1nr | awk '{print $2"\t"$3"\t"$1}' > crosstab.tsv
# per new consensus: copies, and where they came from (top 3 old labels)
awk -F'\t' '{t[$2]+=$3; if(++k[$2]<=3) s[$2]=s[$2] (s[$2]?", ":"") $1" "$3}END{for(n in t)print n"\t"t[n]"\t"s[n]}' crosstab.tsv | sort -t$'\t' -k2,2nr > by_new.tsv
awk -F'\t' '{t[$1]+=$3; if(++k[$1]<=3) s[$1]=s[$1] (s[$1]?", ":"") $2" "$3}END{for(n in t)print n"\t"t[n]"\t"s[n]}' crosstab.tsv | sort -t$'\t' -k2,2nr > by_old.tsv
echo "pairs: $(wc -l < pairs.tsv)"
