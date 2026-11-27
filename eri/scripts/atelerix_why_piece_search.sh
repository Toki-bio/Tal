#!/bin/bash
# New unassigned places of ten11 (no ten10 locus): which search's own hit file produced them?
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
A=/home/toki/sine_runs/Ata/ten10/run_20261126_160931; B=/home/toki/sine_runs/Ata/ten11/run_20261126_214359
W=/home/toki/tmp/ata_why2; mkdir -p $W; cd $W; export TMPDIR=$W LC_ALL=C
loci(){ O=$(ls -dt $1/step2/step2_output* | head -1); awk -F'\t' -v OFS='\t' 'NR>1{x=$1; sub(/\([^()]*\)$/,"",x); i=index(x,":"); split(substr(x,i+1),a,"-"); print substr(x,1,i-1),a[1],a[2],$2,$3,a[2]-a[1]}' $O/loci_verdict.tsv | sort -k1,1 -k2,2n; }
loci $A > a.bed; loci $B > b.bed
bedtools intersect -a b.bed -b a.bed -v | awk '$5!="assigned"' > nu.bed; N=$(wc -l < nu.bed); echo "new unassigned places: $N, mean length $(awk '{s+=$6} END{printf "%.0f", s/NR}' nu.bed) bp"
S=$B/genome.clean_step1/searches; ls $S/*.bed | awk '!/before|labeled|regions/{n=$0; sub(/.*\//,"",n); print n}'
: > ov.tsv
for f in $S/gen-*.bed; do q=$(basename $f .bed); q=${q#gen-}; q=${q#*.part_}; cut -f1-3 $f | sort -k1,1 -k2,2n > q.bed; bedtools intersect -sorted -a nu.bed -b q.bed -u | awk -v q=$q '{print $1":"$2":"$3"\t"q}' >> ov.tsv; echo "$q hits $(wc -l < q.bed), covers new unassigned places: $(awk -v q=$q '$2==q' ov.tsv | wc -l)"; done
awk -F'\t' -v N=$N '{if($2~/piece/) p[$1]=(p[$1]" "$2); else f[$1]=1; all[$1]=1} END{for(k in all){ if(!(k in f)){onlyp++; if(p[k]~/281-349/) t++} else full++} printf "found only by a piece search: %d (%.1f%%), of them by the 281-349 piece: %d\nfound by a whole-consensus search too: %d (%.1f%%)\n", onlyp, 100*onlyp/N, t, full, 100*full/N}' ov.tsv
echo "== ten10 for comparison: places of ten10 unassigned"; awk '$5!="assigned"{n++; s+=$6} END{printf "%d, mean length %.0f bp\n", n, s/n}' a.bed
cd /home/toki/tmp && rm -rf $W && echo "removed $W"
