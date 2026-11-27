#!/bin/bash
# A: why did unassigned loci go from 227,559 (ten10) to 472,641 (ten11)?  B: where did the UnS1 / UnS2 hits go?
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
A=/home/toki/sine_runs/Ata/ten10/run_20261126_160931; B=/home/toki/sine_runs/Ata/ten11/run_20261126_214359; U=/home/toki/sine_runs/Ata/uns4/run_add_20261127_005624
W=/home/toki/tmp/ata_why; mkdir -p $W; cd $W; export TMPDIR=$W LC_ALL=C
loci(){ # run -> chrom start end subfam status length
  O=$(ls -dt $1/step2/step2_output* | head -1)
  awk -F'\t' -v OFS='\t' 'NR>1{x=$1; sub(/\([^()]*\)$/,"",x); i=index(x,":"); split(substr(x,i+1),a,"-"); st=($3=="assigned"?"assigned":"unassigned"); print substr(x,1,i-1),a[1],a[2],(st=="assigned"?$2:"-"),st,a[2]-a[1]}' $O/loci_verdict.tsv | sort -k1,1 -k2,2n; }
loci $A > a.bed; loci $B > b.bed
awk '{n[$5]++} END{print "ten10 loci", NR, "unassigned", n["unassigned"]}' a.bed; awk '{n[$5]++} END{print "ten11 loci", NR, "unassigned", n["unassigned"]}' b.bed
awk '$5=="unassigned"' b.bed > bu.bed
echo "== ten11 unassigned loci: what covered the same place in ten10 (at least half of the ten11 locus)"
bedtools intersect -a bu.bed -b a.bed -f 0.5 -wao | awk -F'\t' '{k=$1":"$2":"$3; if($NF>0){ if($11=="assigned") x[k]="a"; else if(x[k]!="a") x[k]="u"} else if(!(k in x)) x[k]="n"} END{for(k in x) c[x[k]]++; printf "assigned in ten10 %d, unassigned in ten10 %d, no ten10 locus %d\n", c["a"], c["u"], c["n"]}'
echo "== ten11 unassigned loci: which queries hit them (query = consensus or piece)"
awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' $B/genome.clean_step1/searches/all_hits.labeled.bed | sort -k1,1 -k2,2n > hits.bed
awk 'NR<=3' hits.bed
bedtools intersect -sorted -a bu.bed -b hits.bed -s -wa -wb 2>/dev/null | head -1 >/dev/null || true
bedtools intersect -sorted -a bu.bed -b hits.bed -wa -wb | awk -F'\t' '{k=$1":"$2":"$3; q=$10; len[k]=$6; if(index(","s[k]",", ","q",")==0) s[k]=(s[k]==""?q:s[k]","q)} END{for(k in s){n=split(s[k],a,","); onlyp=1; p280=0; p1=0; other=0
   for(i=1;i<=n;i++){ if(a[i]!~/piece/) onlyp=0; if(a[i]~/piece_281-349/) p280=1; else if(a[i]~/piece/) p1=1; else other=1}
   cls=(p280 && !p1 && !other)?"only piece 281-349 of long B":(onlyp?"only other pieces":(p280?"piece 281-349 plus others":"no tail piece")); C[cls]++; L[cls]+=len[k]; T++}
 for(c in C) printf "%s\t%d\t%.1f%%\tmean length %.0f bp\n", c, C[c], 100*C[c]/T, L[c]/C[c]}'
echo "== the same split for ten11 loci with no ten10 locus"
bedtools intersect -a b.bed -b a.bed -v > new.bed; awk '{n[$5]++} END{printf "new places %d: assigned %d, unassigned %d\n", NR, n["assigned"], n["unassigned"]}' new.bed
bedtools intersect -sorted -a new.bed -b hits.bed -wa -wb | awk -F'\t' '{k=$1":"$2":"$3; if($10=="e2-3longB__piece_281-349"||$10~/piece_281-349/) t[k]=1; else o[k]=1} END{for(k in t) if(!(k in o)) a++; for(k in t) if(k in o) b++; printf "new places hit only by the tail piece %d, tail piece and others %d\n", a, b}'
echo "== piece labels present"; awk '$4~/piece/{n[$4]++} END{for(k in n) print k, n[k]}' hits.bed
echo "== B: fate of the UnS1 and UnS2 hits in the --add run (locus covering at least half of the hit)"
loci $U > u.bed
for n in UnS1 UnS2 UnS3 UnS4; do f=$(ls $U/genome.clean_step1/searches/gen-*part_$n.bed | head -1); cut -f1-3 $f | sort -k1,1 -k2,2n > q.bed
  echo "-- $n: $(wc -l < q.bed) hits"
  bedtools intersect -a q.bed -b u.bed -f 0.5 -wao | awk -F'\t' '{k=$1":"$2":"$3; if(k in d) next; d[k]=1; if($NF==0) c["no locus"]++; else if($8=="assigned") c["assigned to " $7]++; else c["unassigned"]++; T++} END{for(x in c) printf "%d\t%.1f%%\t%s\n", c[x], 100*c[x]/T, x}' | sort -k1,1nr | head -8
  bedtools intersect -a q.bed -b u.bed -f 0.5 -wa -wb | awk -F'\t' '{h+=$3-$2; l+=$9; n++} END{if(n) printf "mean hit %.0f bp, mean covering locus %.0f bp\n", h/n, l/n}'
done
cd /home/toki/tmp && rm -rf $W && echo "removed $W"
