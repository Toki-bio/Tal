#!/bin/bash
# Longer e2-3 consensus, two sets of 100 best bit scores:
#   A  copies with unique flanks (verdict SINE)         B  grey copies of length <= 213 bp with >= 3 continuation words
# Window 50 bp before and 300 bp after the copy. Copies whose DISTANT flank (160-300 bp after the copy end, beyond the
# continuation signal) shares >= 5 twelve-letter words with another copy of the set are removed first. Then one mafft
# L-INS-i alignment; consensus rows by EMBOSS cons at 35, 50 and 60 % of the copies (the end is the owner's call).
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
R=/home/toki/sine_runs/Ata/ten10/run_20261126_160931; G=$R/genome.clean.fa
W=/home/toki/sine_runs/Ata/ten10/e23_long; cd $W; mkdir -p tmp; export TMPDIR=$W/tmp LC_ALL=C
bedtools flank -i short.bed -g $G.fai -l 0 -r 300 -s | bedtools getfasta -fi $G -bed - -s -tab -name 2>/dev/null |
 awk -F'\t' '{s=toupper($2); n=0; m=split("ATGGCCTCCAGGAG GATTCATGGTGCAG CCCAGCAATAACCC AATAATAACCACAA AGGCTACAACAACA ACAAGGGCAACAAA",a," "); for(i=1;i<=m;i++) if(index(s,a[i])) n++; if(n>=3){id=$1; sub(/::.*$/,"",id); print id}}' > cont.ids
awk 'FILENAME==ARGV[1]{w[$1];next} ($4 in w)' cont.ids short.bed | sort -k5,5nr -k4,4 | awk 'NR<=100' > B.bed
echo "short grey copies $(wc -l < short.bed); with >= 3 continuation words $(wc -l < cont.ids); set B $(wc -l < B.bed)"
cp e23.fa e2-3.fa
for X in A B; do
  bedtools flank -i $X.bed -g $G.fai -l 0 -r 300 -s | bedtools getfasta -fi $G -bed - -s -tab -name 2>/dev/null | awk -F'\t' '{id=$1; sub(/::.*$/,"",id); print id"\t"toupper(substr($2,161))}' > $X.far.tsv
  awk -F'\t' '{id[NR]=$1; q[NR]=$2} END{for(i=1;i<=NR;i++){split("",w); for(p=1;p<=length(q[i])-11;p++) w[substr(q[i],p,12)]
     for(j=1;j<=NR;j++) if(j!=i){split("",d); h=0; for(p=1;p<=length(q[j])-11;p++){m=substr(q[j],p,12); if((m in w) && !(m in d)){d[m]; h++}} if(h>=5){print id[i]; break}}}}' $X.far.tsv > $X.shared.ids
  awk 'FILENAME==ARGV[1]{x[$1];next} !($4 in x)' $X.shared.ids $X.bed > $X.u.bed
  echo "set $X: distant flank shared with another copy of the set: $(wc -l < $X.shared.ids); kept $(wc -l < $X.u.bed)"
  bedtools slop -i $X.u.bed -g $G.fai -l 50 -r 300 -s | bedtools getfasta -fi $G -bed - -s -name 2>/dev/null | awk '/^>/{sub(/::.*$/,""); print; next}{print toupper($0)}' > $X.fa
  mafft --thread 16 --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet $X.fa > $X.aln.fa
  n=$(grep -c '^>' $X.aln.fa); : > $X.rows.fa
  for p in 35 50 60; do
    cons -sequence $X.aln.fa -plurality $(awk -v n=$n -v p=$p 'BEGIN{print p*n/100}') -outseq t.cons -auto
    awk -v h=">e2-3_long${X}_${p}pct" '/^>/{next}{s=s $0} END{gsub(/[Nn]/,"-",s); print h; print toupper(s)}' t.cons >> $X.rows.fa
  done
  awk '/^>/{h=$0; next}{t=$0; gsub(/-/,"",t); print h; print t}' $X.rows.fa > $X.cons.fa
  cat $X.rows.fa $X.aln.fa > $X.stack.fa
  mafft --thread 16 --nuc --preservecase --quiet --add e2-3.fa $X.stack.fa > $X.add.fa
  # order: 3 consensus rows, original e2-3, copies; outside the 35 % consensus span: lower case, degapped, packed
  awk '/^>/{k++; h[k]=$0; next}{s[k]=s[k] $0}
    END{L=length(s[1]); for(j=1;j<=L;j++) if(substr(s[1],j,1)!="-"){if(!a)a=j; z=j}
      m=0; for(i=1;i<=3;i++) ord[++m]=i; ord[++m]=k; for(i=4;i<k;i++) ord[++m]=i
      for(r=1;r<=m;r++){i=ord[r]; l=tolower(substr(s[i],1,a-1)); gsub(/-/,"",l); t=tolower(substr(s[i],z+1)); gsub(/-/,"",t)
        pl=""; for(j=length(l);j<a-1;j++) pl=pl "-"; pr=""; for(j=length(t);j<L-z;j++) pr=pr "-"
        print h[i]; print pl l toupper(substr(s[i],a,z-a+1)) t pr}}' $X.add.fa > aal_e2-3_long$X.aln.fa
  echo "set $X: plate rows $(grep -c '^>' aal_e2-3_long$X.aln.fa); consensus lengths: $(awk '/^>/{h=substr($0,2);next}{printf "%s=%d ",h,length($0)}' $X.cons.fa)"
  awk -F'\t' -v X=$X '{n[int(($3-$2)/20)*20]++} END{printf "set %s copy lengths:", X; for(k in n) printf " %d-%d:%d", k,k+19,n[k]; print ""}' $X.u.bed
done
cat e2-3.fa A.cons.fa B.cons.fa > e23_versions.fa
mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --preservecase --quiet e23_versions.fa > aal_e2-3_versions.aln.fa
awk '/^>/{print; next}{print}' e23_versions.fa | awk '/^>/{printf "%s ",$0; next}{print substr($0,1,20) " ... " substr($0,length($0)-19) }'
rm -rf tmp t.cons A.in.tmp A.plate.raw bits.tsv v.tsv
md5sum aal_e2-3_*.aln.fa
echo DONE
