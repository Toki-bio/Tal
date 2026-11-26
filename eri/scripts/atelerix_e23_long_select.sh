#!/bin/bash
# Longer e2-3 consensus from the genome copies. Two sets, 100 best bit scores each:
#   A  copies with unique flanks (verdict SINE), end already moved by step 1
#   B  grey copies that run on past their end (length <= 213 bp and >= 3 of the continuation words in the next 300 bp)
# Window: 50 bp before, 300 bp after. mafft L-INS-i; EMBOSS cons with plurality 35 % of the copies (uncalled = gap), as the
# owner's `conse`. Then copies whose 100 bp after the new end are shared with another copy of the set (>= 5 common 12-mers)
# are removed (for set B the genome-wide repeat mask result is also reported), and alignment and consensus are built once more.
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
R=/home/toki/sine_runs/Ata/ten10/run_20261126_160931; G=$R/genome.clean.fa; S2=$R/step2/step2_output
S=/home/toki/SINEderella_main_20261126/flankscan; K=$R/genome.clean.fa.k20gt20.tsv
W=/home/toki/sine_runs/Ata/ten10/e23_long; mkdir -p $W/tmp; cd $W; export TMPDIR=$W/tmp LC_ALL=C
awk -F'\t' '$2=="e2-3"{print $1"\t"$3}' $S2/assignment_full.tsv > bits.tsv
awk -F'\t' 'NR==FNR{b[$1]=$2; next} $2=="e2-3" && ($4=="SINE"||$4=="grey") && ($1 in b){print $1"\t"$4"\t"b[$1]}' bits.tsv $S2/loci_verdict.tsv > v.tsv
awk -F'\t' -v OFS='\t' '{x=$1; st=substr(x,length(x)-1,1); sub(/\([+-]\)$/,"",x); i=index(x,":"); split(substr(x,i+1),a,"-"); print substr(x,1,i-1),a[1],a[2],$1,$3,st > ($2".bed")}' v.tsv
sort -k5,5nr -k4,4 SINE.bed | awk 'NR<=100' > A.bed
awk '$3-$2<=213' grey.bed > short.bed
bedtools flank -i short.bed -g $G.fai -l 0 -r 300 -s | bedtools getfasta -fi $G -bed - -s -tab -name 2>/dev/null |
 awk -F'\t' '{s=toupper($2); n=0; w="ATGGCCTCCAGGAG GATTCATGGTGCAG CCCAGCAATAACCC AATAATAACCACAA AGGCTACAACAACA ACAAGGGCAACAAA"; m=split(w,a," "); for(i=1;i<=m;i++) if(index(s,a[i])) n++; if(n>=3){id=$1; sub(/\([+-]\)$/,"",id); print id}}' > cont.ids
awk 'NR==FNR{w[$1];next} ($4 in w)' cont.ids short.bed | sort -k5,5nr -k4,4 | awk 'NR<=100' > B.bed
echo "short grey copies $(wc -l < short.bed), of them with >= 3 continuation words $(wc -l < cont.ids)"
awk '/^>/{p=($1==">e2-3")} p' $R/consensuses.clean.fa > e23.fa
build(){   # name bed
  bedtools slop -i $2 -g $G.fai -l 50 -r 300 -s | bedtools getfasta -fi $G -bed - -s -name 2>/dev/null | awk '/^>/{sub(/\([+-]\)$/,""); print; next}{print toupper($0)}' > $1.fa
  mafft --thread 16 --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet $1.fa > $1.aln.fa
  n=$(grep -c '^>' $1.aln.fa); p=$(awk -v n=$n 'BEGIN{printf "%.2f", 0.35*n}')
  cons -sequence $1.aln.fa -plurality $p -name e2-3_long$1 -outseq $1.cons.raw -auto
  awk -v n=e2-3_long$1 '/^>/{next}{s=s $0} END{gsub(/[Nn-]/,"",s); print ">" n; print toupper(s)}' $1.cons.raw > $1.cons.fa
  echo "set $1: copies $n, consensus $(awk 'NR==2{print length($0)}' $1.cons.fa) bp"
}
for X in A B; do
  build $X $X.bed
  # 100 bp after the new consensus end, per copy: the column of the last consensus base, then the next 100 bases of the copy
  awk -v n=$(grep -c '^>' $X.aln.fa) '/^>/{k++; id[k]=substr($0,2); next}{s[k]=s[k] toupper($0)}
    END{L=length(s[1]); last=0; for(j=1;j<=L;j++){split("",c); for(i=1;i<=k;i++){b=substr(s[i],j,1); if(b!="-")c[b]++} for(b in c) if(c[b]>=0.35*n && j>last) last=j}
      for(i=1;i<=k;i++){t=substr(s[i],last+1); gsub(/-/,"",t); print id[i]"\t3\t"substr(t,1,100)}}' $X.aln.fa > $X.after.tsv
  awk -F'\t' '{id[NR]=$1; q[NR]=$3} END{for(i=1;i<=NR;i++){split("",w); for(p=1;p<=length(q[i])-11;p++) w[substr(q[i],p,12)]
      for(j=1;j<=NR;j++) if(j!=i){split("",d); h=0; for(p=1;p<=length(q[j])-11;p++){m=substr(q[j],p,12); if((m in w) && !(m in d)){d[m]; h++}} if(h>=5){print id[i]; break}}}}' $X.after.tsv > $X.shared.ids
  python3 -I $S/fs9_repmask.py rep $X.after.tsv $K 20 > $X.rep
  masked=$(python3 -I $S/fs9_repmask.py mask $X.after.tsv $X.rep | awk -F'\t' '{s=$3; gsub(/N/,"",s); if(length(s)<50) n++} END{print n+0}')
  echo "set $X: 100 bp after the end shared with another copy of the set: $(wc -l < $X.shared.ids); mostly repeated genome-wide (< 50 unmasked bp): $masked; short flank (< 50 bp left in window): $(awk -F'\t' 'length($3)<50' $X.after.tsv | wc -l)"
  awk 'NR==FNR{x[$1];next} !($4 in x)' $X.shared.ids $X.bed > $X.u.bed
  cp $X.cons.fa $X.cons.first.fa; cp $X.aln.fa $X.first.aln.fa
  [ $(wc -l < $X.u.bed) -ge 10 ] && build $X $X.u.bed || echo "set $X: fewer than 10 copies with unshared flanks; consensus from all copies kept"
  # plate: new consensus, original e2-3, copies; outside the consensus span flanks are lower case, degapped and packed
  cat $X.cons.fa e23.fa $X.aln.fa | awk '/^>/{print; next}{print}' > $X.in.tmp
  mafft --thread 16 --localpair --maxiterate 1000 --ep 0.123 --nuc --preservecase --quiet --add <(cat $X.cons.fa e23.fa) $X.aln.fa > $X.plate.raw
  awk '/^>/{k++; h[k]=$0; next}{s[k]=s[k] $0}
    END{for(i=1;i<=k;i++) if(h[i]~/e2-3_long/) c=i; for(i=1;i<=k;i++) if(h[i]==">e2-3") o=i
      L=length(s[c]); for(j=1;j<=L;j++) if(substr(s[c],j,1)!="-"){if(!a)a=j; z=j}
      order[++m]=c; order[++m]=o; for(i=1;i<=k;i++) if(i!=c && i!=o) order[++m]=i
      for(r=1;r<=m;r++){i=order[r]; l=tolower(substr(s[i],1,a-1)); gsub(/-/,"",l); t=tolower(substr(s[i],z+1)); gsub(/-/,"",t)
        pl=""; for(j=length(l);j<a-1;j++) pl=pl "-"; pr=""; for(j=length(t);j<L-z;j++) pr=pr "-"
        print h[i]; print pl l toupper(substr(s[i],a,z-a+1)) t pr}}' $X.plate.raw > aal_e2-3_long$X.aln.fa
  echo "set $X plate: $(grep -c '^>' aal_e2-3_long$X.aln.fa) rows"
done
echo "== consensuses"; cat A.cons.fa B.cons.fa e23.fa | tee e23_versions.fa
echo "== pairwise to the original and to each other (ssearch36)"
for X in A B; do cat $X.cons.fa; done > q.fa; cat e23.fa B.cons.fa > l.fa
ssearch36 -n -m 8CB -E 1 -3 q.fa l.fa 2>/dev/null | awk '!/^#/ && $1!=$2{print $1, $2, "identity "$3, "aligned "$4, "q "$7"-"$8, "s "$9"-"$10}'
mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --preservecase --quiet e23_versions.fa > aal_e2-3_versions.aln.fa
rm -rf tmp *.tmp
echo DONE
