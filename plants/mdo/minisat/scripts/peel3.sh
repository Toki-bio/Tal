set -euo pipefail
export PATH=/staging/tmp/apple_minisat/env/bin:$PATH
cd /staging/tmp/apple_minisat/peel || exit 1
T=${T:-0.70}
rm -rf grp; mkdir grp
# 1. groups: connected components at score >= T, named G01.. by size
awk -F'\t' -v T=$T 'function f(x){while(p[x]!=x)x=p[x]=p[p[x]];return x} FNR==NR{p[$1]=$1; ids[++n]=$1; len[$1]=length($3); next} {if($3>=T){a=f($1);b=f($2);if(a!=b)p[a]=b}} END{for(i=1;i<=n;i++){r=f(ids[i]); c[r]++; mem[r]=mem[r]" "ids[i]} for(r in c) print c[r]"\t"r"\t"mem[r]}' motifs.tsv pair_sym.tsv | sort -k1,1nr -k2,2 | awk -F'\t' '{g=sprintf("G%02d",NR); n=split($3,a," "); for(i=1;i<=n;i++) print a[i]"\t"g"\t"$1}' > groups.tsv
cut -f2,3 groups.tsv | sort -u | sort -k2,2nr | head -30 | tr '\n' ' '; echo
# 2. per group with n>=3: medoid, orient, rotate to medoid phase
for g in $(awk -F'\t' '$3>=3{print $2}' groups.tsv | sort -u); do
  d=grp/$g; mkdir -p $d
  awk -F'\t' -v g=$g '$2==g{print $1}' groups.tsv > $d/ids.txt
  # medoid = member with highest sum of pair scores to other members
  awk -F'\t' 'FNR==NR{in_[$1]=1; next} ($1 in in_) && ($2 in in_){s[$1]+=$3; s[$2]+=$3} END{for(i in in_) print (i in s? s[i]:0), i}' $d/ids.txt pair_sym.tsv | sort -k1,1nr -k2,2 | head -1 | cut -d' ' -f2 > $d/ref.id
  ref=$(cat $d/ref.id)
  awk -F'\t' -v r=$ref '$1==r{print ">"$1"\n"$3}' motifs.tsv > $d/ref.fa
  awk -F'\t' 'FNR==NR{in_[$1]=1; next} ($1 in in_){print ">"$1"\n"$3$3}' $d/ids.txt motifs.tsv > $d/q2x.fa
  makeblastdb -in $d/ref.fa -dbtype nucl -out $d/refdb >/dev/null
  # pass 1: strand of best HSP
  blastn -task blastn -query $d/q2x.fa -db $d/refdb -strand both -dust no -word_size 7 -evalue 1e-3 -max_hsps 5 -outfmt '6 qseqid length sstrand' | sort -k1,1 -k2,2nr | awk '!s[$1]++{print $1"\t"$3}' > $d/strand.tsv
  # oriented single sequences
  awk -F'\t' 'function rc(s,  i,o,c){o="";for(i=length(s);i>=1;i--){c=substr(s,i,1);o=o (c=="A"?"T":c=="C"?"G":c=="G"?"C":c=="T"?"A":"N")}return o} FILENAME==ARGV[1]{st[$1]=$2; next} {if(($1 in st)){s=$3; if(st[$1]=="minus")s=rc(s); print $1"\t"s}}' $d/strand.tsv motifs.tsv > $d/oriented.tsv
  awk -F'\t' '{print ">"$1"\n"$2$2}' $d/oriented.tsv > $d/o2x.fa
  # pass 2: plus-only offset
  blastn -task blastn -query $d/o2x.fa -db $d/refdb -strand plus -dust no -word_size 7 -evalue 1e-3 -max_hsps 5 -outfmt '6 qseqid length qstart sstart' | sort -k1,1 -k2,2nr | awk '!s[$1]++{print $1"\t"$3"\t"$4}' > $d/offset.tsv
  awk -F'\t' 'FILENAME==ARGV[1]{q[$1]=$2; r[$1]=$3; next} {id=$1; s=$2; if(!(id in q)){print id > "'$d'/dropped.txt"; next} L=length(s); off=((q[id]-1)-(r[id]-1))%L; if(off<0)off+=L; print ">"id"\n"substr(s s, off+1, L)}' $d/offset.tsv $d/oriented.tsv > $d/rotated.fa
  n=$(grep -c '>' $d/rotated.fa)
  if [ "$n" -ge 3 ]; then mafft --thread 8 --localpair --maxiterate 1000 $d/rotated.fa > $d/aln.fa 2>$d/mafft.log; fi
done
ls grp | wc -l
