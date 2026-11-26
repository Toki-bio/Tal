#!/bin/bash
# Do long A / long B (called on Atelerix copies) exist and end in the same place in E. europaeus? Copies are picked by the
# OLD 203 bp e2-3 only (the head), so the choice does not force any particular continuation. Then they are sorted by what
# follows, and each type gets the same alignment and consensus as on Atelerix.
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
SD=/home/toki/SINEderella_main_20261127; A=/home/toki/sine_runs/Ata/ten10/e23_long
G=/home/toki/tmp/gymnure_keratin/genomes/EEU_pri.fa; W=/home/toki/sine_runs/Eeu/e23_check
mkdir -p $W/tmp $W/search; cd $W/search; export TMPDIR=$W/tmp LC_ALL=C
cp $A/e2-3.fa e23.fa; ln -sf $G genome.fa
if ! ls *-e23.bed >/dev/null 2>&1; then THREADS=32 taskset -c 0-31 bash $SD/sear -k e23.fa genome.fa 0.8 65 50 > sear.log 2>&1; fi
B=$(ls *-e23.bed); cd $W
echo "e2-3 head hits: $(wc -l < search/$B)"; head -2 search/$B
[ -s $G.fai ] || samtools faidx $G
awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3"("$6")",$5,$6}' search/$B | sort -k5,5nr -k4,4 > hits.bed
bedtools flank -i hits.bed -g $G.fai -l 0 -r 300 -s | bedtools getfasta -fi $G -bed - -s -tab -name 2>/dev/null |
 awk -F'\t' -v OFS='\t' '{id=$1; sub(/::.*$/,"",id); sub(/\([+-]\)\([+-]\)$/, substr(id,length(id)-2), id); s=toupper($2)
   nb=(index(s,"ATGGCCTCCAGGAG")>0)+(index(s,"GATTCATGGTGCAG")>0)+(index(s,"CCCAGCAATAACCC")>0)
   a=(index(s,"AACAAGGGCAACAA")>0)+(index(s,"ACAACAATAATAAC")>0)+(index(s,"AATAATAACCACAA")>0)
   print id, (nb>=2?"B":(nb==0&&a>=1?"A":(nb==1?"oneB":"neither")))}' > type.tsv
echo "== what follows the e2-3 head in E. europaeus (all hits)"; awk -F'\t' '{n[$2]++} END{for(k in n) printf "%s\t%d\t%.1f%%\n",k,n[k],100*n[k]/NR}' type.tsv | sort
for X in A B; do
  awk -F'\t' -v X=$X 'FILENAME==ARGV[1]{if($2==X)w[$1];next} ($4 in w)' type.tsv hits.bed | awk 'NR<=200' > $X.all.bed
  echo "set $X: best-scoring candidates $(wc -l < $X.all.bed)"
  [ $(wc -l < $X.all.bed) -ge 10 ] || continue
  bedtools flank -i $X.all.bed -g $G.fai -l 0 -r 300 -s | bedtools getfasta -fi $G -bed - -s -tab -name 2>/dev/null | awk -F'\t' '{id=$1; sub(/::.*$/,"",id); print id"\t"toupper(substr($2,161))}' > $X.far.tsv
  awk -F'\t' '{id[NR]=$1; q[NR]=$2} END{for(i=1;i<=NR;i++){split("",w); for(p=1;p<=length(q[i])-11;p++) w[substr(q[i],p,12)]
     for(j=1;j<=NR;j++) if(j!=i){split("",d); h=0; for(p=1;p<=length(q[j])-11;p++){m=substr(q[j],p,12); if((m in w) && !(m in d)){d[m]; h++}} if(h>=5){print id[i]; break}}}}' $X.far.tsv > $X.shared.ids
  awk -F'\t' 'FILENAME==ARGV[1]{x[$1];next} !($4 in x)' <(cat $X.shared.ids; echo -) $X.far.tsv >/dev/null
  # ids in far.tsv carry the strand twice after getfasta; compare by coordinates
  awk 'FILENAME==ARGV[1]{k=$0; sub(/\(.*$/,"",k); x[k];next} {k=$4; sub(/\(.*$/,"",k); if(!(k in x)) print}' $X.shared.ids $X.all.bed | awk 'NR<=100' > $X.bed
  echo "set $X: candidates with a shared distant flank $(wc -l < $X.shared.ids) of $(wc -l < $X.all.bed); aligned $(wc -l < $X.bed)"
  bedtools slop -i $X.bed -g $G.fai -l 50 -r 300 -s | bedtools getfasta -fi $G -bed - -s -name 2>/dev/null | awk '/^>/{sub(/::.*$/,""); print; next}{print toupper($0)}' > $X.fa
  mafft --thread 16 --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet $X.fa > $X.aln.fa
  n=$(grep -c '^>' $X.aln.fa); : > $X.rows.fa
  for p in 35 50 60; do cons -sequence $X.aln.fa -plurality $(awk -v n=$n -v p=$p 'BEGIN{print p*n/100}') -outseq t.cons -auto
    awk -v h=">Eeu_long${X}_${p}pct" '/^>/{next}{s=s $0} END{gsub(/[Nn]/,"-",s); print h; print toupper(s)}' t.cons >> $X.rows.fa; done
  awk -v X=$X '$0==">e2-3longA"||$0==">e2-3longB"{p=(substr($0,length($0))==X); if(p){print ">Ata_long" X; next}} /^>/{p=0} p' /home/toki/sine_runs/Ata/ten11/bank11.fa > $X.ata.fa
  cat $X.rows.fa $X.aln.fa > $X.stack.fa
  mafft --thread 16 --nuc --preservecase --quiet --add <(cat $X.ata.fa $A/e2-3.fa) $X.stack.fa > $X.add.fa
  awk '/^>/{k++; h[k]=$0; next}{s[k]=s[k] $0}
    END{L=length(s[1]); for(j=1;j<=L;j++) if(substr(s[1],j,1)!="-"){if(!a)a=j; z=j}
      m=0; for(i=1;i<=3;i++) ord[++m]=i; ord[++m]=k-1; ord[++m]=k; for(i=4;i<k-1;i++) ord[++m]=i
      for(r=1;r<=m;r++){i=ord[r]; l=tolower(substr(s[i],1,a-1)); gsub(/-/,"",l); t=tolower(substr(s[i],z+1)); gsub(/-/,"",t)
        pl=""; for(j=length(l);j<a-1;j++) pl=pl "-"; pr=""; for(j=length(t);j<L-z;j++) pr=pr "-"
        print h[i]; print pl l toupper(substr(s[i],a,z-a+1)) t pr}}' $X.add.fa > eeu_e2-3_long$X.aln.fa
  echo "set $X: plate rows $(grep -c '^>' eeu_e2-3_long$X.aln.fa); consensus lengths: $(awk '/^>/{h=substr($0,2);next}{t=$0; gsub(/-/,"",t); printf "%s=%d ",h,length(t)}' $X.rows.fa)"
  echo "set $X 60% consensus:"; awk 'NR==6{gsub(/-/,""); print}' $X.rows.fa
  awk -v X=$X '/^>/{k++; h[k]=substr($0,2); next}{s[k]=s[k] $0}
    END{for(i=1;i<=k;i++){if(h[i]~/_60pct/)c=i; if(h[i]=="Ata_long" X)o=i} n=0; same=0; L=length(s[c])
      for(j=1;j<=L;j++){a=toupper(substr(s[c],j,1)); b=toupper(substr(s[o],j,1)); if(a!="-"&&b!="-"){n++; if(a==b)same++}}
      f=0;lz=0; for(j=1;j<=L;j++){if(substr(s[o],j,1)!="-"){if(!f)f=j; lz=j}} cf=0;cz=0; for(j=1;j<=L;j++){if(substr(s[c],j,1)!="-"){if(!cf)cf=j; cz=j}}
      printf "set %s: Eeu 60%% vs Atelerix long%s: identical %d of %d aligned positions (%.1f%%); plate columns Atelerix %d-%d, Eeu %d-%d\n",X,X,same,n,100*same/n,f,lz,cf,cz}' eeu_e2-3_long$X.aln.fa
done
rm -rf tmp t.cons search/genome.search* search/*.bnk
md5sum eeu_e2-3_long*.aln.fa
echo DONE
