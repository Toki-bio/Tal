set -euo pipefail
cd /staging/tmp/apple_minisat/run_fam/catalogue || exit 1
FAI=/staging/tmp/apple_sines/genome/apple.fna.fai
awk -F'\t' 'FILENAME==ARGV[1]{L[$1]=$2; next} FILENAME==ARGV[2]{f[$1]=$8; st[$1]=$7; next} {id=$4; len=$3-$2; mid=($2+$3)/2; d=mid; e=L[$1]-mid; w="start"; if(e<d){d=e;w="end"} print id"\t"f[id]"\t"$1"\t"$2"\t"$3"\t"len"\t"st[id]"\t"d"\t"w"\t"L[$1]}' $FAI MASTER.final_assignment.clean.tsv MASTER.loci.segmented.bed > /staging/tmp/apple_minisat/final_loci.tsv
cd /staging/tmp/apple_minisat
echo "== final catalogue: loci, total bp, chromosomes, within 1 Mb of an end (expected if uniform), median/max locus length"
for g in $(cut -f2 final_loci.tsv | sort -u); do awk -F'\t' -v g=$g '$2==g{n++; bp+=$6; ch[$3]++; if($8<=1000000)t++; ex+=2000000/$10; ll[n]=$6; if($6>mx)mx=$6} END{asort(ll); nc=0; for(c in ch)nc++; printf "%s\tloci=%d\tbp=%d\tchr=%d\tnear_end=%d (%.0f%%; uniform %.1f%%)\tmedian=%d\tmax=%d\n",g,n,bp,nc,t,100*t/n,100*ex/n,ll[int((n+1)/2)],mx}' final_loci.tsv; done
echo "== assignment status"; cut -f7 final_loci.tsv | sort | uniq -c
echo "== total catalogued bp, fraction of genome"; awk -F'\t' '{b+=$6} END{printf "%d bp = %.3f%% of 651,332,481\n", b, 100*b/651332481}' final_loci.tsv
echo "== longest loci"; sort -t$'\t' -k6,6nr final_loci.tsv | head -8 | cut -f1-6,8,9
echo "== per chromosome x family"; awk -F'\t' '{c[$3,$2]++; k[$3]=1; f[$2]=1} END{for(x in k){printf "%s",x; n=asorti(f,fs); for(i=1;i<=n;i++)printf " %s=%d",fs[i],c[x,fs[i]]+0; print ""}}' final_loci.tsv | sort
