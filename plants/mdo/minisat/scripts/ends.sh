set -euo pipefail
cd /staging/tmp/apple_minisat/peel || exit 1
FAI=/staging/tmp/apple_sines/genome/apple.fna.fai
C=/staging/tmp/apple_minisat/run_fam/large/mdo/NC_091675.1/candidates.tsv
# all loci: group family, chr, mid, len
awk -F'\t' 'FNR==NR{L[$1]=$2; next} {print $2"\t"$3"\t"($4+$5)/2"\t"L[$3]}' $FAI loci.tsv | awk -F'\t' '$1 ~ /^G0[1-7]$/' > L1.tsv
awk -F'\t' 'FNR==NR{L[$1]=$2; next} NR>1{f=substr($5,1,3); print f"\tNC_091675.1\t"($1+$2)/2"\t"L["NC_091675.1"]}' $FAI $C >> L1.tsv
echo "family n  | within 1 Mb of a chromosome end | expected if uniform (2*1Mb/L, mean over loci) | median distance to nearest end (Mb) | n chromosomes | n loci per chromosome (max)"
for g in G01 G02 G03 G04 G05 G06 G07; do awk -F'\t' -v g=$g '$1==g{n++; d=$3; e=$4-$3; if(e<d)d=e; dist[n]=d; if(d<=1000000)t++; ex+=2000000/$4; ch[$2]++} END{asort(dist); nc=0; mx=0; for(c in ch){nc++; if(ch[c]>mx)mx=ch[c]} printf "%s n=%d | %d (%.0f%%) | %.1f%% | %.1f | %d chr | max %d\n", g,n,t,100*t/n,100*ex/n,dist[int((n+1)/2)]/1e6,nc,mx}' L1.tsv; done
echo "-- terminal loci per family: chr, distance to nearest end (kb), which end"
awk -F'\t' '{d=$3; w="start"; e=$4-$3; if(e<d){d=e; w="end"} if(d<=1000000) printf "%s\t%s\t%s\t%.0f\n",$1,$2,w,d/1000}' L1.tsv | sort -k1,1 -k2,2 -k3,3 -k4,4n | awk -F'\t' '{k=$1" "$2" "$3; n[k]++; if(!(k in mn))mn[k]=$4; mx[k]=$4} END{for(k in n)printf "%s: %d loci, %d-%d kb from the end\n",k,n[k],mn[k],mx[k]}' | sort
