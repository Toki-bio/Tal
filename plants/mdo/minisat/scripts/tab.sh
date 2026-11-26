set -euo pipefail
cd /staging/tmp/apple_minisat || exit 1
R=run_fam/cross/cross_clusters
# group -> family name
awk -F'\t' '{split($1,a,"_"); print a[1]"\t"$1}' cfg/fam/families.tsv > cfg/g2fam.tsv
# locus key -> group from the pass-1 peel
awk -F'\t' '{print $3":"$4"-"$5"\t"$2}' peel/loci.tsv > cfg/locus2group.tsv
for f in $R/clust_*; do n=${f##*_};
  awk -v n=$n 'FILENAME==ARGV[1]{g2f[$1]=$2; next} FILENAME==ARGV[2]{lg[$1]=$2; next} /^>/{h=substr($0,2); if(match(h,/^[A-Za-z0-9]+\|[^|:]+:[0-9]+-[0-9]+/)){k=substr(h,RSTART,RLENGTH); sub(/^[^|]*\|/,"",k); tot++; if(k in lg){c[lg[k]]++}else un++} else {split(h,a," "); r=a[1]; for(g in g2f) if(g2f[g]==r){c[g]++; ref=r}; tot++}} END{best="";bm=0; for(g in c)if(c[g]>bm){bm=c[g];best=g} nf=0; for(g in c)nf++; printf "%s\t%s\t%d\t%d\t%d\t%s\n", n, (best==""?"-":g2f[best]), bm, tot, nf, (un?"unmapped="un:"")}' cfg/g2fam.tsv cfg/locus2group.tsv $f
done | sort -n > cfg/cluster_map_draft.tsv
cat cfg/cluster_map_draft.tsv
