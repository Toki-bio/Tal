set -euo pipefail
cd /staging/tmp/apple_minisat/peel || exit 1
FAI=/staging/tmp/apple_sines/genome/apple.fna.fai
# loci table: id, group, chr, start, end, len, copies
awk -F'\t' 'FNR==NR{g[$1]=$2; next} ($1 in g){h=$2; split(h,a,"|"); split(a[2],b,":"); split(b[2],c,"-"); l=length($3); cp=""; if(match(h,/copies=[0-9.]+/))cp=substr(h,RSTART+7,RLENGTH-7); ll=""; if(match(h,/locuslen=[0-9]+/))ll=substr(h,RSTART+9,RLENGTH-9); print $1"\t"g[$1]"\t"b[1]"\t"c[1]"\t"c[2]"\t"l"\t"cp"\t"ll}' groups.tsv motifs.tsv > loci.tsv
awk -F'\t' 'BEGIN{OFS="\t"} FNR==NR{L[$1]=$2; next} {pos=($4+$5)/2/L[$3]; d=int(pos*10); if(d>9)d=9; print $0, sprintf("%.3f",pos)}' $FAI loci.tsv > loci_pos.tsv
echo "group n chrs total_locus_bp median_locuslen max_locuslen  rel-position deciles(0-9)"
for g in G01 G02 G03 G04 G05 G06 G07; do awk -F'\t' -v g=$g '$2==g{n++; ch[$3]=1; t+=$8; ll[n]=$8; if($8>mx)mx=$8; d=int($9*10); if(d>9)d=9; dc[d]++} END{nc=0; for(c in ch)nc++; asort(ll); s=""; for(i=0;i<10;i++)s=s" "(dc[i]+0); printf "%s n=%d chrs=%d bp=%d med=%d max=%d |%s\n", g,n,nc,t,ll[int((n+1)/2)],mx,s}' loci_pos.tsv; done
echo "-- chromosome x group counts"
awk -F'\t' '{k[$3]=1; c[$3,$2]++} END{for(x in k){printf "%s", x; split("G01 G02 G03 G04 G05 G06 G07",G," "); for(i=1;i<=7;i++)printf " %s=%d",G[i],c[x,G[i]]+0; print ""}}' loci_pos.tsv | sort
echo "-- chr lengths"; head -20 $FAI | cut -f1,2 | tr '\n' ' '
