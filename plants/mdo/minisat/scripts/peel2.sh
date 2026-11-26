set -euo pipefail
export PATH=/staging/tmp/apple_minisat/env/bin:$PATH
cd /staging/tmp/apple_minisat/peel || exit 1
# pair score: best over HSPs of pident/100 * min(1, alnlen/min(slen,qlen)) ; symmetrised by min of the two directions
awk -F'\t' 'BEGIN{OFS="\t"} {if($1==$2)next; m=($5<$6)?$5:$6; cov=$4/m; if(cov>1)cov=1; sc=$3/100*cov; k=$1 SUBSEP $2; if(!(k in b)||sc>b[k]){b[k]=sc; st[k]=$11}} END{for(k in b){split(k,a,SUBSEP); print a[1],a[2],b[k],st[k]}}' allvsall.tsv > pair_dir.tsv
awk -F'\t' 'BEGIN{OFS="\t"} {k=$1<$2?$1 SUBSEP $2:$2 SUBSEP $1; if(!(k in v)){v[k]=$3; s[k]=$4} else if($3<v[k]){v[k]=$3; s[k]=$4}} END{for(k in v){split(k,a,SUBSEP); print a[1],a[2],v[k],s[k]}}' pair_dir.tsv > pair_sym.tsv
wc -l pair_sym.tsv
for T in 0.90 0.80 0.70 0.60; do
awk -F'\t' -v T=$T 'function f(x){while(p[x]!=x)x=p[x]=p[p[x]];return x} FNR==NR{id[FNR]=$1; len[$1]=length($3); p[$1]=$1; next} {if($3>=T){a=f($1);b=f($2);if(a!=b)p[a]=b}} END{for(i in len){r=f(i);n[r]++; if(!(r in mn)||len[i]<mn[r])mn[r]=len[i]; if(len[i]>mx[r])mx[r]=len[i]} c=0; printf "T=%s components:",T; for(r in n){if(n[r]>=2)c++} ; s=""; for(r in n)s=s" "n[r]"("mn[r]"-"mx[r]")"; print c" with n>=2;"; print s}' motifs.tsv pair_sym.tsv
done
