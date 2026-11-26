set -euo pipefail
export PATH=/staging/tmp/apple_minisat/env/bin:$PATH
cd /staging/tmp/apple_minisat/peel || exit 1
echo "== A. shuffle control: same motifs, bases shuffled within each motif, same grouping rule (T=0.70, n>=3)"
mkdir -p ctl; cd ctl
awk -F'\t' 'BEGIN{srand(42)} {s=$3; n=length(s); for(i=1;i<=n;i++)a[i]=substr(s,i,1); for(i=n;i>1;i--){j=int(rand()*i)+1; t=a[i];a[i]=a[j];a[j]=t} o=""; for(i=1;i<=n;i++)o=o a[i]; print $1"\t"$2"\t"o}' ../motifs.tsv > motifs.tsv
awk -F'\t' '{print ">"$1"\n"$3}' motifs.tsv > subj.fa; awk -F'\t' '{print ">"$1"\n"$3$3}' motifs.tsv > q2x.fa
makeblastdb -in subj.fa -dbtype nucl -out subj >/dev/null
blastn -task blastn -query q2x.fa -db subj -strand both -dust no -word_size 7 -evalue 1e-3 -max_target_seqs 1000 -max_hsps 3 -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send sstrand' -num_threads 16 > a.tsv
awk -F'\t' 'BEGIN{OFS="\t"} {if($1==$2)next; m=($5<$6)?$5:$6; cov=$4/m; if(cov>1)cov=1; sc=$3/100*cov; k=$1 SUBSEP $2; if(!(k in b)||sc>b[k])b[k]=sc} END{for(k in b){split(k,x,SUBSEP); print x[1],x[2],b[k]}}' a.tsv > d.tsv
awk -F'\t' '{k=$1<$2?$1 SUBSEP $2:$2 SUBSEP $1; if(!(k in v)||$3<v[k])v[k]=$3} END{for(k in v){split(k,x,SUBSEP);print x[1],x[2],v[k]}}' OFS='\t' d.tsv > s.tsv
awk -F'\t' -v T=0.70 'function f(x){while(p[x]!=x)x=p[x]=p[p[x]];return x} FNR==NR{p[$1]=$1;ids[++n]=$1;next}{if($3>=T){a=f($1);b=f($2);if(a!=b)p[a]=b}} END{for(i=1;i<=n;i++)c[f(ids[i])]++; big=0;mx=0; for(r in c){if(c[r]>=3)big++; if(c[r]>mx)mx=c[r]} print "shuffled: groups with n>=3:",big," largest:",mx}' motifs.tsv s.tsv
cd ..
echo "== B. leftovers: motifs not in a group of >=3"
awk -F'\t' 'FNR==NR{n[$2]=$3; next} {print $1"\t"n[$2]}' groups.tsv groups.tsv | awk -F'\t' '$2<3{print $1}' > left.ids
wc -l < left.ids
# best score of each leftover to any member of groups G01..G07 (using pair_sym), and its length
awk -F'\t' 'FNR==NR{g[$1]=$2; sz[$1]=$3; next} {a=$1;b=$2; if(sz[a]>=3&&sz[b]<3){x=b;y=a} else if(sz[b]>=3&&sz[a]<3){x=a;y=b} else next; if($3>best[x]){best[x]=$3; bg[x]=g[y]}} END{for(i in best)print i"\t"best[i]"\t"bg[i]}' groups.tsv pair_sym.tsv > left_best.tsv
awk -F'\t' 'FNR==NR{L[$1]=length($3); next} {print $1,L[$1]}' motifs.tsv left.ids > /dev/null
awk -F'\t' 'FILENAME==ARGV[1]{L[$1]=length($3); next} FILENAME==ARGV[2]{l[$1]=1; next} FILENAME==ARGV[3]{b[$1]=$2; g[$1]=$3; next} END{ for(i in l){ sc=(i in b)?b[i]:0; bin=(sc>=0.6)?">=0.60":(sc>=0.5)?"0.50-0.59":(sc>=0.4)?"0.40-0.49":"<0.40/none"; cnt[bin]++; ln[bin]=ln[bin]" "L[i]} for(k in cnt)print k, cnt[k]}' motifs.tsv left.ids left_best.tsv | sort
echo "leftover lengths:"; awk -F'\t' 'FILENAME==ARGV[1]{L[$1]=length($3); next}{print L[$1]}' motifs.tsv left.ids | sort -n | uniq -c | awk '{printf "%s(x%s) ",$2,$1}'; echo
echo "== C. G02 (170 bp, interior?) positions"
awk -F'\t' '$2=="G02"{printf "%s:%d-%d (rel %.2f, %d bp)\n",$3,$4,$5,$9,$8}' loci_pos.tsv | sort
echo "== D. G03 and others: spans"
for g in G03 G06 G07 G05 G04; do awk -F'\t' -v g=$g '$2==g{if(!(($3) in mn)||$4<mn[$3])mn[$3]=$4; if($5>mx[$3])mx[$3]=$5; n[$3]++} END{for(c in n)printf "%s %s n=%d span %d-%d (%.0f kb)\n",g,c,n[c],mn[c],mx[c],(mx[c]-mn[c])/1000}' loci_pos.tsv; done | sort
