set -euo pipefail
export PATH=/staging/tmp/apple_minisat/env/bin:$PATH
cd /staging/tmp/apple_minisat/peel || exit 1
: > consensus_all.fa; printf "group\tn_aligned\tref\tmedian_member_len\tcons_len\tmean_id_to_cons\tmin_id\tcols_kept\tcols_total\n" > consensus_stats.tsv
for d in grp/G*; do g=$(basename $d); [ -f $d/aln.fa ] || continue
awk -v g=$g -v out=$d/cons.fa 'function flush(){ if(id!=""){seq[++n]=s; nm[n]=id} } /^>/{flush(); id=substr($0,2); s=""; next}{s=s toupper($0)} END{flush(); W=length(seq[1]);
 for(c=1;c<=W;c++){delete cnt; ng=0; for(i=1;i<=n;i++){b=substr(seq[i],c,1); if(b!="-")cnt[b]++; else ng++} best="";bm=0; for(b in cnt) if(cnt[b]>bm||(cnt[b]==bm&&b<best)){bm=cnt[b];best=b}
   keep[c]=((n-ng)/n>=0.5); cb[c]=best}
 cons=""; kept=0; for(c=1;c<=W;c++) if(keep[c]){cons=cons cb[c]; kept++}
 print ">cons_"g"_n"n"\n"cons > out
 # identity of each member to consensus over kept columns (non-gap member bases only)
 tot=0; mn=1; for(i=1;i<=n;i++){m=0;t=0; k=0; for(c=1;c<=W;c++) if(keep[c]){k++; b=substr(seq[i],c,1); if(b!="-"){t++; if(b==cb[c])m++}} id=(t>0?m/t:0); tot+=id; if(id<mn)mn=id; lens[i]=t}
 asort(lens); printf "%s\t%d\t-\t%d\t%d\t%.3f\t%.3f\t%d\t%d\n", g, n, lens[int((n+1)/2)], length(cons), tot/n, mn, kept, W}' $d/aln.fa >> consensus_stats.tsv
cat $d/cons.fa >> consensus_all.fa
done
sed -i 's/\t-\t/\t'"$(true)"'\t/' consensus_stats.tsv 2>/dev/null || true
cat consensus_stats.tsv
# cross-group relation of consensuses: doubled query vs subject, both strands
awk '/^>/{if(h!="")print h"\t"s; h=substr($0,2); s=""; next}{s=s $0} END{print h"\t"s}' consensus_all.fa > cons.tsv
awk -F'\t' '{print ">"$1"\n"$2}' cons.tsv > csubj.fa; awk -F'\t' '{print ">"$1"\n"$2$2}' cons.tsv > cq2x.fa
makeblastdb -in csubj.fa -dbtype nucl -out csubj >/dev/null
echo "-- cross-consensus (best HSP: pident, alnlen, strand)"
blastn -task blastn -query cq2x.fa -db csubj -strand both -dust no -word_size 7 -evalue 1e-2 -outfmt '6 qseqid sseqid pident length qlen slen sstrand' | awk '$1!=$2' | sort -k1,1 -k2,2 -k4,4nr | awk '!k[$1,$2]++' | awk '{printf "%s vs %s pid=%.1f aln=%d slen=%d %s\n",$1,$2,$3,$4,$6,$7}'
echo "-- self-periodicity of each consensus (doubled vs itself, plus, excluding the trivial diagonal)"
