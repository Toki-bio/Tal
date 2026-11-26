#!/bin/bash
# Compare Tal eri consensuses with the Atelerix albiventris SINE consensuses of
# Genes 16:397 (PMC12026660; github.com/SystemBio-Sdut/Ata_TEs data/).
# Run in a directory holding: family_sine.fasta family_sine_info1.txt
# eri_candidate_consensuses_vs_sinebase.aln.fa
# Outputs: ours.fa their_short.fa control.tsv hits.tsv null.tsv matches.tsv
#          eri_atelerix_matched.aln.fa eri_atelerix_all_short.aln.fa
set -euo pipefail
MAXLEN=600
MINALN=80

# our ten, ungapped, first word of the header (E. europaeus candidates e*, literature ERI-*)
awk '{sub(/\r$/,"")}
/^>/{n=substr($1,2); n=(n~/^ERI/)?n"_literature":"Eeu_"n; print ">"n; next}
{gsub(/-/,""); print}' eri_candidate_consensuses_vs_sinebase.aln.fa > ours.fa

# their consensuses <= MAXLEN bp: Ata_<id>_<copies>copies; the four Figure A6 ones carry UnS-n
awk -v max=$MAXLEN 'BEGIN{FS="\t"; u["DR0187830"]="UnS-1_"; u["DR0188050"]="UnS-2_";
  u["rnd-1_family-337"]="UnS-3_"; u["rnd-1_family-322"]="UnS-4_"}
{sub(/\r$/,"")}
FNR==NR{if(FNR>1){k=$1; sub(/;$/,"",k); copies[k]=$3}; next}
/^>/{out(); id=substr($0,2); sub(/;$/,"",id); seq=""; next}
{seq=seq $0}
END{out()}
function out(){ if(id!="" && length(seq)<=max){
  if(!(id in copies)){print "no copy number for " id > "/dev/stderr"; exit 1}
  print ">Ata_" u[id] id "_" copies[id] "copies\n" seq } }' family_sine_info1.txt family_sine.fasta > their_short.fa

# null: the same sequences shuffled (fixed seed)
awk 'BEGIN{srand(1)} /^>/{print ">shuf_" substr($0,2); next}
{n=split($0,a,""); for(i=n;i>1;i--){j=int(rand()*i)+1; t=a[i];a[i]=a[j];a[j]=t}
 s=""; for(i=1;i<=n;i++)s=s a[i]; print s}' their_short.fa > their_shuf.fa

S="ssearch36 -n -m 8CB -E 100000 -T 4"
# -m 8CB comment lines are dropped; column 12 is the bit score
$S ours.fa ours.fa 2>/dev/null | awk '!/^#/' > control.tsv
$S their_short.fa ours.fa 2>/dev/null | awk '!/^#/' > hits.tsv
$S their_shuf.fa ours.fa 2>/dev/null | awk '!/^#/' > null.tsv

NULLMAX=$(awk -v m=$MINALN 'BEGIN{x=0} $4>=m && $12+0>x{x=$12+0} END{print x}' null.tsv)
echo "null max bit score (aligned >= $MINALN bp): $NULLMAX"

# best passing hit per Atelerix consensus
awk -v m=$MINALN -v t="$NULLMAX" -v OFS='\t' '
$4>=m && $12+0>t+0 { if(!($1 in b) || $12+0>b[$1]){b[$1]=$12+0;
  r[$1]=$1 OFS $2 OFS $3 OFS $4 OFS $7 OFS $8 OFS $9 OFS $10 OFS $12} c[$1 SUBSEP $2]=1 }
END{for(q in r){n=0; for(k in c){split(k,p,SUBSEP); if(p[1]==q)n++}; print r[q], n}}' hits.tsv |
 sort -t$'\t' -k9,9nr > matches.body
{ printf 'atelerix\tbest_ours\tidentity_pct\taligned_bp\tq_start\tq_end\ts_start\ts_end\tbits\tours_passing\n'; cat matches.body; } > matches.tsv
rm matches.body

# plate A: matched + UnS-1..4 + ours; plate B: all short + ours
awk 'FNR==NR{if(FNR>1)k[$1]=1; next} /^>/{n=substr($0,2); p=(n in k)||n~/UnS-/} p' matches.tsv their_short.fa | cat ours.fa - > in_matched.fa
cat ours.fa their_short.fa > in_all_short.fa
for x in matched all_short; do
  mafft --thread 8 --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --preservecase --adjustdirection --quiet in_$x.fa > eri_atelerix_$x.aln.fa
done

# check: every plate holds exactly its input sequences (or their reverse complement)
for x in matched all_short; do
  awk 'function rc(s,   i,o,c){o="";for(i=length(s);i>0;i--){c=substr(s,i,1);
        o=o (c in C?C[c]:c)}return o}
  BEGIN{split("a t c g g c t a r y y r k m m k s s w w n n b v v b d h h d",z," ");for(i=1;i<48;i+=2)C[z[i]]=z[i+1]}
  FNR==1{f++}
  /^>/{n=substr($0,2); if(f==2)sub(/^_R_/,"",n); next}
  {s=tolower($0); gsub(/-/,"",s); if(f==1)A[n]=A[n] s; else B[n]=B[n] s}
  END{for(n in A){na++; if(B[n]!=A[n] && B[n]!=rc(A[n])){print "MISMATCH " n; bad++}}
      for(n in B)nb++; print FILENAME, "in=" na, "out=" nb, "bad=" bad+0; exit (bad||na!=nb)}' in_$x.fa eri_atelerix_$x.aln.fa
done
echo "reversed by mafft:"; awk '/^>_R_/' eri_atelerix_*.aln.fa
wc -l matches.tsv
echo DONE
