#!/usr/bin/env bash
# Made the seed plate of 2026-11-26 on Biotech2024 (run as: bash make_seed_plate.sh).
set -euo pipefail
D=/staging/tmp/apple_sines
mkdir -p $D/publish
cd $D/publish || exit 1
export TMPDIR=$D/tmp
E=$D/env/bin
S=$D/anno_full/out/Seed_SINE.fa
# rename: SINE_k -> mdo_seedk
awk '/^>/{split(substr($1,2),a,"_"); print ">mdo_seed"a[2]; next}{print}' $S > mdo_seeds.fa
$E/ssearch36 -Q -n -T 32 -E 1e-3 -m 8 mdo_seeds.fa $D/repbase_plant84.fa > mdo_vs_bank.m8
awk -F'\t' -v OFS='\t' '
FILENAME==ARGV[1]{if($0~/^>/){n=substr($0,2)} else L[n]+=length($0); next}
FILENAME==ARGV[2]{if($3>=80 && (($8>$7?$8-$7:$7-$8)+1)>=0.8*L[$1] && !($1 in M)){M[$1]=$2} next}
/^>/{h=substr($0,2); split(h,w," "); split(w[1],a,"_"); name="mdo_seed"a[2]
  m=split(h,p,"|"); split(p[m],t," "); loc=w[2]" "w[3]; sub(/\|.*/,"",loc); bc=""
  if(match(h,/blast_count:[0-9]+/)) bc=substr(h,RSTART+12,RLENGTH-12)
  print name,w[1],loc,L[name],t[1],bc,(name in M?M[name]:"-")}
BEGIN{print "name","annosine_id","annosine_location","length_bp","head_tag","blast_count","apple_repbase_match_80_80"}
' mdo_seeds.fa mdo_vs_bank.m8 $S > mdo_seeds.names.tsv
rm mdo_vs_bank.m8
$E/mafft --quiet --localpair --maxiterate 1000 --threadit 0 --thread 32 --adjustdirection --reorder mdo_seeds.fa > mdo_consensuses.aln.fa
# check: every seed present once, bases equal to the seed or its reverse complement
awk 'function rc(s,   i,o,c){o="";for(i=length(s);i>0;i--){c=substr(s,i,1); o=o (c in C?C[c]:"N")}return o}
BEGIN{C["A"]="T";C["T"]="A";C["G"]="C";C["C"]="G"}
/^>/{n=substr($1,2); r=sub(/^_R_/,"",n); if(NR!=FNR){N++; R+=r}; next}
NR==FNR{s[n]=s[n] toupper($0); next}
{x=toupper($0); gsub(/-/,"",x); a[n]=a[n] x; if(W==0)W=0; len[n]+=length($0)}
END{for(k in a){if(a[k]==s[k])f++; else if(a[k]==rc(s[k]))v++; else bad++}; for(k in len){c=len[k]; break}
 print "plate_sequences",N; print "same_strand",f+0; print "reverse_complemented",v+0; print "mismatch",bad+0; print "columns",c}' mdo_seeds.fa mdo_consensuses.aln.fa > counts.txt
awk -F'\t' 'NR>1{n++; L[n]=$4; t[$5]++; if($7!="-"){m++; h[$7]++}}
END{asort(L); print "seeds",n; print "len_min",L[1]; print "len_median",L[int((n+1)/2)]; print "len_max",L[n]; for(k in t) print "tag_"k,t[k]; print "matched",m; for(k in h) print "match_"k,h[k]}' mdo_seeds.names.tsv >> counts.txt
cat counts.txt; head -4 mdo_seeds.names.tsv; ls -l
