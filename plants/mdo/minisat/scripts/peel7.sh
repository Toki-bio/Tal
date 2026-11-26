set -euo pipefail
export PATH=/staging/tmp/apple_minisat/env/bin:$PATH
cd /staging/tmp/apple_minisat/peel || exit 1
cat consensus_all.fa | awk '/^>/{print; next}{print}' | cut -c1-250
echo "== pairwise ssearch36 (doubled query, both strands) best bits / identity / overlap"
awk '/^>/{if(h!="")print h"\t"s; h=substr($0,2); s=""; next}{s=s $0} END{print h"\t"s}' consensus_all.fa > c.tsv
n=$(wc -l < c.tsv)
mkdir -p ss; rm -f ss/*
awk -F'\t' '{split($1,a,"_"); print ">"a[2]"\n"$2 > "ss/"a[2]".fa"; print ">"a[2]"\n"$2$2 > "ss/"a[2]".2x.fa"}' c.tsv
for a in $(ls ss/*.2x.fa | xargs -n1 basename | sed 's/.2x.fa//'); do for b in $(ls ss/G??.fa | xargs -n1 basename | sed 's/.fa//'); do [ $a = $b ] && continue; 
 ssearch36 -m 8C -n -3 -E 1e-3 -T 4 ss/$a.2x.fa ss/$b.fa 2>/dev/null | awk -v a=$a -v b=$b '!/^#/ && NF>=12{printf "%s vs %s: %s%% over %s, bits %s\n",a,b,$3,$4,$12; exit}'; done; done | head -60
