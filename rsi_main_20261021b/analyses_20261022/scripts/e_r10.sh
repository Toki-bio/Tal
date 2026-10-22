#!/bin/bash
# task E: every r10_19seqs copy with verdict SINE, plus its grey copies ([grey]); not_SINE copies left out (shared flanks are never shown)
set -euo pipefail
cd /home/toki/tmp/rsi_aln_20261022 || exit 1
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
R=/home/toki/rhin/rsi_main_20261021b/run_20261021_192808
SF=r10_19seqs; rm -rf E; mkdir -p E; cd E || exit 1
awk -F'\t' -v s=$SF 'NR>1 && $2==s && $3=="assigned"{c[$4]++} END{for(k in c) print k, c[k]}' $R/step2/step2_output/loci_verdict.tsv | sort > counts.txt
# copy id cN in the BED name; names.tsv maps it back to the locus (+ [grey])
gawk -F'\t' -v s=$SF 'NR>1 && $2==s && $3=="assigned" && ($4=="SINE" || $4=="grey") {
    if (match($1, /^(.+):([0-9]+)-([0-9]+)\(([+-])\)$/, m)) { k++
        print m[1] "\t" m[2] "\t" m[3] "\tc" k "\t0\t" m[4] > "copies.bed"
        print "c" k "\t" $1 ($4=="grey" ? " [grey]" : "") > "names.tsv" } }' $R/step2/step2_output/loci_verdict.tsv
G=$R/genome.clean.fa
tab(){ gawk -F'\t' -v U="$1" '{ n = $1; sub(/::.*$/, "", n); sub(/\([+-]\)$/, "", n); print n "\t" (U ? toupper($2) : tolower($2)) }'; }
bedtools getfasta -fi $G -bed copies.bed -s -name -tab | tab 1 > el.tsv
bedtools flank -i copies.bed -g $G.fai -l 50 -r 0 -s | bedtools getfasta -fi $G -bed - -s -name -tab | tab 0 > f5.tsv
bedtools flank -i copies.bed -g $G.fai -l 0 -r 70 -s | bedtools getfasta -fi $G -bed - -s -name -tab | tab 0 > f3.tsv
head -2 el.tsv | cut -c1-60; head -1 f5.tsv
gawk -v s=$SF '/^>/{p=($1==">"s)} p' $R/consensuses.clean.fa > cons.fa
{ cat cons.fa; gawk -F'\t' '{print ">" $1 "\n" $2}' el.tsv; } > el.fa
mafft --auto --threadit 0 --thread 16 --quiet el.fa > el.aln.fa
python3 - <<'PY'
import collections
def rd(p):
    d, n = collections.OrderedDict(), None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"): n = l[1:]; d[n] = []
        elif n: d[n].append(l)
    return collections.OrderedDict((k, "".join(v)) for k, v in d.items())
aln = rd("el.aln.fa")
f5 = dict(l.rstrip("\n").split("\t") for l in open("f5.tsv"))
f3 = dict(l.rstrip("\n").split("\t") for l in open("f3.tsv"))
names = dict(l.rstrip("\n").split("\t") for l in open("names.tsv"))
assert len(f5) == len(f3) == len(names), (len(f5), len(f3), len(names))
with open("rsi_r10_19seqs_all_SINE_grey.aln.fa", "w") as fo:
    for n, s in aln.items():
        a, b = f5.get(n, ""), f3.get(n, "")
        fo.write(">%s\n%s%s%s\n" % (names.get(n, n), "-" * (50 - len(a)) + a, s.upper(), b + "-" * (70 - len(b))))
print(len(aln), "rows")
PY
cat counts.txt
