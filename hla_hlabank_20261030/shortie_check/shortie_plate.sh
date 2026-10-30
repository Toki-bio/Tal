#!/bin/bash
# shortie check: 100 random SINE-verdict copies each of consensus_47seqs and consensus_28seqs with 150 bp genome flanks
# (element orientation), aligned with the shorties and the long consensuses whose head they match
set -euo pipefail
cd /home/toki/tmp/hla_cmp || exit 1
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
R=/home/toki/rhin/hla_hlabank_20261030/run_20261029_220442; B=/home/toki/rhin/hla_bank_20261030/bank30_hla.fa
G=$R/genome.clean.fa
for c in consensus_47seqs consensus_28seqs; do
  zcat -f $R/step2/step2_output/loci_verdict.tsv* | awk -F'\t' -v c=$c '$2==c && $3=="assigned" && $4=="SINE"{print $1}' \
    | shuf --random-source=<(yes) -n 100 \
    | awk '{n=$1; st=substr(n,length(n)-1,1); sub(/\([+-]\)$/,"",n); i=match(n,/:[0-9]+-[0-9]+$/);
        split(substr(n,i+1),a,"-"); print substr(n,1,i-1)"\t"a[1]"\t"a[2]"\t"$1"\t0\t"st}' > $c.sel.bed
  bedtools slop -i $c.sel.bed -g $G.fai -b 150 | bedtools getfasta -fi $G -bed - -s -nameOnly | sed 's/([+-])$//' > $c.flanked.fa
  # copy in upper case, flanks in lower case
  python3 - $c.flanked.fa $c.sel.bed <<'P'
import sys
fa, bed = sys.argv[1:3]
L = {}
for l in open(bed):
    f = l.split("\t"); L[f[3]] = int(f[2]) - int(f[1])
name, seq, out = None, [], []
def flush():
    if name:
        s = "".join(seq); n = L[name]; a = (len(s) - n) // 2 if len(s) >= n + 300 else None
        if a is None: a = 150   # clipped at a contig end: keep 150 bp marking
        out.append(">%s\n%s%s%s\n" % (name, s[:a].lower(), s[a:a + n].upper(), s[a + n:].lower()))
for l in open(fa):
    l = l.strip()
    if l.startswith(">"):
        flush(); name, seq = l[1:], []
    else:
        seq.append(l)
flush()
open(fa, "w").write("".join(out))
P
done
seqkit grep -n -p consensus_47seqs -p consensus_28seqs -p consensus_10seqs_b -p consensus_29seqs -p consensus_12seqs $B > refs.fa
cat refs.fa consensus_47seqs.flanked.fa consensus_28seqs.flanked.fa > shortie_in.fa
mafft --quiet --localpair --maxiterate 1000 --thread 16 --preservecase shortie_in.fa > hla_shorties_flank150.aln.fa
grep -c '^>' hla_shorties_flank150.aln.fa
