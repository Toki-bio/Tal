#!/bin/bash
# task D: genomic Ile tRNA genes of rsi outside SINE copies, + 100 bp downstream, against the SINE consensuses
set -euo pipefail
cd /home/toki/tmp/rsi_aln_20261022 || exit 1
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
R=/home/toki/rhin/rsi_main_20261021b/run_20261021_192808
G=/home/toki/rhin/genomes/rsi.fna
rm -rf D; mkdir -p D; cd D || exit 1
# Ile tRNA genes of the RefSeq annotation (GCF_057876585.1_RSI-T2T genomic.gff: tRNA features, product tRNA-Ile)
zcat ../refseq.gff.gz | gawk -F'\t' -v OFS='\t' '$3 == "tRNA" && /product=tRNA-Ile/ {
      match($9, /ID=rna-([^;]+)/, m); match($9, /anticodon ([ACGU]{3})\)/, x)
      print $1, $4 - 1, $5, "Ile-" x[1] "_" m[1] "_" $1 ":" $4 "-" $5 "(" $7 ")", 0, $7 }' | sort -k1,1 -k2,2n > ile_genes.bed
echo "RefSeq Ile tRNA genes: $(wc -l < ile_genes.bed)"
# SINE hits of the run (sanitized names: _ -> @U@), back to genome names
HITS=$R/genome.clean_step1/searches/all_hits.labeled.bed
[ -s $HITS ] || HITS=$(ls $R/results/all_hits.labeled.bed)
# SINE copies: hits longer than 110 bp (a tRNA-headed consensus also hits genuine tRNA genes, over ~70-90 bp)
gawk -F'\t' -v OFS='\t' '{ gsub(/@U@/, "_", $1); if ($3 - $2 > 110) print $1, $2, $3 }' $HITS | sort -k1,1 -k2,2n > sine_hits.bed
bedtools intersect -v -a ile_genes.bed -b sine_hits.bed > ile_free.bed
bedtools intersect -u -a ile_genes.bed -b sine_hits.bed > ile_in_sine.bed
echo "outside SINE copies: $(wc -l < ile_free.bed); inside SINE copies: $(wc -l < ile_in_sine.bed)"
samtools faidx $G 2> /dev/null || true
# gene + 100 bp downstream (strand-aware); gene part upper case, downstream lower case
bedtools slop -i ile_free.bed -g $G.fai -l 0 -r 100 -s > ile_free_dn.bed
bedtools getfasta -fi $G -bed ile_free_dn.bed -s -name -tab | gawk -F'\t' '{ n = $1; sub(/::.*$/, "", n); sub(/\([+-]\)$/, "", n); print n "\t" $2 }' > genes_dn.tsv
python3 - "$R/consensuses.clean.fa" <<'PY'
import collections, re, sys
seq, n = collections.OrderedDict(), None
for l in open(sys.argv[1]):
    l = l.strip()
    if l.startswith(">"): n = l[1:].split()[0]; seq[n] = ""
    elif n: seq[n] += l.upper()
glen = {}
with open("genes_dn.fa", "w") as fo:
    for l in open("genes_dn.tsv"):
        nm, s = l.rstrip("\n").split("\t")
        m = re.search(r":(\d+)-(\d+)(\([+-]\))?$", nm); L = int(m.group(2)) - int(m.group(1)) + 1   # GFF 1-based gene span
        glen[nm] = L
        fo.write(">%s\n%s%s\n" % (nm, s[:L].upper(), s[L:].lower()))
open("glen.tsv", "w").write("".join("%s\t%d\n" % kv for kv in glen.items()))
# the Rhin consensuses whose head is tRNA-Ile: r-units and composites
with open("sines.fa", "w") as fo:
    for k, s in seq.items():
        if k.startswith("r"):
            fo.write(">%s\n%s\n" % (k, s))
PY
# every SINE consensus against every gene+100 (forward strand of both: the genes are already in tRNA orientation)
ssearch36 -m 8 -3 -E 1e-3 -T 16 sines.fa genes_dn.fa > hits.m8 2> /dev/null
# null: the same with the downstream 100 bp shuffled (gene kept), 5 shuffles per gene
python3 - <<'PY'
import random
random.seed(7)
g, n = {}, None
for l in open("genes_dn.fa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:]; g[n] = ""
    else: g[n] += l
L = dict(l.split("\t") for l in open("glen.tsv"))
with open("genes_dn_shuf.fa", "w") as fo:
    for k, s in g.items():
        a, b = s[:int(L[k])], list(s[int(L[k]):])
        for i in range(5):
            random.shuffle(b); fo.write(">%s_shuf%d\n%s%s\n" % (k, i + 1, a, "".join(b)))
PY
ssearch36 -m 8 -3 -E 1e-3 -T 16 sines.fa genes_dn_shuf.fa > hits_shuf.m8 2> /dev/null
python3 - <<'PY'
import collections, statistics
L = {k: int(v) for k, v in (l.split("\t") for l in open("glen.tsv"))}
def best(path, strip=False):
    b = {}
    for l in open(path):
        f = l.split("\t"); q, s = f[0], f[1]
        g = s.rsplit("_shuf", 1)[0] if strip else s
        key = (q, s)
        if key not in b or float(f[11]) > b[key][0]:
            b[key] = (float(f[11]), float(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8]), int(f[9]), g)
    return b
real, shuf = best("hits.m8"), best("hits_shuf.m8", True)
ext = lambda h: h[6] - L[h[7]]          # bp past the tRNA 3' end reached by the alignment on the genomic sequence
rows = []
for (q, s), h in real.items():
    sh = [ext(v) for (q2, s2), v in shuf.items() if q2 == q and v[7] == s]
    rows.append((q, s, h[0], h[1], h[2], h[3], h[4], h[5], h[6], ext(h), max(sh) if sh else 0))
with open("extension.tsv", "w") as fo:
    fo.write("sine\tgene\tbits\tident\taln_len\tq_start\tq_end\ts_start\ts_end\tbp_past_tRNA_end\tshuffled_max_bp_past\n")
    for r in sorted(rows, key=lambda r: (r[0], -r[9])):
        fo.write("\t".join(map(str, r)) + "\n")
per_sine = collections.defaultdict(list)
for r in rows: per_sine[r[0]].append(r)
print("sine\tgenes_hit\tmedian_bp_past\tmax_bp_past (gene)\tshuffled_median_max")
for q, rs in sorted(per_sine.items()):
    rs.sort(key=lambda r: -r[9])
    print("%s\t%d\t%s\t%d (%s, %.1f bits)\t%s" % (q, len(rs), statistics.median(r[9] for r in rs), rs[0][9], rs[0][1], rs[0][2],
          statistics.median(r[10] for r in rs)))
PY

# alignment: SINE consensuses (Rhin) and the free Ile genes + 100 bp downstream (gene upper case, downstream lower case)
cat sines.fa genes_dn.fa > d_input.fa
mafft --localpair --maxiterate 1000 --threadit 0 --thread 16 --quiet --preservecase d_input.fa > rsi_ile_tRNA_genes_dn100_vs_SINEs.aln.fa
grep -c ">" rsi_ile_tRNA_genes_dn100_vs_SINEs.aln.fa
