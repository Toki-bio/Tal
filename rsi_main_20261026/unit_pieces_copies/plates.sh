#!/bin/bash
# plates of the unit pieces (unit_pieces_copies.py output in $1): one plate per unit label (all its groups, consensus rows
# first), the r8 question plate, and identity of every piece to the candidate consensuses (median per group)
set -euo pipefail
D=$(readlink -f "$1"); cd "$D" || exit 1
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
R=/home/toki/rhin/rsi_main_20261026/run_20261026_203020
mkdir -p plates
# per unit label: consensus rows of every group first, then the copies group by group
for lab in $(cut -f2 groups.tsv | tail -n +2 | sort -u); do
  fs=$(awk -F'\t' -v l="$lab" 'NR>1 && $2==l{print "pieces/"$1".fa"}' groups.tsv)
  { for f in $fs; do seqkit grep -r -n -p '\|consensus$' "$f"; done
    for f in $fs; do seqkit grep -r -n -v -p '\|consensus$' "$f"; done; } > plates/$lab.in.fa
  mafft --quiet --localpair --maxiterate 1000 --thread 16 plates/$lab.in.fa > plates/$lab.aln.fa
done
# r8 question: r8 alone, r8 piece of P2, P18's second unit in P18 and in P18del; r7 / r10 consensuses as reference rows
fs="pieces/r8_83seqs@alone.fa pieces/r8_83seqs@r1_r3_P2:2-85.fa pieces/P18_2nd@r10_r8_P18:290-426.fa pieces/P18_2nd@r10_r8_P18del:152-255.fa"
{ seqkit grep -n -p r7_133seqs -p r10_19seqs $R/consensuses.clean.fa | sed 's/^>\(.*\)$/>\1|consensus/'
  for f in $fs; do seqkit grep -r -n -p '\|consensus$' "$f"; done
  for f in $fs; do seqkit grep -r -n -v -p '\|consensus$' "$f"; done; } > plates/r8_question.in.fa
mafft --quiet --localpair --maxiterate 1000 --thread 16 plates/r8_question.in.fa > plates/r8_question.aln.fa
# identity of every piece to the candidate consensuses (best local hit, ssearch36), median per group
seqkit grep -n -p r7_133seqs -p r8_83seqs -p r10_19seqs $R/consensuses.clean.fa > refs.fa
seqkit grep -r -n -p '^r10_r8_P18:290-426\|consensus$' "pieces/P18_2nd@r10_r8_P18:290-426.fa" | sed 's/|consensus//' >> refs.fa
cat pieces/r8_83seqs@*.fa pieces/P18_2nd@*.fa | seqkit grep -r -n -v -p '\|consensus$' > q8.fa
ssearch36 -3 -m 8 -E 1 -T 16 q8.fa refs.fa > q8.m8
python3 - <<'P'
import collections, statistics
best = {}
for l in open("q8.m8"):
    f = l.split("\t")
    k = (f[0], f[1])
    if k not in best or float(f[11]) > best[k][1]:
        best[k] = (float(f[2]), float(f[11]))
refs = ["r8_83seqs", "r7_133seqs", "r10_19seqs", "r10_r8_P18:290-426"]
grp = collections.defaultdict(lambda: collections.defaultdict(list))
win = collections.defaultdict(collections.Counter)
names = [l[1:].strip() for l in open("q8.fa") if l.startswith(">")]
for n in names:
    g = n.split("|", 1)[1]
    sc = {r: best.get((n, r), (0, 0)) for r in refs}
    for r in refs:
        grp[g][r].append(sc[r][0])
    win[g][max(refs, key=lambda r: sc[r][1])] += 1
with open("r8_question_identity.tsv", "w") as fo:
    fo.write("group\tpieces\t" + "\t".join("median_ident_" + r for r in refs) + "\t" + "\t".join("best_" + r for r in refs) + "\n")
    for g in grp:
        fo.write("%s\t%d\t%s\t%s\n" % (g, len(grp[g][refs[0]]), "\t".join("%.1f" % statistics.median(grp[g][r]) for r in refs),
                                      "\t".join(str(win[g][r]) for r in refs)))
P
column -t r8_question_identity.tsv
grep -c '>' plates/*.aln.fa
