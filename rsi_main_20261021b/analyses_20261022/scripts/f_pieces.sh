#!/bin/bash
# task F: every unit piece of the Rhin consensuses: stand-alone units (r1 ... r10) and every stretch cut from a composite
# (step2_output/composite_units.tsv), named unit@consensus:start-end; unlabelled composite stretches >= 30 bp as unassigned@...
set -euo pipefail
cd /home/toki/tmp/rsi_aln_20261022 || exit 1
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
R=/home/toki/rhin/rsi_main_20261021b/run_20261021_192808
rm -rf F; mkdir -p F/by_unit; cd F || exit 1
python3 - $R/consensuses.clean.fa $R/step2/step2_output/composite_units.tsv <<'PY'
import csv, collections, sys
seq, n = {}, None
for l in open(sys.argv[1]):
    l = l.strip()
    if l.startswith(">"): n = l[1:].split()[0]; seq[n] = ""
    elif n: seq[n] += l.upper()
by = collections.defaultdict(list)
for r in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
    c = r["consensus"]
    if c.startswith("MEG"):
        continue
    a, b = int(r["start"]), int(r["end"])
    u = r["unit"] if r["unit"] != "-" else "unassigned"
    if u == "unassigned" and b - a + 1 < 30:
        continue
    nm = ("%s@%s" % (u, c)) if r["kind"] == "unit" else ("%s@%s:%d-%d" % (u, c, a, b))
    by[u].append((nm, seq[c][a - 1:b], r["kind"], r["identity"]))
order = sorted(by, key=lambda u: (u == "unassigned", int(u[1:].split("_")[0]) if u[0] == "r" and u[1:].split("_")[0].isdigit() else 99, u))
with open("all_pieces.fa", "w") as fo:
    for u in order:
        for nm, s, k, idt in sorted(by[u], key=lambda x: x[2] != "unit"):
            fo.write(">%s\n%s\n" % (nm, s))
for u in order:
    if len(by[u]) < 2:
        continue
    with open("by_unit/%s.fa" % u, "w") as fo:
        for nm, s, k, idt in sorted(by[u], key=lambda x: x[2] != "unit"):
            fo.write(">%s\n%s\n" % (nm, s))
for u in order:
    print(u, len(by[u]), " ".join(x[0] for x in by[u]))
PY
mafft --localpair --maxiterate 1000 --threadit 0 --thread 16 --quiet all_pieces.fa > rsi_all_unit_pieces.aln.fa
for f in by_unit/*.fa; do mafft --localpair --maxiterate 1000 --threadit 0 --thread 8 --quiet $f > ${f%.fa}.aln.fa; rm -f $f; done
grep -c ">" rsi_all_unit_pieces.aln.fa; ls by_unit
