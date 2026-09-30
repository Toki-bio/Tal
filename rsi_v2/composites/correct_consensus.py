"""correct_consensus.py - rebuild each candidate consensus from its flanked alignment (flanked/*.flanked.aln.fa).

Row 1 of the input is the candidate consensus as built by flankscan stage 5, rows 2.. the copies with 100 bp
lowercase flanks. Correction, from the copies only:
  - inside the candidate's span: a column whose top base is carried by >= SHARE of all copies and differs
    from the candidate base replaces it (a candidate gap there becomes that base);
  - past each end: walk outward column by column; a column occupied by < SKIP of the copies is passed over,
    a column whose top base is carried by >= SHARE of all copies extends the consensus, the first column
    that is neither stops the walk.
Output: corrected/<name>.aln.fa = corrected row (uppercase, name <cand>_corrected) + the candidate as built
+ the copies, same columns (nothing realigned); corrected/<name>.fa; corrected/changes.tsv.
"""
import glob
import os
from collections import Counter

SHARE, SKIP = 0.6, 0.3
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "corrected")
os.makedirs(OUT, exist_ok=True)


def read(p):
    recs = []
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            recs.append([line[1:], ""])
        elif recs:
            recs[-1][1] += line
    return recs


rows = ["name\tlen_before\tlen_after\tadded5\tadded3\tseq5\tseq3\tsubstitutions"]
for p in sorted(glob.glob(os.path.join(HERE, "flanked", "*.flanked.aln.fa"))):
    name = os.path.basename(p)[: -len(".flanked.aln.fa")]
    recs = read(p)
    cname, cons = recs[0]
    cop = [s for _, s in recs[1:]]
    N, L = len(cop), len(cons)

    def call(i):
        bases = [s[i].upper() for s in cop if s[i] != "-"]
        occ = len(bases) / N
        if not bases:
            return occ, "-", 0.0
        b, c = Counter(bases).most_common(1)[0]
        return occ, b, c / N

    cols = [i for i, c in enumerate(cons) if c != "-"]
    a, b = cols[0], cols[-1]
    new = list(cons.upper())
    subs = 0
    for i in range(a, b + 1):
        occ, top, sh = call(i)
        if sh >= SHARE and top != new[i]:
            new[i] = top
            subs += 1
    ext = {}
    for side, rng in (("5", range(a - 1, -1, -1)), ("3", range(b + 1, L))):
        got = []
        for i in rng:
            occ, top, sh = call(i)
            if occ < SKIP:
                continue
            if sh >= SHARE:
                new[i] = top
                got.append(i)
            else:
                break
        ext[side] = sorted(got)
    seq5 = "".join(new[i] for i in ext["5"])
    seq3 = "".join(new[i] for i in ext["3"])
    row = "".join(c if c != "-" else "-" for c in new)
    ungapped = row.replace("-", "")
    with open(os.path.join(OUT, name + ".aln.fa"), "w") as fh:
        fh.write(">%s_corrected\n%s\n>%s_as_built\n%s\n" % (name, row, name, cons))
        for n, s in recs[1:]:
            fh.write(">%s\n%s\n" % (n, s))
    with open(os.path.join(OUT, name + ".fa"), "w") as fh:
        fh.write(">%s_corrected\n%s\n" % (name, ungapped))
    rows.append("%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d" % (name, len(cons.replace("-", "")), len(ungapped),
                                                    len(seq5), len(seq3), seq5 or "-", seq3 or "-", subs))
open(os.path.join(OUT, "changes.tsv"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
