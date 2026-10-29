"""Plurality consensus of chosen rows of a plate (columns occupied in >= half of the rows; ties -> first in ACGT).
Usage: rows_consensus.py <plate.aln.fa> <name> <first row name> <last row name>   (rows taken in plate order, inclusive)"""
import collections, sys

plate, name, first, last = sys.argv[1:5]
rows, n = [], None
for l in open(plate):
    l = l.strip()
    if l.startswith(">"):
        n = l[1:].split()[0]; rows.append([n, ""])
    elif rows:
        rows[-1][1] += l.upper()
names = [r[0] for r in rows]
i, j = names.index(first), names.index(last)
sel = [r[1] for r in rows[i:j + 1]]
cons = []
for c in range(len(sel[0])):
    col = [s[c] for s in sel if s[c] in "ACGT"]
    if len(col) * 2 >= len(sel):
        cnt = collections.Counter(col)
        cons.append(max("ACGT", key=lambda b: cnt[b]))
print(">%s rows %s..%s (%d rows) of %s" % (name, first, last, len(sel), plate.rsplit("/", 1)[-1]))
print("".join(cons))
