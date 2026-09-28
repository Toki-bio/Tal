#!/usr/bin/env python3
"""groups.py ALN OUT.tsv [cut...] - pairwise identity of the rebuilt consensuses (and references) over the
columns both have a letter (>= 60 shared, and >= 60 % of the shorter one), then single-linkage groups at
each identity cut. Output: per row, its group at each cut; printed: the groups at each cut."""
import sys

names, seqs, cur = [], [], []
for l in open(sys.argv[1]):
    l = l.rstrip("\n")
    if l.startswith(">"):
        if names:
            seqs.append("".join(cur))
        n = l[1:].split()[0]
        names.append(n[3:] if n.startswith("_R_") else n)
        cur = []
    else:
        cur.append(l)
seqs.append("".join(cur))
cuts = [float(x) for x in sys.argv[3:]] or [0.90, 0.80, 0.70]
L = [sum(1 for c in s if c not in "-.") for s in seqs]
n = len(seqs)
ident = {}
for i in range(n):
    for j in range(i + 1, n):
        both = [(a.upper(), b.upper()) for a, b in zip(seqs[i], seqs[j]) if a not in "-." and b not in "-."]
        if len(both) >= 60 and len(both) >= 0.6 * min(L[i], L[j]):
            ident[(i, j)] = sum(a == b for a, b in both) / float(len(both))


def groups(cut):
    par = list(range(n))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for (i, j), v in ident.items():
        if v >= cut:
            par[f(i)] = f(j)
    g = {}
    for i in range(n):
        g.setdefault(f(i), []).append(i)
    return sorted(g.values(), key=lambda m: (-len(m), names[m[0]]))


res = {}
with open(sys.argv[2], "w") as out:
    out.write("name\t" + "\t".join("group_%.2f" % c for c in cuts) + "\n")
    for c in cuts:
        for k, m in enumerate(groups(c), 1):
            for i in m:
                res.setdefault(i, []).append("G%d" % k if len(m) > 1 else "-")
        print("== identity >= %.2f: %d groups with >1 member" % (c, sum(1 for m in groups(c) if len(m) > 1)))
        for k, m in enumerate(groups(c), 1):
            if len(m) > 1:
                print("  G%d (%d): %s" % (k, len(m), " ".join(sorted(names[i] for i in m))))
    for i in range(n):
        out.write(names[i] + "\t" + "\t".join(res[i]) + "\n")
