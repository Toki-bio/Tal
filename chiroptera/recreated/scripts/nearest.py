#!/usr/bin/env python3
"""nearest.py ALN MEMBERS.tsv OUT.tsv - for every rebuilt consensus, identity to each REF_ row in the
alignment (columns where both have a letter; >= 40 shared columns), the closest reference, and how much
of that reference it covers. A low-copy "Rhin-1" that is really another element shows up as a different
closest reference, or a low identity to every one."""
import sys
from collections import OrderedDict

names, seqs, cur = [], [], []
for l in open(sys.argv[1]):
    l = l.rstrip("\n")
    if l.startswith(">"):
        if names:
            seqs.append("".join(cur))
        names.append(l[1:].split()[0].replace("_R_", "", 1) if l[1:].startswith("_R_") else l[1:].split()[0])
        cur = []
    else:
        cur.append(l)
seqs.append("".join(cur))
refs = [(n, s) for n, s in zip(names, seqs) if n.startswith("REF_")]
meta = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(sys.argv[2]) if not l.startswith("name\t")}
out = open(sys.argv[3], "w")
out.write("name\tcode\tbat_family\tsine_family\tn\tscore\tlen\tclosest_ref\tidentity\tref_covered\tidentity_to_own_query\tsecond_ref\tsecond_identity\n")
for n, s in zip(names, seqs):
    if n.startswith("REF_") or n not in meta:
        continue
    res = []
    for rn, rs in refs:
        both = [(a.upper(), b.upper()) for a, b in zip(s, rs) if a not in "-." and b not in "-."]
        rlen = sum(1 for b in rs if b not in "-.")
        if len(both) >= 40:
            res.append((sum(a == b for a, b in both) / float(len(both)), rn[4:], len(both) / float(rlen)))
    res.sort(reverse=True)
    m = meta[n]
    own = m[3] if not m[3].startswith("r") else m[3]
    own_id = next(("%.2f" % r[0] for r in res if r[1] == own), "")
    best = res[0] if res else (0, "-", 0)
    sec = res[1] if len(res) > 1 else (0, "-", 0)
    out.write("\t".join([n, m[1], m[2], m[3], m[4], m[5], m[7], best[1], "%.2f" % best[0], "%.2f" % best[2],
                         own_id, sec[1], "%.2f" % sec[0]]) + "\n")
