#!/usr/bin/env python3
"""Add soft-call columns to a tribes summary: for each tribe, the Soft_Subfamily of its members
that are in results/unassigned.tsv (copies that failed 10/10 unanimity or the bitscore bar).

usage: add_soft.py TRIBE_DIR CODE RUN_ROOT      (rewrites <code>_tribes_summary.tsv in place)
Adds: soft_count, soft_dominant, soft_breakdown. Idempotent (existing soft_* columns replaced).
"""
import csv
import os
import sys
from collections import Counter, defaultdict

tdir, code, run = sys.argv[1:4]
summ = os.path.join(tdir, f"{code}_tribes_summary.tsv")
uc = os.path.join(tdir, "_work", "clusters.uc")

soft = {}
with open(os.path.join(run, "results", "unassigned.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        soft[r["SeqID"]] = r["Soft_Subfamily"]

members = defaultdict(list)
with open(uc) as fh:
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if c[0] in ("S", "H"):
            members[c[1]].append(c[8].split()[0])

with open(summ, newline="") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    fields = [f for f in rd.fieldnames if not f.startswith("soft_")]
    rows = list(rd)

for r in rows:
    cnt = Counter(soft[m] for m in members[r["cluster_id"]] if m in soft and soft[m])
    r["soft_count"] = sum(cnt.values())
    r["soft_dominant"] = cnt.most_common(1)[0][0] if cnt else ""
    r["soft_breakdown"] = ",".join(f"{k}:{v}" for k, v in cnt.most_common())

tmp = summ + ".tmp"
with open(tmp, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields + ["soft_count", "soft_dominant", "soft_breakdown"],
                       delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
os.replace(tmp, summ)
n = sum(1 for r in rows if int(r["soft_count"]))
print(f"{code}: {len(rows)} tribes, {n} with soft-called members")
