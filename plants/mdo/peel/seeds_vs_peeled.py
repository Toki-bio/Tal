#!/usr/bin/env python3
"""Compare every seed of the published plate with the owner's seven peeled consensuses.

Writes seeds_vs_peeled.tsv (best and second consensus per seed) and prints how the seeds fall onto the
consensuses at several identity cuts, next to the member counts in the consensus names. The group members
themselves were not given, so any membership read from this table is a reconstruction, not the owner's call.
Local alignment: match 5, mismatch -4, gap open -10, extend -1, both strands. Identity = identical bases /
alignment columns; coverage = aligned seed bases / seed length.
"""
from pathlib import Path

from Bio.Align import PairwiseAligner

H = Path(__file__).resolve().parent


def fasta(p):
    d, n = {}, None
    for line in p.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            n = line[1:].split()[0]
            d[n] = ""
        elif n:
            d[n] += line.strip()
    return d


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


cons = {k: v.upper() for k, v in fasta(H / "peel_round1_owner_consensuses.fa").items()}
seeds = {k: v.replace("-", "").upper() for k, v in fasta(H.parent / "alignments" / "mdo_consensuses.aln.fa").items()}
al = PairwiseAligner(mode="local", match_score=5, mismatch_score=-4, open_gap_score=-10, extend_gap_score=-1)


def compare(q, c):
    s_plus, s_minus = al.score(q, c), al.score(rc(q), c)
    qq, strand = (q, "+") if s_plus >= s_minus else (rc(q), "-")
    a = al.align(qq, c)[0]
    same = cols = used = 0
    prev = None
    for (qs, qe), (cs, ce) in zip(*a.aligned):
        qs, qe, cs, ce = int(qs), int(qe), int(cs), int(ce)
        if prev:
            cols += max(qs - prev[0], cs - prev[1])
        same += sum(x == y for x, y in zip(qq[qs:qe], c[cs:ce]))
        cols += qe - qs
        used += qe - qs
        prev = (qe, ce)
    return max(s_plus, s_minus), same / cols if cols else 0.0, used / len(q), strand


rows = []
for name, q in seeds.items():
    r = sorted(((compare(q, c), k) for k, c in cons.items()), key=lambda x: -x[0][0])
    (s1, i1, v1, d1), k1 = r[0]
    (s2, i2, v2, _), k2 = r[1]
    rows.append((name, len(q), k1, s1, i1, v1, d1, k2, s2, i2, v2))
with open(H / "seeds_vs_peeled.tsv", "w", encoding="utf-8", newline="\n") as f:
    f.write("seed\tlength\tbest\tbest_score\tbest_identity\tbest_seed_coverage\tstrand\tsecond\tsecond_score\t"
            "second_identity\tsecond_seed_coverage\n")
    for r in rows:
        f.write("\t".join(f"{x:.3f}" if isinstance(x, float) and i in (4, 5, 9, 10) else str(x)
                          for i, x in enumerate(r)) + "\n")

order = list(cons)
named = {k: int(k.split("_")[1].replace("seqs", "")) for k in order}
print("seeds", len(rows), "sum of named counts", sum(named.values()))
print("consensus\tnamed\t" + "\t".join(f"id>={c}&cov>=0.5" for c in (0.0, 0.6, 0.7, 0.75, 0.8)))
for k in order:
    print(k, named[k], *[sum(r[2] == k and r[4] >= c and r[5] >= 0.5 for r in rows) for c in (0.0, 0.6, 0.7, 0.75, 0.8)], sep="\t")
