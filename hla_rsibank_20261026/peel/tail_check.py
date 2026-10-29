"""Shortie tail check for the hla 30k SubFam plate peel (owner's rows 321-347 = G47, 526-598 = G28).
For every copy in those chunks: extend 50 bp 5' and 300 bp 3' in the genome (strand-aware), align the owner's
consensuses to it (ssearch36 -z 11), and report where on each consensus the copy ends and whether the genome copy
runs past the step-1 hit end.
Usage: tail_check.py <chunks.tsv> <groups.tsv> <genome.fa> <owner.fa> <out dir>"""
import collections, os, subprocess, sys

chunks, groups, genome, owner, out = sys.argv[1:6]
EXT5, EXT3 = 50, 300
os.makedirs(out, exist_ok=True)
grp = dict(l.rstrip("\n").split("\t") for l in open(groups))
size = dict((l.split("\t")[0], int(l.split("\t")[1])) for l in open(genome + ".fai"))
bed, meta = [], {}
for l in open(chunks):
    cid, ch = l.rstrip("\n").split("\t")[:2]           # column 3 is SubFam's orientation, not the genome strand
    if ch not in grp:
        continue
    c, rest = cid.rsplit(":", 1)
    st = rest[-2]                                       # genome strand from the copy name, e.g. ...(-)
    s, e = rest.split("(")[0].split("-")
    s, e = int(s), int(e)                      # bedtools getfasta style: 0-based start, end exclusive
    if st == "+":
        a, b = max(0, s - EXT5), min(size[c], e + EXT3); off5 = s - a
    else:
        a, b = max(0, s - EXT3), min(size[c], e + EXT5); off5 = b - e
    name = "%s|%s|%d|%d" % (cid, ch, off5, e - s)
    bed.append("%s\t%d\t%d\t%s\t0\t%s" % (c, a, b, name, st))
open(out + "/ext.bed", "w").write("\n".join(bed) + "\n")
subprocess.run("bedtools getfasta -s -nameOnly -fi %s -bed %s/ext.bed > %s/ext.fa" % (genome, out, out), shell=True, check=True)
subprocess.run("ssearch36 -z 11 -T 16 -m 8 -b 100000 -d 0 -E 1e-3 -3 %s %s/ext.fa > %s/hits.m8" % (owner, out, out),
               shell=True, check=True)
qlen, q = {}, None
for l in open(owner):
    if l.startswith(">"):
        q = l[1:].split()[0]; qlen[q] = 0
    else:
        qlen[q] += len(l.strip())
best = {}
for l in open(out + "/hits.m8"):
    f = l.split("\t")
    qi, si = f[0], f[1]
    if si.endswith("(+)") or si.endswith("(-)"):          # bedtools -s -nameOnly appends the strand
        si = si[:-3]
    bs = float(f[11])
    if (qi, si) not in best or bs > best[(qi, si)][0]:
        best[(qi, si)] = (bs, int(f[6]), int(f[7]), int(f[8]), int(f[9]), float(f[2]))
rows = []
for b in bed:
    name = b.split("\t")[3]
    cid, ch, off5, hl = name.split("|"); off5, hl = int(off5), int(hl)
    hit_end = off5 + hl                         # 1-based end of the step-1 hit inside the extended sequence
    r = [cid, ch, grp[ch], str(hl)]
    for qn in ("consensus_28seqs", "consensus_47seqs", "r5_r6_long_hla"):
        h = best.get((qn, name))
        if h is None:
            r += ["-", "-", "-", "-"]
        else:
            r += ["%d-%d" % (h[1], h[2]), "%.1f" % h[5], str(h[4] - hit_end), "%.0f" % h[0]]
    rows.append(r)
hdr = ["copy", "chunk", "group", "hit_len"]
for qn in ("c28", "c47", "long"):
    hdr += [qn + "_span", qn + "_pid", qn + "_end_minus_hit_end", qn + "_bits"]
with open(out + "/tail_check.tsv", "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
# summary per group: c28 end on consensus (tail = 103-134 on consensus_28seqs), past-hit extension, long continuation
summ = collections.defaultdict(collections.Counter)
for r in rows:
    g = r[2]; summ[g]["copies"] += 1
    sp = r[4]
    if sp == "-":
        summ[g]["c28_none"] += 1; continue
    qe = int(sp.split("-")[1])
    summ[g]["c28_end>=128 (tail full)" if qe >= 128 else ("c28_end104-127 (tail partial)" if qe >= 104 else "c28_end<104 (no tail)")] += 1
    if int(r[6]) > 10:
        summ[g]["genome copy runs >10 bp past the step-1 hit end"] += 1
    ls = r[12]
    if ls != "-" and int(ls.split("-")[1]) > 150:
        summ[g]["aligns past long_hla pos 150 (second monomer)"] += 1
with open(out + "/tail_summary.txt", "w") as fh:
    for g in sorted(summ):
        fh.write(g + "\t" + "; ".join("%s: %d" % kv for kv in sorted(summ[g].items())) + "\n")
print(open(out + "/tail_summary.txt").read())
