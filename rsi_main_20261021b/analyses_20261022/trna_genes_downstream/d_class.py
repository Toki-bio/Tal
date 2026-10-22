"""d_class.py (in D/): classify each RefSeq Ile tRNA gene outside long SINE hits as genuine (Pol III terminator: a run of >= 4 T
within 40 bp downstream of the gene end) or SINE-like (no T-run, and an A-rich stretch of >= 8 with >= 7 A within 40 bp), and
summarize how far the SINE consensuses align past the tRNA 3' end for each class (extension.tsv from d_trna.sh)."""
import collections, re, statistics
L = {k: int(v) for k, v in (l.rstrip("\n").split("\t") for l in open("glen.tsv"))}
seq, n = {}, None
for l in open("genes_dn.fa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:]; seq[n] = ""
    else: seq[n] += l
cls = {}
for g, s in seq.items():
    dn = s[L[g]:L[g] + 40].upper()
    trun = re.search(r"T{4,}", dn)
    arich = any(dn[i:i + 8].count("A") >= 7 for i in range(0, max(1, len(dn) - 7)))
    cls[g] = "genuine (T-run %d bp after the end)" % (trun.start() + 1) if trun else ("SINE-like (A-rich tail, no T-run)" if arich else "unclear (no T-run, no A-tail)")
rows = [l.rstrip("\n").split("\t") for l in open("extension.tsv")][1:]
by = collections.defaultdict(lambda: collections.defaultdict(list))
best_per_gene = {}
for r in rows:
    sine, g, bits, ext, shuf = r[0], r[1], float(r[2]), int(r[9]), int(r[10])
    k = cls[g].split(" (")[0]
    by[k][sine].append(ext)
    if g not in best_per_gene or bits > best_per_gene[g][1]:
        best_per_gene[g] = (sine, bits, ext)
with open("gene_classes.tsv", "w") as fo:
    fo.write("gene\tclass\tbest SINE (bits, bp past the tRNA end)\tdownstream 40 bp\n")
    for g in sorted(seq):
        b = best_per_gene.get(g)
        fo.write("%s\t%s\t%s\t%s\n" % (g, cls[g], "%s (%.1f, %d)" % b if b else "-", seq[g][L[g]:L[g] + 40]))
cnt = collections.Counter(c.split(" (")[0] for c in cls.values())
print("classes:", dict(cnt))
for k in by:
    allv = [v for s in by[k] for v in by[k][s]]
    print("%s: %d gene x SINE pairs, median bp past the tRNA end %s, max %d" % (k, len(allv), statistics.median(allv), max(allv)))
    for s in ("r9_15seqs", "r7_133seqs", "r8_83seqs", "r5_27seqs", "r6_210seqs", "r10_19seqs"):
        if by[k].get(s):
            print("   %s median %s max %d (n=%d)" % (s, statistics.median(by[k][s]), max(by[k][s]), len(by[k][s])))
