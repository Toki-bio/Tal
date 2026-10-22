"""c_heads.py TRNA_DIR OUT_DIR - task C: tRNA heads of the rsi consensuses with their best and close tRNA genes.

From the run's tools/trna_heads.sh output (queries.fa = every consensus + every composite stretch, ss.m8 = ssearch36 -3 -m 8
against TRNA_REF hg19 GtRNAdb, trna_heads.tsv with null_max_bits from 20 dinucleotide shuffles).
A query has a tRNA head when its best tRNA bits exceed the best shuffle (the report's rule). For it:
  best   = the best-scoring gene (ties: the first in the file);
  close  = every other gene with bits >= CLOSE x best bits whose sequence differs from the genes already kept (identical
           paralogues collapsed), at most MAXCLOSE, preferring other anticodons / amino acids first.
Rows per amino-acid group: each query's head (q_span start .. end + 30 bp), then the genes, named
  <gene>|best:<queries>|close:<queries>.
Writes OUT/trna_heads_<aa>.fa (unaligned) and OUT/trna_heads_table.tsv."""
import collections, os, re, sys

TD, OUT = sys.argv[1], sys.argv[2]
CLOSE, MAXCLOSE, PAST = 0.90, 3, 30


def rd(p):
    d, n = collections.OrderedDict(), None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = []
        elif n:
            d[n].append(l.strip())
    return collections.OrderedDict((k, "".join(v).upper()) for k, v in d.items())


q, ref = rd(os.path.join(TD, "queries.fa")), rd(os.path.join(TD, "ref.fa"))
iso = lambda g: re.search(r"tRNA-([A-Za-z]+)-([A-Z]{3})", g).groups()
head = {}
for l in open(os.path.join(TD, "trna_heads.tsv")):
    f = l.rstrip("\n").split("\t")
    if f[0] == "query":
        H = f; continue
    head[f[0]] = dict(zip(H, f))
hits = collections.defaultdict(list)
for l in open(os.path.join(TD, "ss.m8")):
    f = l.split("\t")
    if "tRNA" in f[1]:
        hits[f[0]].append((float(f[11]), f[1], int(f[6]), int(f[7]), float(f[2])))

groups = collections.defaultdict(list)       # aa -> [(query, a, b)]
genes = collections.OrderedDict()            # gene -> {"best": [], "close": []}
tab = []
for qn in q:
    r = head.get(qn)
    if not r or r["best_tRNA"] == "-":
        continue
    hs = sorted(hits[qn], key=lambda h: -h[0])
    if not hs or hs[0][0] <= float(r["null_max_bits"]):
        tab.append((qn, "no tRNA head (best %s bits <= shuffle %s)" % (hs[0][0] if hs else 0, r["null_max_bits"]), "", ""))
        continue
    b = hs[0]
    aa = iso(b[1])[0]
    kept_seqs = {ref[b[1]]}
    genes.setdefault(b[1], {"best": [], "close": []})["best"].append(qn)
    # close: other anticodon / amino acid first, then other genes of the same anticodon
    cand = [h for h in hs[1:] if h[0] >= CLOSE * b[0]]
    cand.sort(key=lambda h: (iso(h[1]) == iso(b[1]), -h[0]))
    close = []
    for h in cand:
        if ref[h[1]] in kept_seqs:
            continue
        kept_seqs.add(ref[h[1]]); close.append(h)
        genes.setdefault(h[1], {"best": [], "close": []})["close"].append(qn)
        if len(close) >= MAXCLOSE:
            break
    other_aa = next((h for h in hs if iso(h[1])[0] != aa), None)
    a0, b0 = b[2], b[3]
    groups[aa].append((qn, a0, min(len(q[qn]), b0 + PAST)))
    tab.append((qn, "%s %.1f bits (%.1f %%), q %d-%d" % (b[1].replace("Homo_sapiens_", ""), b[0], b[4], a0, b0),
                "; ".join("%s %.1f (%.0f %% of best)" % (h[1].replace("Homo_sapiens_", ""), h[0], 100 * h[0] / b[0]) for h in close) or "-",
                ("%s %.1f (%.0f %% of best)" % (other_aa[1].replace("Homo_sapiens_", ""), other_aa[0], 100 * other_aa[0] / b[0])) if other_aa else "-"))

os.makedirs(OUT, exist_ok=True)
for aa, qs in groups.items():
    with open(os.path.join(OUT, "trna_heads_%s.fa" % aa), "w") as fo:
        for qn, a, b in qs:
            fo.write(">%s_head:%d-%d\n%s\n" % (qn.replace("|", "~"), a, b, q[qn][a - 1:b]))
        for g, d in genes.items():
            if iso(g)[0] != aa and not any(x in [y[0] for y in qs] for x in d["best"] + d["close"]):
                continue
            if not any(x in [y[0] for y in qs] for x in d["best"] + d["close"]):
                continue
            nm = g.replace("Homo_sapiens_", "") + "_human"
            if d["best"]:
                nm += "|best:" + ",".join(x.split(":")[0] if ":" not in x else x.replace("|", "~") for x in d["best"][:6]) + ("..." if len(d["best"]) > 6 else "")
            if d["close"]:
                nm += "|close:" + ",".join(x.replace("|", "~") for x in d["close"][:6]) + ("..." if len(d["close"]) > 6 else "")
            fo.write(">%s\n%s\n" % (nm, ref[g]))
with open(os.path.join(OUT, "trna_heads_table.tsv"), "w") as fo:
    fo.write("query\tbest tRNA gene\tclose genes (>= 90 %% of best bits, distinct sequences)\tbest other amino acid\n")
    for t in tab:
        fo.write("\t".join(t) + "\n")
print("groups", {k: len(v) for k, v in groups.items()}, "genes", len(genes))
