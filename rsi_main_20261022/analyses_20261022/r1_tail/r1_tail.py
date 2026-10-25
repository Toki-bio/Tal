"""What follows the 3' end of r1_9seqs copies (rsi_main_20261022)? r1's consensus (155 bp) ends at ...CTGCAACTAACAA where
r4 (and r3) go on (r1 81-155 = r4 30-106). Assigned r1 copies, padded 60 bp 5' and 1 kb 3' (as tsdscan.py for r10):
1. own r1 match; complete = r1 1-10 to >= 145;
2. the 300 bp past r1's end searched (ssearch36, E <= 1e-3) against every bank consensus: which consensus, at which
   positions, does the copy go on into;
3. nearest TSD (5' flank ending <= 4 bp before the element, repeated 10-20 bp, <= 10 % mismatches) up to 300 bp past
   r1's end, null = 5' side of copy i against 3' side of copy i+1;
4. plate: r1, r4, r3 consensuses + every complete copy from r1's start to its 3' TSD (or 250 bp past r1's end), MAFFT."""
import subprocess, sys
RUN, BANK, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
DN, UP, MAXD = 1000, 60, 300
run = lambda a: subprocess.run(a, capture_output=True, text=True, check=True).stdout
rows = [l.rstrip("\n").split("\t") for l in open(RUN + "/results/assignment_full.tsv")][1:]
cop = [r[0] for r in rows if r[1] == "r1_9seqs" and r[4] == "assigned"]
ver = {}
for l in open(RUN + "/step2/step2_output/loci_verdict.tsv"):
    f = l.rstrip("\n").split("\t"); ver[f[0]] = f[3]
fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(RUN + "/genome.clean.fa.fai")}
bed = []
for i, c in enumerate(cop):
    ch, rest = c.rsplit(":", 1); co, st = rest[:-3], rest[-2]
    s, e = map(int, co.split("-"))
    a, b = (s - UP, e + DN) if st == "+" else (s - DN, e + UP)
    bed.append("%s\t%d\t%d\tc%d\t0\t%s" % (ch, max(0, a), min(fai[ch], b), i, st))
open(OUT + ".bed", "w").write("\n".join(bed) + "\n")
pad = {}
for l in run(["bedtools", "getfasta", "-s", "-nameOnly", "-tab", "-fi", RUN + "/genome.clean.fa", "-bed", OUT + ".bed"]).splitlines():
    n, s = l.split("\t"); pad[n.split("(")[0]] = s.upper()
with open(OUT + ".lib.fa", "w") as fh:
    for n, s in pad.items(): fh.write(">%s\n%s\n" % (n, s))
cons, n = {}, None
for l in open(BANK):
    l = l.strip()
    if l.startswith(">"): n = l[1:].split()[0]; cons[n] = ""
    else: cons[n] += l.upper()
open(OUT + ".q.fa", "w").write(">r1_9seqs\n%s\n" % cons["r1_9seqs"])
own = {}
for l in run(["ssearch36", "-m", "8", "-3", "-E", "1e-5", "-b", str(len(pad) + 10), "-d", "0", OUT + ".q.fa", OUT + ".lib.fa"]).splitlines():
    if l.startswith("#") or not l.strip(): continue
    f = l.split("\t")
    if f[1] not in own or float(f[11]) > own[f[1]][2]: own[f[1]] = (int(f[8]), int(f[9]), float(f[11]), int(f[6]), int(f[7]))
L1 = len(cons["r1_9seqs"])
cs = [k for k in sorted(own, key=lambda x: int(x[1:])) if own[k][3] <= 10 and own[k][4] >= L1 - 10]
print("r1 assigned copies %d, own match %d, complete %d" % (len(cop), len(own), len(cs)))
# 2. what the 300 bp past r1's end match
with open(OUT + ".post.fa", "w") as fh:
    for k in cs: fh.write(">%s\n%s\n" % (k, pad[k][own[k][1]:own[k][1] + MAXD]))
with open(OUT + ".bank.fa", "w") as fh:
    for c, s in cons.items(): fh.write(">%s\n%s\n" % (c, s))
post = {}
for l in run(["ssearch36", "-m", "8", "-3", "-E", "1e-3", "-b", "40", "-d", "0", OUT + ".post.fa", OUT + ".bank.fa"]).splitlines():
    if l.startswith("#") or not l.strip(): continue
    f = l.split("\t"); post.setdefault(f[0], []).append((f[1], int(f[6]), int(f[7]), int(f[8]), int(f[9]), float(f[2]), float(f[11])))
def tsd(up, dn):
    for d3 in range(0, min(MAXD, len(dn) - 10) + 1):
        for L in range(20, 9, -1):
            if d3 + L > len(dn): continue
            for d5 in range(5):
                a = up[len(up) - L - d5:len(up) - d5]; b = dn[d3:d3 + L]
                mm = sum(x != y for x, y in zip(a, b))
                if mm <= int(L * 0.1): return d3, L, d5, mm
    return None
shuf = 0; res = []
for i, k in enumerate(cs):
    o = own[k]; up, dn = pad[k][:o[0] - 1], pad[k][o[1]:]
    nk = cs[(i + 1) % len(cs)]
    if tsd(up, pad[nk][own[nk][1]:]): shuf += 1
    res.append((k, tsd(up, dn)))
print("TSD within %d bp past r1's end: real %d of %d, shuffled pairs %d" % (MAXD, sum(1 for _, t in res if t), len(cs), shuf))
print("copy\tlocus\tverdict\tr1_match\tTSD_start_past_r1_end\tTSD_len\tbest hits of the 300 bp past r1's end (consensus:cons_start-cons_end q_start-q_end identity bits)")
segs = []
for k, t in res:
    loc = cop[int(k[1:])]; o = own[k]
    hits = sorted(post.get(k, []), key=lambda h: -h[6])
    best = {}
    for h in hits:
        if h[0] not in best: best[h[0]] = h
    hs = "; ".join("%s:%d-%d q%d-%d %.0f%% %.0f" % (h[0], h[3], h[4], h[1], h[2], h[5], h[6]) for h in sorted(best.values(), key=lambda h: -h[6])[:4])
    print("%s\t%s\t%s\t%d-%d\t%s\t%s\t%s" % (k, loc, ver.get(loc, "?"), o[3], o[4], t[0] if t else "-", t[1] if t else "-", hs or "none"))
    end = o[1] + (t[0] if t else 250)
    segs.append(("%s %s%s %s" % (loc, ver.get(loc, "?"), "", ("TSD %d bp at +%d" % (t[1], t[0])) if t else "no TSD within 300 bp, 250 bp shown"),
                 pad[k][o[0] - 1:end] + (pad[k][end:end + (t[1] if t else 0)].lower())))
with open(OUT + ".plate.in.fa", "w") as fh:
    for c in ("r1_9seqs", "r4_32seqs", "r3_58seqs", "r1_r3_P1"): fh.write(">%s (consensus)\n%s\n" % (c, cons[c]))
    for nm, s in segs: fh.write(">%s\n%s\n" % (nm, s))
open(OUT + ".plate.aln.fa", "w").write(run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "8", "--quiet", OUT + ".plate.in.fa"]))
