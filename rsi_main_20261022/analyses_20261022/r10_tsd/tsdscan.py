"""Where do r10 copies end? Direct repeat of the 5' flank searched at every position up to 1 kb past r10's end.
Matching as SINEderella tools/length_variants.has_tsd (5' copy ends <= 4 bp before the element), length 10-20, <= 10 %
mismatches; null = 5' side of copy i with the 3' side of copy i+1."""
import subprocess, sys, collections, random
RUN, CONS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
DN, UP = 1000, 60
rows = [l.rstrip("\n").split("\t") for l in open(RUN + "/results/assignment_full.tsv")][1:]
cop = [r[0] for r in rows if r[1] == "r10_19seqs" and r[4] == "assigned"]
fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(RUN + "/genome.clean.fa.fai")}
bed = []
for i, c in enumerate(cop):
    ch, rest = c.rsplit(":", 1); co, st = rest[:-3], rest[-2]
    s, e = map(int, co.split("-"))
    a, b = (s - UP, e + DN) if st == "+" else (s - DN, e + UP)
    bed.append("%s\t%d\t%d\tc%d\t0\t%s" % (ch, max(0, a), min(fai[ch], b), i, st))
open(OUT + ".bed", "w").write("\n".join(bed) + "\n")
tab = subprocess.run(["bedtools", "getfasta", "-s", "-nameOnly", "-tab", "-fi", RUN + "/genome.clean.fa", "-bed", OUT + ".bed"],
                     capture_output=True, text=True, check=True).stdout
pad = {}
for l in tab.splitlines():
    n, s = l.split("\t"); pad[n.split("(")[0]] = s.upper()
with open(OUT + ".lib.fa", "w") as fh:
    for n, s in pad.items(): fh.write(">%s\n%s\n" % (n, s))
q = subprocess.run(["bash", "-c", "awk '/^>/{p=($1==\">r10_19seqs\")} p' " + CONS], capture_output=True, text=True).stdout
open(OUT + ".q.fa", "w").write(q)
m8 = subprocess.run(["ssearch36", "-m", "8", "-3", "-E", "1e-5", "-b", str(len(pad) + 10), "-d", "0", "-T", "16", OUT + ".q.fa", OUT + ".lib.fa"],
                    capture_output=True, text=True, check=True).stdout
own = {}
for l in m8.splitlines():
    f = l.split("\t")
    if f[1] not in own or float(f[11]) > own[f[1]][2]: own[f[1]] = (int(f[8]), int(f[9]), float(f[11]), int(f[6]), int(f[7]))
def tsds(up, dn):
    """all (d3, L) with a TSD: 5' copy within 4 bp before the element, 10-20 bp, <= 10 % mismatches; longest L per d3"""
    hit = {}
    for L in range(20, 9, -1):
        mm = int(L * 0.1)
        for d5 in range(5):
            a = up[len(up) - L - d5:len(up) - d5]
            for d3 in range(0, len(dn) - L + 1):
                if d3 in hit: continue
                if sum(x != y for x, y in zip(a, dn[d3:d3 + L])) <= mm: hit[d3] = L
    return hit
cs = [(n, pad[n][:o[0] - 1], pad[n][o[1]:], o) for n, o in own.items() if o[3] <= 15 and o[4] >= 130]   # 5'- and 3'-complete r10 matches
print("r10 copies", len(cop), "own match found", len(own), "complete r10 match", len(cs))
real = collections.Counter(); shuf = collections.Counter(); first = []
for i, (n, up, dn, o) in enumerate(cs):
    h = tsds(up, dn); hs = tsds(up, cs[(i + 1) % len(cs)][2])
    for d in h: real[d // 25] += 1
    for d in hs: shuf[d // 25] += 1
    first.append((n, min(h) if h else None, h.get(min(h)) if h else 0, len(dn)))
print("bin(bp past r10 end)\treal\tshuffled")
for b in range(0, DN // 25):
    if real[b] or shuf[b]: print("%d-%d\t%d\t%d%s" % (b * 25, b * 25 + 24, real[b], shuf[b], "  <==" if real[b] >= 3 * max(1, shuf[b]) and real[b] >= 5 else ""))
with open(OUT + ".first.tsv", "w") as fh:
    for r in first: fh.write("%s\t%s\t%s\t%s\n" % r)
