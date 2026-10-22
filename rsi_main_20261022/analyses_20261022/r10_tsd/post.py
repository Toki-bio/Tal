"""Do copies of other subfamilies share what follows r10's end? Query = 100 bp right after r10's own match in each r10 copy
(control: 100 bp of the same copies' 5' flank, 60-160 bp before the element). Library = up to 300 copies of every other
subfamily (assigned), each read 50 bp 5' / 400 bp 3' in its own orientation, own consensus located by ssearch36. Hits are
placed relative to the library copy's own consensus end: inside its element, or past its end."""
import subprocess, random, collections, sys, os
RUN = sys.argv[1]; T = 16
rnd = random.Random(1)
cons = {}; n = None
for l in open(RUN + "/consensuses.clean.fa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:].split()[0]; cons[n] = ""
    else: cons[n] += l.upper()
fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(RUN + "/genome.clean.fa.fai")}
by = collections.defaultdict(list)
for l in list(open(RUN + "/results/assignment_full.tsv"))[1:]:
    f = l.split("\t")
    if f[4] == "assigned": by[f[1]].append(f[0])
def getpad(locs, up, dn, tag):
    bed = []
    for i, c in enumerate(locs):
        ch, rest = c.rsplit(":", 1); co, st = rest[:-3], rest[-2]; s, e = map(int, co.split("-"))
        a, b = (s - up, e + dn) if st == "+" else (s - dn, e + up)
        bed.append("%s\t%d\t%d\t%s%d\t0\t%s" % (ch, max(0, a), min(fai[ch], b), tag, i, st))
    open(tag + ".bed", "w").write("\n".join(bed) + "\n")
    out = subprocess.run(["bedtools", "getfasta", "-s", "-nameOnly", "-tab", "-fi", RUN + "/genome.clean.fa", "-bed", tag + ".bed"], capture_output=True, text=True, check=True).stdout
    return {l.split("\t")[0].split("(")[0]: l.split("\t")[1].upper() for l in out.splitlines()}
def ss(q, lib, b):
    out = subprocess.run(["ssearch36", "-m", "8", "-3", "-E", "1e-3", "-b", str(b), "-d", "0", "-T", str(T), q, lib], capture_output=True, text=True)
    assert out.returncode == 0 and "***" not in out.stderr, out.stderr[:300]
    return [l.split("\t") for l in out.stdout.splitlines() if l and not l.startswith("#")]
def ownmatch(name, pad):
    open("own_q.fa", "w").write(">%s\n%s\n" % (name, cons[name]))
    open("own_l.fa", "w").write("".join(">%s\n%s\n" % kv for kv in pad.items()))
    o = {}
    for f in ss("own_q.fa", "own_l.fa", len(pad) + 10):
        if int(f[9]) > int(f[8]) and (f[1] not in o or float(f[11]) > o[f[1]][2]): o[f[1]] = (int(f[8]), int(f[9]), float(f[11]), int(f[6]), int(f[7]))
    return o
# queries
r10 = getpad(by["r10_19seqs"], 200, 300, "q")
o10 = ownmatch("r10_19seqs", r10)
with open("query.fa", "w") as fh, open("ctrl.fa", "w") as fc:
    nq = 0
    for k, (s, e, _, qs, qe) in o10.items():
        if qe < 130: continue
        t = r10[k][e:e + 100]
        if len(t) == 100: fh.write(">%s\n%s\n" % (k, t)); nq += 1
        c = r10[k][max(0, s - 161):s - 61]
        if len(c) == 100: fc.write(">%s\n%s\n" % (k, c))
print("queries (100 bp past r10 end):", nq)
# library: other subfamilies
lib = {}; info = {}
for sf, locs in sorted(by.items()):
    if sf == "r10_19seqs": continue
    smp = locs if len(locs) <= 300 else rnd.sample(locs, 300)
    pad = getpad(smp, 50, 400, "L_" + sf.replace("-", "_") + "_")
    o = ownmatch(sf, pad)
    for k, (s, e, _, qs, qe) in o.items():
        lib[k] = pad[k]; info[k] = (sf, s, e, qe, len(cons[sf]))
open("lib.fa", "w").write("".join(">%s\n%s\n" % kv for kv in lib.items()))
print("library copies:", len(lib), "from", len({v[0] for v in info.values()}), "subfamilies")
for q, lab in (("query.fa", "r10 3' continuation"), ("ctrl.fa", "control: r10 5' genomic flank")):
    hits = ss(q, "lib.fa", 50)
    seen = set(); c = collections.Counter(); pos = collections.defaultdict(list)
    for f in hits:
        if int(f[3]) < 40: continue
        key = (f[0], f[1])
        if key in seen: continue
        seen.add(key)
        sf, s, e, qe, L = info[f[1]]
        ls, le = int(f[8]), int(f[9])
        where = "past its consensus end" if ls > e + 5 else ("inside its element" if le <= e + 5 else "across its end")
        c[(sf, where)] += 1
        if where != "inside its element": pos[sf].append(ls - e)
    qn = sum(1 for l in open(q) if l.startswith(">"))
    print("\n== %s: %d queries; query-copy pairs with a hit >= 40 bp, E <= 1e-3" % (lab, qn))
    tot = collections.Counter()
    for (sf, w), v in c.items(): tot[sf] += v
    for sf, v in tot.most_common(12):
        d = sorted(pos[sf])
        print("  %-18s %5d   inside element %4d | across end %4d | past end %4d%s" % (sf, v, c[(sf, "inside its element")], c[(sf, "across its end")],
              c[(sf, "past its consensus end")], ("   (past/across: median %d bp after its end)" % d[len(d)//2]) if d else ""))
