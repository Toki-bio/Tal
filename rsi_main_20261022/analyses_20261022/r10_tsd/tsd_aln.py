"""r10 copies with a TSD up to 250 bp past r10's end: plate with the TSDs in shared columns + a BED track (alignment columns).
Rows: left flank (lower) | 5' TSD (UPPER) | bases between it and the element (lower) | element = r10 match + everything to
the 3' TSD (UPPER, MAFFT) | 3' TSD (UPPER) | right flank (lower). Copies with verdict not_SINE are left out."""
import subprocess, sys, random
RUN, OUT = sys.argv[1], sys.argv[2]
MAXD, LF, RF, TW, XW = 250, 40, 60, 20, 4
pad = {}; n = None
for l in open("r10_22.lib.fa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:]; pad[n] = ""
    else: pad[n] += l
bed = {l.split("\t")[3]: l.split("\t") for l in open("r10_22.bed").read().split("\n") if l}
ver = {}
for l in open(RUN + "/step2/step2_output/loci_verdict.tsv"):
    f = l.rstrip("\n").split("\t"); ver[f[0]] = f[3]
rows = [l.rstrip("\n").split("\t") for l in open(RUN + "/results/assignment_full.tsv")][1:]
cop = [r[0] for r in rows if r[1] == "r10_19seqs" and r[4] == "assigned"]
m8 = subprocess.run(["ssearch36", "-m", "8", "-3", "-E", "1e-5", "-b", "400", "-d", "0", "-T", "16", "r10_22.q.fa", "r10_22.lib.fa"], capture_output=True, text=True).stdout
own = {}
for l in m8.splitlines():
    f = l.split("\t")
    if f[1] not in own or float(f[11]) > own[f[1]][2]: own[f[1]] = (int(f[8]), int(f[9]), float(f[11]), int(f[6]), int(f[7]))
def nearest(up, dn):
    for d3 in range(0, min(MAXD, len(dn) - 10) + 1):
        for L in range(20, 9, -1):
            if d3 + L > len(dn): continue
            for d5 in range(5):
                a = up[len(up) - L - d5:len(up) - d5]; b = dn[d3:d3 + L]
                mm = sum(x != y for x, y in zip(a, b))
                if mm <= int(L * 0.1): return d3, L, d5, mm
    return None
cs = [k for k, o in own.items() if o[3] <= 15 and o[4] >= 130]
real = []; shuf = 0
for i, k in enumerate(cs):
    o = own[k]; up, dn = pad[k][:o[0] - 1], pad[k][o[1]:]
    t = nearest(up, dn)
    if nearest(up, pad[cs[(i + 1) % len(cs)]][own[cs[(i + 1) % len(cs)]][1]:]): shuf += 1
    if t: real.append((k, t))
print("complete r10 matches %d; TSD within %d bp: real %d, shuffled pairs %d" % (len(cs), MAXD, len(real), shuf))
segs = []
for k, (d3, L, d5, mm) in sorted(real, key=lambda x: x[1][0]):
    loc = cop[int(k[1:])]
    v = ver.get(loc, "?")
    if v == "not_SINE": continue
    o = own[k]; up, dn = pad[k][:o[0] - 1], pad[k][o[1]:]
    E = pad[k][o[0] - 1:o[1] + d3]
    name = "%s%s TSD %d bp, %d mismatch, starts %d bp past r10 end" % (loc, " [grey]" if v == "grey" else "", L, mm, d3)
    segs.append((name, up[:len(up) - L - d5][-LF:].lower(), up[len(up) - L - d5:len(up) - d5].upper(), up[len(up) - d5:].lower(),
                 E.upper(), dn[d3:d3 + L].upper(), dn[d3 + L:d3 + L + RF].lower()))
cons = open("r10_22.q.fa").read().split("\n", 1)[1].replace("\n", "").upper()
with open(OUT + ".E.fa", "w") as fh:
    fh.write(">cons\n%s\n" % cons)
    for i, s in enumerate(segs): fh.write(">e%d\n%s\n" % (i, s[4]))
al = subprocess.run(["mafft", "--localpair", "--maxiterate", "100", "--thread", "8", "--quiet", OUT + ".E.fa"], capture_output=True, text=True, check=True).stdout
A = {}; n = None
for l in al.splitlines():
    if l.startswith(">"): n = l[1:]; A[n] = ""
    else: A[n] += l.strip().upper()
W = len(A["cons"])
# r10 consensus columns inside the element block: the columns where the consensus row has a base
cc = [i for i, c in enumerate(A["cons"]) if c != "-"]
with open(OUT + ".aln.fa", "w") as fh:
    fh.write(">r10_19seqs consensus (146 bp)\n%s\n" % ("-" * (LF + TW + XW) + A["cons"] + "-" * (TW + RF)))
    for i, s in enumerate(segs):
        row = s[1].rjust(LF, "-") + s[2].rjust(TW, "-") + s[3].rjust(XW, "-") + A["e%d" % i] + s[5].ljust(TW, "-") + s[6].ljust(RF, "-")
        fh.write(">%s\n%s\n" % (s[0], row))
o5 = LF; e0 = LF + TW + XW
with open(OUT + ".bed", "w") as fh:
    fh.write('track name="TSD"\n')
    fh.write("alignment\t%d\t%d\t5' TSD (10-20 bp, right-justified)\t0\t+\t%d\t%d\t220,60,60\n" % (o5, o5 + TW, o5, o5))
    fh.write("alignment\t%d\t%d\tr10 consensus\t0\t+\t%d\t%d\t60,110,200\n" % (e0 + cc[0], e0 + cc[-1] + 1, e0, e0))
    if e0 + cc[-1] + 1 < e0 + W:
        fh.write("alignment\t%d\t%d\tpast r10 end, up to the 3' TSD\t0\t+\t%d\t%d\t120,160,120\n" % (e0 + cc[-1] + 1, e0 + W, e0, e0))
    fh.write("alignment\t%d\t%d\t3' TSD (left-justified)\t0\t+\t%d\t%d\t220,60,60\n" % (e0 + W, e0 + W + TW, e0, e0))
print("rows %d (not_SINE left out), width %d" % (len(segs), LF + TW + XW + W + TW + RF))
