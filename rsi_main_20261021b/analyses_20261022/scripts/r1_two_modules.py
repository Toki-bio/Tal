"""All unit pieces aligned to r1 with its two modules as two column blocks (2026-10-25). r1 1-80 is an r10-type head and
r1 81-155 an r2/r4-type element (Tal rsi_main_20261022/README.md, "Which r1 module is r2's?"), so every piece is put
under the r1 module it is like: its best local score against r1 27-80 (module 1) and r1 81-133 (module 2) decides
(equal scores: module 2, header says "tie"; both below 18: a piece without a module, aligned only in the tail columns).
Module-1 pieces are aligned against r1 without module 2 (1-80 + 134-155), module-2 pieces against r1 without module 1
(1-26 + 76-155), each set with mafft --localpair. Each piece is then cut into 5' part | module | tail at r1's bases; the
5' parts of all pieces and the tails of all pieces are realigned together (mafft --localpair with r1's own 5' part /
tail), the two module blocks are kept side by side. The other r1 pieces are added with mafft --add --keeplength; every
other piece is checked to keep all its bases. Input order kept. Python 3.5 (KIT).
Usage: r1_two_modules.py rsi_all_unit_pieces.aln.fa OUT.aln.fa"""
import subprocess, sys
src, out = sys.argv[1], sys.argv[2]
seqs = []
for l in open(src):
    l = l.strip()
    if l.startswith(">"): seqs.append([l[1:], ""])
    elif seqs: seqs[-1][1] += l
seqs = [(h, s.replace("-", "").upper()) for h, s in seqs]
r1 = dict(seqs)["r1_9seqs@r1_9seqs"]
assert r1[75:80] == "GTGGG" and r1[80:89] == "CTCTCAACC" and r1[124:133] == "GATGGTGGG", r1
M1, M2 = r1[26:80], r1[80:133]

def sw(a, b, ma=2, mi=-3, go=-5, ge=-2):
    """best local score, affine gaps"""
    n = len(b); H = [0] * (n + 1); E = [-10**9] * (n + 1); best = 0
    for i in range(1, len(a) + 1):
        prev, F, Hn = 0, -10**9, [0] * (n + 1)
        for j in range(1, n + 1):
            E[j] = max(E[j] + ge, H[j] + go + ge)
            F = max(F + ge, Hn[j - 1] + go + ge)
            h = max(0, H[j - 1] + (ma if a[i - 1] == b[j - 1] else mi), E[j], F)
            Hn[j] = h; best = max(best, h)
        H = Hn
    return best

r1rows = [(h, s) for h, s in seqs if h.startswith("r1_9seqs@") and h != "r1_9seqs@r1_9seqs"]
others = [(h, s) for h, s in seqs if not h.startswith("r1_9seqs@")]
side = {}; tie = set()   # equal scores (r6: its module is half one type, half the other) go under module 2
for h, s in others:
    if not s: side[h] = "empty"; continue
    a, b = sw(s, M1), sw(s, M2)
    side[h] = "nomod" if max(a, b) < 18 else "m1" if a > b else "m2"   # < 18: a tail piece without a module (scores 10-15)
    if a == b and side[h] != "nomod": tie.add(h)
    print("%-40s m1 %3d  m2 %3d  -> %s" % (h, a, b, side[h]))

run = lambda args: subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True).stdout
def mafft(ref, rows, tag):
    with open(tag + ".in.fa", "w") as fh:
        fh.write(">REF\n%s\n" % ref)
        for i, (h, s) in enumerate(rows): fh.write(">%s%d\n%s\n" % (tag, i, s))
    A, n = {}, None
    for l in run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "8", "--quiet", tag + ".in.fa"]).splitlines():
        if l.startswith(">"): n = l[1:]; A[n] = ""
        else: A[n] += l.strip().upper()
    return A

def cut(A, starts):
    """split every row of A at the columns of r1 bases starts (indexes in r1), gap columns before a base go with it"""
    ref = A["REF"]; pos = []; t = 0
    for c, ch in enumerate(ref):
        if ch != "-": pos.append(c)
    # reference base number t in A -> its r1 index is refidx[t]
    bounds = [0]
    for st in starts:
        t = refidx[A is AA["b"]].index(st)
        c = pos[t]
        while c > 0 and ref[c - 1] == "-" and c - 1 >= bounds[-1]: c -= 1   # gap columns just before the base join it
        bounds.append(c)
    bounds.append(len(ref))
    return {n: [s[bounds[k]:bounds[k + 1]] for k in range(len(bounds) - 1)] for n, s in A.items()}

def squeeze(rows):
    """drop columns that are gaps in every row"""
    if not rows: return rows
    keep = [c for c in range(len(rows[0])) if any(r[c] != "-" for r in rows)]
    return ["".join(r[c] for c in keep) for r in rows]

p1 = list(range(0, 80)) + list(range(133, len(r1)))          # r1 without module 2
p2 = list(range(0, 26)) + list(range(75, len(r1)))           # r1 without module 1
refidx = {False: p1, True: p2}
g1 = [(h, s) for h, s in others if side[h] == "m1"]
g2 = [(h, s) for h, s in others if side[h] == "m2"]
AA = {"a": mafft("".join(r1[i] for i in p1), g1, "a"), "b": mafft("".join(r1[i] for i in p2), g2, "b")}
# segments: 5' part | module | tail. Module 1 rows: module = from r1 base 27 to r1 base 80, tail from r1 134 on.
# Module 2 rows: module = from r1 base 76 (the module's first bases, GCGGG in r2, lie over r1 76-80) to 133, tail 134 on.
C = {"a": cut(AA["a"], [26, 133]), "b": cut(AA["b"], [75, 133])}
key = {}
for t, g in (("a", g1), ("b", g2)):
    for i, (h, _) in enumerate(g): key[h] = (t, t + str(i))
seg5 = [(h, C[key[h][0]][key[h][1]][0].replace("-", "")) for h, _ in others if h in key]
segT = [(h, C[key[h][0]][key[h][1]][2].replace("-", "") if h in key else s) for h, s in others if h in key or side[h] == "nomod"]
def joint(ref, items, tag):
    """mafft --localpair of ref + the non-empty pieces; returns {name: aligned row} with REF"""
    ne = [(h, s) for h, s in items if s]
    A = mafft(ref, ne, tag)
    out = {"REF": A["REF"]}
    for i, (h, _) in enumerate(ne): out[h] = A[tag + str(i)]
    w = len(A["REF"])
    for h, s in items:
        if not s: out[h] = "-" * w
    return out
J5, JT = joint(r1[:26], seg5, "p"), joint(r1[133:], segT, "t")
m1rows = squeeze([C["a"]["REF"][1]] + [C["a"][key[h][1]][1] for h, _ in g1])
m2rows = squeeze([C["b"]["REF"][1]] + [C["b"][key[h][1]][1] for h, _ in g2])
# r1's row shows r1 76-80 once, in module 1: blank its first 5 bases in the module 2 block
r1m2, k = list(m2rows[0]), 0
for c, ch in enumerate(r1m2):
    if ch != "-" and k < 5: r1m2[c] = "-"; k += 1
m2rows[0] = "".join(r1m2)
assert m1rows[0].replace("-", "") == r1[26:80] and m2rows[0].replace("-", "") == r1[80:133]
w1, w2 = len(m1rows[0]), len(m2rows[0])
M1r = {h: m1rows[i + 1] for i, (h, _) in enumerate(g1)}
M2r = {h: m2rows[i + 1] for i, (h, _) in enumerate(g2)}
R1 = J5["REF"] + m1rows[0] + m2rows[0] + JT["REF"]
rows = {}
for h, s in others:
    if side[h] == "nomod": rows[h] = "-" * (len(J5["REF"]) + w1 + w2) + JT[h]; continue
    if h not in key: rows[h] = "-" * len(R1); continue
    rows[h] = J5[h] + M1r.get(h, "-" * w1) + M2r.get(h, "-" * w2) + JT[h]
with open("seed.fa", "w") as fh:
    fh.write(">R1\n%s\n" % R1)
    for i, (h, s) in enumerate(others): fh.write(">o%d\n%s\n" % (i, rows[h]))
with open("add.fa", "w") as fh:
    for i, (h, s) in enumerate(r1rows): fh.write(">x%d\n%s\n" % (i, s))
B, n = {}, None
for l in run(["mafft", "--add", "add.fa", "--keeplength", "--thread", "8", "--quiet", "seed.fa"]).splitlines():
    if l.startswith(">"): n = l[1:]; B[n] = ""
    else: B[n] += l.strip().upper()
for i, (h, s) in enumerate(others):                 # the merge must not change any piece
    assert B["o%d" % i].replace("-", "") == s, h
lab = {"m1": " [module 1, r10 type]", "m2": " [module 2, r2/r4 type]", "nomod": " [no module]", "empty": ""}
with open(out, "w") as fh:
    fh.write(">r1_9seqs@r1_9seqs (module 1 = r1 27-80, module 2 = r1 81-133)\n%s\n" % B["R1"])
    for i, (h, s) in enumerate(r1rows):
        kept = len(B["x%d" % i].replace("-", ""))
        fh.write(">%s%s\n%s\n" % (h, "" if kept == len(s) else " (%d of %d bases kept)" % (kept, len(s)), B["x%d" % i]))
    for i, (h, s) in enumerate(others): fh.write(">%s%s\n%s\n" % (h, lab[side[h]], B["o%d" % i]))
# track in r1's own coordinates (row positions)
open(out.replace(".aln.fa", ".bed"), "w").write('track name="r1 parts"\n'
    "r1_9seqs@r1_9seqs\t0\t26\t5' part (A box side)\t0\t+\t0\t0\t60,110,200\n"
    "r1_9seqs@r1_9seqs\t26\t80\tmodule 1 (r10 type)\t0\t+\t0\t0\t220,60,60\n"
    "r1_9seqs@r1_9seqs\t80\t133\tmodule 2 (r2/r4 type)\t0\t+\t0\t0\t230,140,40\n")
print("rows", 1 + len(r1rows) + len(others), "width", len(B["R1"]), "m1", len(g1), "m2", len(g2))
