"""All unit pieces aligned with r1's first B-box module (r1 27-75) set apart: every other piece is aligned against r1 without
that module (r1 1-26 + 76-155), so their single module sits under r1's second one; the module is then put back as its own
column block (gaps in all other rows). The other r1 pieces are added with mafft --add --keeplength (bases they carry beyond
r1's columns are dropped; none expected, they are r1 copies). Input order kept."""
import subprocess, sys
src, out = sys.argv[1], sys.argv[2]
seqs = []; n = None
for l in open(src):
    l = l.strip()
    if l.startswith(">"): seqs.append([l[1:], ""])
    else: seqs[-1][1] += l
seqs = [(h, s.replace("-", "").upper()) for h, s in seqs]
r1 = dict(seqs)["r1_9seqs@r1_9seqs"]
assert r1[80:89] == "CTCTCAACC" and r1[75:80] == "GTGGG", r1
mod = r1[26:75]                     # r1 27-75
r1d = r1[:26] + r1[75:]
others = [(h, s) for h, s in seqs if not h.startswith("r1_9seqs@") and s]
r1rows = [(h, s) for h, s in seqs if h.startswith("r1_9seqs@") and h != "r1_9seqs@r1_9seqs"]
with open("in.fa", "w") as fh:
    fh.write(">R1D\n%s\n" % r1d)
    for i, (h, s) in enumerate(others): fh.write(">o%d\n%s\n" % (i, s))
run = lambda args: subprocess.run(args, capture_output=True, text=True, check=True).stdout
A = {}; order = []
for l in run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "8", "--quiet", "in.fa"]).splitlines():
    if l.startswith(">"): n = l[1:]; A[n] = ""; order.append(n)
    else: A[n] += l.strip().upper()
# column after r1d's 26th base
k = 0
for c, ch in enumerate(A["R1D"]):
    if ch != "-":
        k += 1
        if k == 26: cut = c + 1; break
W = len(mod)
full = {n: s[:cut] + ("-" * W if n != "R1D" else mod) + s[cut:] for n, s in A.items()}
with open("seed.fa", "w") as fh:
    fh.write(">R1\n%s\n" % full["R1D"])
    for i in range(len(others)): fh.write(">o%d\n%s\n" % (i, full["o%d" % i]))
with open("add.fa", "w") as fh:
    for i, (h, s) in enumerate(r1rows): fh.write(">a%d\n%s\n" % (i, s))
B = {}
for l in run(["mafft", "--add", "add.fa", "--keeplength", "--thread", "8", "--quiet", "seed.fa"]).splitlines():
    if l.startswith(">"): n = l[1:]; B[n] = ""
    else: B[n] += l.strip().upper()
with open(out, "w") as fh:
    fh.write(">r1_9seqs@r1_9seqs (first B-box module r1 27-75 set apart)\n%s\n" % B["R1"])
    for i, (h, s) in enumerate(r1rows):
        kept = len(B["a%d" % i].replace("-", ""))
        fh.write(">%s%s\n%s\n" % (h, "" if kept == len(s) else " (%d of %d bases kept)" % (kept, len(s)), B["a%d" % i]))
    for i, (h, s) in enumerate(others): fh.write(">%s\n%s\n" % (h, B["o%d" % i]))
open(out.replace(".aln.fa", ".bed"), "w").write('track name="r1 head parts"\n'
    "r1_9seqs@r1_9seqs\t0\t26\t5' part (A box side)\t0\t+\t0\t0\t60,110,200\n"
    "r1_9seqs@r1_9seqs\t26\t75\tB-box module 1 (r1 only)\t0\t+\t0\t0\t220,60,60\n"
    "r1_9seqs@r1_9seqs\t75\t133\tB-box module 2 (all units)\t0\t+\t0\t0\t230,140,40\n")
print("rows", 1 + len(r1rows) + len(others), "width", len(B["R1"]), "cut column", cut + 1, "module", mod)
