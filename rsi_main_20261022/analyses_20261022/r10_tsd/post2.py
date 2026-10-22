import subprocess, collections
info = {}
# rebuild the P18 library copies' own matches: run ownmatch again quickly
cons = {}; n = None
RUN = "/home/toki/rhin/rsi_main_20261022/run_20261022_174009"
for l in open(RUN + "/consensuses.clean.fa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:].split()[0]; cons[n] = ""
    else: cons[n] += l.upper()
lib = {}; n = None
for l in open("lib.fa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:]; lib[n] = ""
    else: lib[n] += l
p18 = {k: v for k, v in lib.items() if k.startswith("L_r10_r8_P18_")}
open("p18q.fa", "w").write(">P18\n%s\n" % cons["r10_r8_P18"])
open("p18l.fa", "w").write("".join(">%s\n%s\n" % kv for kv in p18.items()))
def m8(q, l, b):
    return [x.split("\t") for x in subprocess.run(["ssearch36", "-m", "8", "-3", "-E", "1e-3", "-b", str(b), "-d", "0", "-T", "16", q, l], capture_output=True, text=True).stdout.splitlines() if x and not x.startswith("#")]
# map library-copy position -> P18 consensus position with a global-ish local alignment (ssearch -m 8 gives only ends): use
# the aligned segment ends of hits of the r10 continuation and of P18 itself, interpolate linearly inside the P18 match
own = {}
for f in m8("p18q.fa", "p18l.fa", 400):
    if f[1] not in own or float(f[11]) > float(own[f[1]][11]): own[f[1]] = f
pos = collections.Counter(); ends = collections.Counter(); n = 0
for f in m8("query.fa", "p18l.fa", 50):
    if int(f[3]) < 40 or f[1] not in own: continue
    o = own[f[1]]; cs, ce, ls, le = int(o[6]), int(o[7]), int(o[8]), int(o[9])
    if not (ls <= int(f[8]) and int(f[9]) <= le + 5): continue
    k = lambda x: cs + (x - ls) * (ce - cs) / float(max(1, le - ls))
    pos[int(k(int(f[8])) // 20) * 20] += 1; ends[int(k(int(f[9])) // 20) * 20] += 1; n += 1
print("r10 continuation hits inside P18 copies:", n)
print("start on P18 consensus (20-bp bins):", sorted(pos.items()))
print("end   on P18 consensus (20-bp bins):", sorted(ends.items()))
# P18 copies: length of their own element and coverage of the middle piece 146-291
L = collections.Counter(min(int(int(o[9]) - int(o[8])) // 50, 12) * 50 for o in own.values())
print("P18 copies: own-match length on the copy (50-bp bins):", sorted(L.items()))
