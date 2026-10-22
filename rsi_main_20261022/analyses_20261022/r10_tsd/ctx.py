"""r10 copies (element + everything to their 3' TSD, from tsd_aln) aligned with 40 r10_r8_P18 copies (SINE verdict, best
bitscore; own match + 40 bp) and the r10, P18 consensuses; BED marks P18's tandem middle piece (146-291) on the P18 row."""
import subprocess
RUN = "/home/toki/rhin/rsi_main_20261022/run_20261022_174009"
cons = {}; n = None
for l in open(RUN + "/consensuses.clean.fa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:].split()[0]; cons[n] = ""
    else: cons[n] += l.upper()
ver = {l.split("\t")[0]: l.split("\t")[3] for l in open(RUN + "/step2/step2_output/loci_verdict.tsv")}
rows = [l.rstrip("\n").split("\t") for l in open(RUN + "/results/assignment_full.tsv")][1:]
p = sorted([r for r in rows if r[1] == "r10_r8_P18" and r[4] == "assigned" and ver.get(r[0]) == "SINE"], key=lambda r: -int(r[2]))[:40]
fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(RUN + "/genome.clean.fa.fai")}
bed = []
for i, r in enumerate(p):
    ch, rest = r[0].rsplit(":", 1); co, st = rest[:-3], rest[-2]; s, e = map(int, co.split("-"))
    bed.append("%s\t%d\t%d\tp%d\t0\t%s" % (ch, max(0, s - 50), min(fai[ch], e + 50), i, st))
open("p.bed", "w").write("\n".join(bed) + "\n")
pad = {l.split("\t")[0].split("(")[0]: l.split("\t")[1].upper() for l in subprocess.run(["bedtools", "getfasta", "-s", "-nameOnly", "-tab", "-fi", RUN + "/genome.clean.fa", "-bed", "p.bed"], capture_output=True, text=True, check=True).stdout.splitlines()}
open("pq.fa", "w").write(">P18\n%s\n" % cons["r10_r8_P18"]); open("pl.fa", "w").write("".join(">%s\n%s\n" % kv for kv in pad.items()))
own = {}
for x in subprocess.run(["ssearch36", "-m", "8", "-3", "-E", "1e-5", "-b", "100", "-d", "0", "-T", "16", "pq.fa", "pl.fa"], capture_output=True, text=True).stdout.splitlines():
    f = x.split("\t")
    if len(f) > 11 and int(f[9]) > int(f[8]) and (f[1] not in own or float(f[11]) > own[f[1]][2]): own[f[1]] = (int(f[8]), int(f[9]), float(f[11]))
E = {}; n = None
for l in open("/home/toki/tmp/tsd_r10/rsi_r10_tsd.E.fa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:]; E[n] = ""
    else: E[n] += l
names = [l[1:].strip() for l in open("/home/toki/tmp/tsd_r10/rsi_r10_tsd.aln.fa") if l.startswith(">")][1:]
with open("ctx_in.fa", "w") as fh:
    fh.write(">r10_19seqs consensus (146 bp)\n%s\n>r10_r8_P18 consensus (426 bp)\n%s\n" % (cons["r10_19seqs"], cons["r10_r8_P18"]))
    for i, nm in enumerate(names):
        fh.write(">r10 copy %s, element to its 3' TSD\n%s\n" % (nm.split(" ")[0] + (" [grey]" if "[grey]" in nm else ""), E["e%d" % i]))
    for i, r in enumerate(p):
        k = "p%d" % i
        if k not in own: continue
        s, e, _ = own[k]
        fh.write(">P18 copy %s, bitscore %s\n%s\n" % (r[0], r[2], pad[k][s - 1:e]))
al = subprocess.run(["mafft", "--localpair", "--maxiterate", "100", "--thread", "16", "--quiet", "ctx_in.fa"], capture_output=True, text=True, check=True).stdout
open("rsi_r10_vs_P18.aln.fa", "w").write(al)
open("rsi_r10_vs_P18.bed", "w").write('track name="P18 parts"\n'
    "r10_r8_P18\t0\t145\tP18 1-145 (r10 part)\t0\t+\t0\t0\t60,110,200\n"
    "r10_r8_P18\t145\t291\tP18 146-291: tandem of ~10-bp GC-rich units\t0\t+\t0\t0\t220,140,40\n"
    "r10_r8_P18\t291\t426\tP18 292-426 (r8 module)\t0\t+\t0\t0\t60,160,90\n")
print("rows", al.count(">"), "width", len(al.split(">")[1].split("\n", 1)[1].replace("\n", "")))
