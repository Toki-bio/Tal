"""unit_pieces_copies.py RUN OUTDIR [N=50] - genome copies of every elementary unit, standalone and cut out of composites.

Groups (one per unit x source):
  <unit>@alone           : N random SINE-verdict (unique flanks) copies assigned to the unit's own consensus, whole copy
  <unit>@<composite>:a-b : N random SINE-verdict copies of the composite; each copy aligned to the composite consensus
                           (ssearch36 -3, BTOP), the stretch a-b projected onto the copy and cut there; kept when the
                           alignment covers >= 80 % of the stretch
Stretches: composite_units.tsv (composites; unassigned stretches >= 30 bp as unassigned@...) + EXTRA_CUTS below.
Out: OUTDIR/pieces/<group>.fa, OUTDIR/groups.tsv (group, unit, source, stretch, copies tried, pieces kept)."""
import csv
import gzip
import os
import random
import re
import subprocess
import sys

run, out = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 50
S2 = sorted(d for d in os.listdir(os.path.join(run, "step2")) if d.startswith("step2_output"))[-1]
S2 = os.path.join(run, "step2", S2)
# P18del is r10 head + the shortened second unit; cut it like P18 (r10 ends at 151)
EXTRA_CUTS = [("r10_r8_P18del", 1, 151, "r10_19seqs"), ("r10_r8_P18del", 152, 255, "P18_2nd")]
P18_2ND = {"r10_r8_P18": (290, 426)}                         # P18's second unit, labelled P18_2nd too


def fasta(p):
    s, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            n = l[1:].split()[0]; s[n] = []
        elif n:
            s[n].append(l)
    return {k: "".join(v).upper() for k, v in s.items()}


cons = fasta(os.path.join(run, "consensuses.clean.fa"))
rows = list(csv.DictReader(open(os.path.join(S2, "composite_units.tsv")), delimiter="\t"))
units = [r["consensus"] for r in rows if r["kind"] == "unit" and not r["consensus"].startswith("MEG")]
cuts = []                                                    # (composite, a, b, unit label)
for r in rows:
    if r["kind"] != "composite":
        continue
    a, b = int(r["start"]), int(r["end"])
    if b - a + 1 < 30:
        continue
    lab = r["unit"] if r["unit"] != "-" else "unassigned"
    if P18_2ND.get(r["consensus"]) == (a, b):
        lab = "P18_2nd"
    cuts.append((r["consensus"], a, b, lab))
cuts += EXTRA_CUTS

op = gzip.open if os.path.exists(os.path.join(S2, "loci_verdict.tsv.gz")) else open
vf = os.path.join(S2, "loci_verdict.tsv.gz") if op is gzip.open else os.path.join(S2, "loci_verdict.tsv")
members = {}
for r in csv.DictReader(op(vf, "rt"), delimiter="\t"):
    if r["assignment_status"] == "assigned" and r["verdict"] == "SINE":
        members.setdefault(r["subfamily"], []).append(r["copy"])
need = set()
rng = random.Random(42)
pick = {}
for c in set(units) | {c for c, _, _, _ in cuts}:
    m = sorted(members.get(c, []))
    pick[c] = rng.sample(m, min(N, len(m)))
    need.update(pick[c])
seqs = {}
n = None
for l in open(os.path.join(run, "genome.clean_step1", "extracted.fasta")):
    if l.startswith(">"):
        n = l[1:].strip().split()[0]
        n = n if n in need else None
        if n:
            seqs[n] = []
    elif n:
        seqs[n].append(l.strip())
seqs = {k: "".join(v).upper() for k, v in seqs.items()}
os.makedirs(os.path.join(out, "pieces"), exist_ok=True)
tmp = os.path.join(out, "tmp"); os.makedirs(tmp, exist_ok=True)


def btop_map(qs, ss, btop):
    """subject position -> query position (1-based) along a BTOP alignment"""
    m, q, s = {}, qs, ss
    for tok in re.findall(r"\d+|[A-Z*-]{2}", btop):
        if tok.isdigit():
            for _ in range(int(tok)):
                m[s] = q; q += 1; s += 1
        else:
            x, y = tok
            if x != "-" and y != "-":
                m[s] = q; q += 1; s += 1
            elif x == "-":                                  # gap in the copy
                m[s] = None; s += 1
            else:                                            # insertion in the copy
                q += 1
    return m


aligned = {}


def align(comp):
    """every picked copy of a composite against its consensus -> {copy: (qs, qe, ss, se, map)}"""
    if comp in aligned:
        return aligned[comp]
    qf, tf = os.path.join(tmp, "q.fa"), os.path.join(tmp, "t.fa")
    with open(qf, "w") as fo:
        for c in pick[comp]:
            fo.write(">%s\n%s\n" % (c, seqs[c]))
    with open(tf, "w") as fo:
        fo.write(">%s\n%s\n" % (comp, cons[comp]))
    res = subprocess.run(["ssearch36", "-3", "-m", "8CB", "-E", "1e-3", "-T", "8", "-d", "0", qf, tf],
                         capture_output=True, text=True, check=True).stdout
    best = {}
    for l in res.splitlines():
        if l.startswith("#") or not l.strip():
            continue
        f = l.split("\t")
        if f[0] not in best or float(f[11]) > best[f[0]][0]:
            best[f[0]] = (float(f[11]), int(f[6]), int(f[7]), int(f[8]), int(f[9]), f[12])
    aligned[comp] = {k: (v[1], v[2], v[3], v[4], btop_map(v[1], v[3], v[5])) for k, v in best.items()}
    return aligned[comp]


with open(os.path.join(out, "groups.tsv"), "w") as gt:
    gt.write("group\tunit\tsource\tstretch\tcopies\tpieces\n")
    for u in units:
        g = "%s@alone" % u
        with open(os.path.join(out, "pieces", g + ".fa"), "w") as fo:
            fo.write(">%s|consensus\n%s\n" % (u, cons[u]))
            for c in pick[u]:
                fo.write(">%s|%s\n%s\n" % (c, g, seqs[c]))
        gt.write("%s\t%s\talone\t1-%d\t%d\t%d\n" % (g, u, len(cons[u]), len(pick[u]), len(pick[u])))
    for comp, a, b, lab in cuts:
        al = align(comp)
        g = "%s@%s:%d-%d" % (lab, comp, a, b)
        kept = 0
        with open(os.path.join(out, "pieces", g + ".fa"), "w") as fo:
            fo.write(">%s:%d-%d|consensus\n%s\n" % (comp, a, b, cons[comp][a - 1:b]))
            for c in pick[comp]:
                if c not in al:
                    continue
                qs, qe, ss, se, m = al[c]
                if min(b, se) - max(a, ss) + 1 < 0.8 * (b - a + 1):
                    continue
                inside = [m[p] for p in range(max(a, ss), min(b, se) + 1) if m.get(p)]
                if not inside:
                    continue
                x, y = min(inside), max(inside)
                fo.write(">%s:%d-%d|%s\n%s\n" % (c, x, y, g, seqs[c][x - 1:y]))
                kept += 1
        gt.write("%s\t%s\t%s\t%d-%d\t%d\t%d\n" % (g, lab, comp, a, b, len(pick[comp]), kept))
print("done", out)
