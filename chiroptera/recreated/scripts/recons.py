#!/usr/bin/env python3
"""Collect the consensus REBUILT in each bat species (row 1 of every top100 plate), for one alignment.

Source: ~/tmp/bw/new/<code>/<code>_<family>_top100.aln.fa - the 312-plate corpus run through the current
chain (2026-09-28). Row 1 = <family>_extended, the consensus rebuilt from the plate's copies; the part
taken is fix_alignments.judged_span (no trim proposals; additions only where the copies carry them),
degapped, case kept (lowercase = a proposed addition past the original).
Skipped: plates whose verdict has NO_ELEMENT (row 1 then is built from flank) or fewer than 5 copies.
References: the queries as searched (chiro_bank.fa Rhin-1/VES, MEG.resolved.fa, rsi r1-r10 row 2).
Writes OUT/recreated.fa, OUT/skipped.tsv, OUT/members.tsv.
"""
import glob, os, sys
sys.path.insert(0, os.path.expanduser("~/tmp/bw/code"))
from fix_alignments import read_fa, judged_span, is_seed, GAPS
from verdict import verdict

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
FAM = {"rle": "Pteropodidae", "tpe": "Rhinonycteridae", "mgi": "Megadermatidae", "cth": "Craseonycteridae",
       "rmi": "Rhinopomatidae", "mau": "Myzopodidae", "nth": "Nycteridae", "rna": "Emballonuridae",
       "tni": "Phyllostomidae", "mme": "Mormoopidae", "nle": "Noctilionidae", "fho": "Furipteridae",
       "hla": "Hipposideridae", "tbr": "Molossidae", "lly": "Megadermatidae", "mev": "Vespertilionidae",
       "vmu": "Vespertilionidae", "rsi": "Rhinolophidae", "rda": "Rhinolophidae", "rre": "Rhinolophidae"}
for l in open(os.path.expanduser("~/chiro/chiro_genomes.tsv")):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 3:
        FAM.setdefault(f[0], f[2])

recs, refs, skipped, members = [], [], [], []
for p in sorted(glob.glob(os.path.expanduser("~/tmp/bw/new/*/*_top100.aln.fa"))):
    base = os.path.basename(p)[:-len("_top100.aln.fa")]
    code, fam = base.split("_", 1)
    names, seqs = read_fa(p)
    v = verdict(p)
    flags = [x["code"] for x in v.get("flags", [])]
    n = v.get("n") or 0
    if "NO_ELEMENT" in flags or n < 5:
        skipped.append((code, fam, n, v.get("score"), ",".join(flags) or "-"))
        continue
    lo, hi = judged_span(names, seqs, 0, p)
    seq = "".join(c for c in seqs[0][lo:hi + 1] if c not in GAPS)
    o = [i for i in range(1, len(names)) if is_seed(names[i], names[0])]
    orig_len = sum(1 for c in seqs[o[0]] if c not in GAPS) if o else len(seq)
    core = ""
    # A capped plate (tandem array / small core, score < 50) whose rebuilt row is > 2.5x the original is
    # mostly the repeated unit's shared flank (cse MEG-T2 970 bp, vmu MEG-TR 1044 bp): keep its uppercase
    # core, the part over the original's span, so it does not stretch the whole alignment.
    if (v.get("score") or 0) < 50 and len(seq) > 2.5 * orig_len:
        seq = "".join(c for c in seqs[0] if c not in GAPS and c.isupper())
        core = "_core"
    name = "%s__%s_%s_n%d_s%d%s" % (fam, code, FAM.get(code, "?"), n, round(v.get("score") or 0), core)
    recs.append((name, seq))
    members.append((name, code, FAM.get(code, "?"), fam, n, v.get("score"), ",".join(flags) or "-", len(seq)))
    if code == "rsi" and fam.startswith("r"):                  # rsi peel-bank originals (row 2)
        o = [i for i in range(1, len(names)) if is_seed(names[i], names[0])]
        if o:
            refs.append(("REF_" + fam, "".join(c for c in seqs[o[0]] if c not in GAPS).upper()))

for f in ("~/chiro/chiro_bank.fa", "~/chiro/MEG.resolved.fa"):
    n, s = read_fa(os.path.expanduser(f))
    refs += [("REF_" + a.split()[0], b.upper()) for a, b in zip(n, s)]
# references FIRST: mafft --adjustdirection orients every sequence against the first one
refs.sort(key=lambda r: (r[0] != "REF_Rhin-1", r[0]))
recs = refs + recs

with open(os.path.join(OUT, "recreated.fa"), "w") as fh:
    for a, b in recs:
        fh.write(">%s\n%s\n" % (a, b))
with open(os.path.join(OUT, "skipped.tsv"), "w") as fh:
    fh.write("code\tfamily\tn\tscore\tflags\n")
    fh.writelines("\t".join(map(str, r)) + "\n" for r in skipped)
with open(os.path.join(OUT, "members.tsv"), "w") as fh:
    fh.write("name\tcode\tbat_family\tsine_family\tn\tscore\tflags\tlen\n")
    fh.writelines("\t".join(map(str, r)) + "\n" for r in members)
print("rebuilt %d, references %d, skipped %d" % (len(members), sum(a.startswith("REF_") for a, _ in recs), len(skipped)))
