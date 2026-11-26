#!/usr/bin/env python3
"""Copy a SINEderella report and its alignments into the public tree without licensed Repbase sequence.

    python strip_repbase.py <report.html> <alignments dir> <repbase.fa> <NAME> <out dir> [OTHER_NAME ...]

NAME is the Repbase query of the run (here SINE-3_Mad); OTHER_NAME are Repbase entries whose rows were carried
over from an earlier run (here SINE-1_Mad). Removed: their consensus rows and NAME_extended in every alignment,
and, in the report, the sequence string of NAME (replaced by N of the same length, so that the drawings keep
their coordinates). Genome copies and the SubFam chunk consensuses made from them stay. The script stops if a
whole removed entry is left anywhere outside genome copies, or any 60 bp window of NAME in the report or in a
consensus row.
"""
import re
import sys
from pathlib import Path

report, aln_dir, repbase, name, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4], Path(sys.argv[5])
also = sys.argv[6:]


def fasta(p):
    d = []
    for line in p.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            d.append([line[1:], ""])
        elif d:
            d[-1][1] += line.strip()
    return d


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


ref = next(s for n, s in fasta(repbase) if n.split()[0] == name).upper()
windows = [ref[i:i + 60] for i in range(0, len(ref) - 59, 5)]
windows += [rc(w) for w in windows]


def leaks(text):
    u = text.upper()
    return sum(w in u for w in windows)


drop = {name, name + "_extended"} | set(also)
whole = {n.split()[0]: s.upper() for n, s in fasta(repbase)}
whole = {n: whole[n] for n in {name} | set(also)}
(out / "alignments").mkdir(parents=True, exist_ok=True)
removed = []
for f in sorted(aln_dir.iterdir()):
    if f.suffix != ".fa":
        if f.stat().st_size and f.suffix == ".tsv":
            (out / "alignments" / f.name).write_bytes(f.read_bytes())
        continue
    rows = fasta(f)
    keep = [r for r in rows if r[0].split()[0].split("|")[0] not in drop]
    if len(keep) != len(rows):
        removed.append((f.name, len(rows) - len(keep)))
        cols = [i for i in range(len(keep[0][1])) if any(r[1][i] not in "-." for r in keep)]
        keep = [[n, "".join(s[i] for i in cols)] for n, s in keep]
        (out / "alignments" / f.name).write_bytes("".join(f">{n}\n{s}\n" for n, s in keep).encode())
    else:
        (out / "alignments" / f.name).write_bytes(f.read_bytes())
    for n, s in keep:
        u = s.replace("-", "").replace(".", "").upper()
        if not re.match(r"(_R_)?NC_", n) and any(v in u or rc(v) in u for v in whole.values()):
            sys.exit(f"whole Repbase entry left in {f.name} row {n}")
        if not re.match(r"(_R_)?(NC_|input_)", n) and leaks(u):
            sys.exit(f"Repbase sequence left in {f.name} row {n}")

html = report.read_text(encoding="utf-8")
n_replaced = 0


def blank(m):
    global n_replaced
    if leaks(m.group(0).replace("-", "")):
        n_replaced += 1
        return "N" * len(m.group(0))
    return m.group(0)


html = re.sub(r"[ACGTNRYKMSWacgtnrykmsw\-]{60,}", blank, html)
if leaks(html) or any(v in html.upper() or rc(v) in html.upper() for v in whole.values()):
    sys.exit("Repbase sequence left in the report")
old_label = "and 9 seeds"
if old_label in html:
    html = html.replace(old_label, "and the query consensuses (Repbase rows not published)")
(out / "report.html").write_bytes(html.encode("utf-8"))
print("alignment rows removed:", removed)
print("report strings blanked:", n_replaced)
print("files:", len(list((out / "alignments").iterdir())))
