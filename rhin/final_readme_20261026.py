"""README for the 2026-10-26 rhinolophid-bank runs (rsi_main_20261026, rfe_main_20261026, hla_rsibank_20261026).
Usage (from the Tal repo root): python rhin/final_readme_20261026.py <dir> "<species line>" "<extra markdown file or ->"
Counts come from the published files in <dir>: loci_verdict.tsv.gz, families.tsv, family_view/, assignment_stats.tsv,
pipeline_stages.txt, manifest.txt."""
import collections
import gzip
import os
import sys

d, species, extra = sys.argv[1], sys.argv[2], sys.argv[3]
fam = {}
for l in open(os.path.join(d, "families.tsv"), encoding="utf-8"):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 2 and f[0] != "subfamily":
        fam[f[0]] = f[1]
rows = [l.rstrip("\n").split("\t") for l in gzip.open(os.path.join(d, "loci_verdict.tsv.gz"), "rt", encoding="utf-8")]
h = rows[0]; rows = rows[1:]
ci, si, ai, vi = h.index("copy"), h.index("subfamily"), h.index("assignment_status"), h.index("verdict")
tot = collections.Counter(r[vi] for r in rows)
assigned = [r for r in rows if r[ai] == "assigned"]
per = collections.defaultdict(collections.Counter)
nsub = collections.defaultdict(set)
for r in assigned:
    F = fam.get(r[si], r[si]); per[F][r[vi]] += 1; nsub[F].add(r[si])
fv = os.path.join(d, "family_view")
pc_path = os.path.join(fv, "partial_copies.tsv")
pc_open = (lambda: gzip.open(pc_path + ".gz", "rt", encoding="utf-8")) if not os.path.exists(pc_path) else (lambda: open(pc_path, encoding="utf-8"))
pc = [l.rstrip("\n").split("\t") for l in pc_open()][1:]
st = [l.rstrip("\n").split("\t") for l in open(os.path.join(fv, "partial_stretches.tsv"), encoding="utf-8")][1:]
fvs = [l.rstrip("\n").split("\t") for l in open(os.path.join(fv, "family_view.summary.tsv"), encoding="utf-8")]
code = ""
for l in open(os.path.join(d, "manifest.txt"), encoding="utf-8"):
    if l.startswith("SINEDERELLA_COMMIT") or l.startswith("CODE"):
        code = l.split("\t", 1)[1].strip(); break
stages = open(os.path.join(d, "pipeline_stages.txt"), encoding="utf-8").read().strip().splitlines()
out = []
w = out.append
w("# %s\n" % os.path.basename(d))
w("%s Bank: the *Rhinolophus sinicus* bank without r1_9seqs (21 consensuses: r2-r10, r10_r8_P18del, seven composites, four "
  "MEG SINEs; therioserver `/home/toki/rhin/final_bank_20261026/bank21_nor1.fa`), families Rhin and MEG (`families.tsv`). "
  "SINEderella main a9d2a51 (2026-10-26): composite stretches of >= 50 bp searched as queries of their own in step 1; "
  "partial copies labelled and grouped into stretches in the family roll-up (stretch tables regrouped with main 9d484fa); "
  "flank scan stages 1-7 and the report rebuild as for rsi_main_20261022.\n" % species)
w("- Server: therioserver `/home/toki/rhin/%s/` (`../final_bank_20261026/launch_final.sh`, `pack_final.sh`)." % os.path.basename(d))
if stages:
    w("- Run: %s; %s." % (stages[0].lstrip("= ").split(" code")[0], stages[-1].lstrip("= ")))
w("")
w("## Verdicts (every copy; `loci_verdict.tsv.gz`, rules SINEderella `docs/VERDICTS.md`)\n")
w("| | copies |\n|---|---|")
w("| copies in the list | %s |" % format(len(rows), ","))
w("| assigned to a subfamily | %s |" % format(len(assigned), ","))
for k in ("SINE", "grey", "not_SINE"):
    w("| %s (all copies) | %s |" % ({"SINE": "SINE loci", "grey": "grey zone", "not_SINE": "not SINE"}[k], format(tot[k], ",")))
w("\nCopies without a subfamily are listed as not SINE in this table, partial copies included (below).\n")
w("Per family (copies assigned to a subfamily):\n")
w("| family | subfamilies with copies | assigned | SINE | grey | not SINE |\n|---|---|---|---|---|---|")
for F in sorted(per):
    c = per[F]
    w("| %s | %d | %s | %s | %s | %s |" % (F, len(nsub[F]), format(sum(c.values()), ","), format(c["SINE"], ","),
                                          format(c["grey"], ","), format(c["not_SINE"], ",")))
w("\n## Family roll-up and partial copies (`family_view/`)\n")
w("| family | copies | subfamily-firm | near tie | family only | not firm |\n|---|---|---|---|---|---|")
for r in fvs[1:]:
    w("| %s | %s | %s | %s | %s | %s |" % (r[0], format(int(r[2]), ","), format(int(r[3]), ","), format(int(r[4]), ","),
                                          format(int(r[5]), ","), format(int(r[7]), ",")))
npart = sum(1 for r in pc if r[7] == "partial")
w("\nPartial copies (family-firm, covering < 60 %% of every competing consensus): %s, plus %s below their stretch bar; "
  "%d stretches (`partial_stretches.tsv`, `partial_copies.tsv`). Largest stretches:\n" % (format(npart, ","), format(len(pc) - npart, ","), len(st)))
w("| stretch | copies | median length | spans |\n|---|---|---|---|")
for r in st[:8]:
    w("| %s | %s | %s | %s |" % (r[0], format(int(r[2]), ","), r[3], r[4].replace(",", ", ")))
if extra != "-" and os.path.exists(extra):
    w("\n" + open(extra, encoding="utf-8").read().strip())
w("\n## Contents\n| Path | What |\n|---|---|")
w("| `report.html` | full report |\n| `loci_verdict.tsv.gz`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |")
w("| `family_view/` | family roll-up per copy, partial copies, stretches, the composite stretches searched in step 1 |")
w("| `alignments/` | plates (SINE loci, topped up with grey copies marked [grey]) |")
w("| `continuation/`, `flankscan/` | continuation check; flank scan stages 1-7 |")
w("| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers |")
w("| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | flank twins, tandem-array share, satellite screen |")
w("| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its copies |")
w("| `composite_units.tsv`, `consensus_blocks.tsv` | composites read as units; the bank as blocks |")
w("| `consensuses.clean.fa`, `manifest.txt`, `pipeline_stages.txt` | the bank; tool versions and command; stage times |")
open(os.path.join(d, "README.md"), "w", encoding="utf-8", newline="\n").write("\n".join(out) + "\n")
print("README", d, len(rows), "copies")
