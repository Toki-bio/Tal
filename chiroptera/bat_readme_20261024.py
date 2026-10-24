"""README.md for each <code>_main_20261024 bat run: species line, code version, verdict counts per subfamily and family,
compared with the copy numbers of the published 2026-09-27 run (<code>/summary.by_subfam.tsv)."""
import csv, gzip, os, re, sys, collections
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = {r["code"]: r for r in csv.DictReader(open(os.path.join(ROOT, "chiroptera", "genomes.tsv"), encoding="utf-8"), delimiter="\t")}
FAM = {"Rhin-1": "Rhin", "VES": "VES", "MEG-RL": "MEG", "MEG-RS": "MEG", "MEG-T2": "MEG", "MEG-TR": "MEG"}
ORDER = ["Rhin-1", "VES", "MEG-RL", "MEG-RS", "MEG-T2", "MEG-TR"]
def old_counts(code):
    p = os.path.join(ROOT, code, "summary.by_subfam.tsv")
    if not os.path.exists(p): return {}
    rows = list(csv.reader(open(p, encoding="utf-8"), delimiter="\t"))
    h = rows[0]; i = h.index("subfamily") if "subfamily" in h else 0
    j = next((k for k, x in enumerate(h) if x.lower() in ("n", "count", "copies", "n_copies")), 1)
    return {r[i]: r[j] for r in rows[1:] if len(r) > j}
for code in sys.argv[1:]:
    d = os.path.join(ROOT, code + "_main_20261024"); g = G[code]
    c = collections.defaultdict(collections.Counter); tot = collections.Counter()
    for r in csv.DictReader(gzip.open(os.path.join(d, "loci_verdict.tsv.gz"), "rt", encoding="utf-8"), delimiter="\t"):
        tot["list"] += 1; tot[r["verdict"]] += 1
        if r["assignment_status"] == "assigned":
            tot["assigned"] += 1; c[r["subfamily"]][r["verdict"]] += 1; c[r["subfamily"]]["n"] += 1
    man = open(os.path.join(d, "manifest.txt"), encoding="utf-8").read()
    stages = open(os.path.join(d, "pipeline_stages.txt"), encoding="utf-8").read().splitlines()
    code_ver = re.search(r"code (\w+)", stages[0]).group(1) if stages else "?"
    t0 = stages[0].split("start ")[1][:19] if stages else "?"; t1 = stages[-1].split("done ")[1][:19] if stages and "done" in stages[-1] else "?"
    old = old_counts(code)
    L = [f"# {code}_main_20261024", "",
         f"*{g['species']}* ({g['family']}; {g['superfamily']}, {g['suborder']}), genome `{g['accession']}`. Complete SINEderella rerun on",
         f"main {code_ver}: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,",
         "MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks",
         "(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.", "",
         f"- Server: therioserver `/home/toki/chiro/{code}_main_20261024/` (run {t0} to {t1}; scripts `/home/toki/chiro/rerun_20261024/`).",
         "- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).", "",
         "| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |", "|---|---|---|---|---|---|---|"]
    fam = collections.defaultdict(collections.Counter)
    for sf in ORDER + sorted(set(c) - set(ORDER)):
        if sf not in c: continue
        k = c[sf]; F = FAM.get(sf, "?")
        for x in ("n", "SINE", "grey", "not_SINE"): fam[F][x] += k[x]
        L.append(f"| {sf} | {F} | {k['n']:,} | {k['SINE']:,} | {k['grey']:,} | {k['not_SINE']:,} | {old.get(sf, '-')} |")
    L += ["", "| family | assigned | SINE | grey | not SINE |", "|---|---|---|---|---|"]
    for F in ("Rhin", "VES", "MEG"):
        if F in fam: L.append(f"| {F} | {fam[F]['n']:,} | {fam[F]['SINE']:,} | {fam[F]['grey']:,} | {fam[F]['not_SINE']:,} |")
    L += ["", f"Copies in the list {tot['list']:,}; assigned {tot['assigned']:,}; not SINE in total (assigned + not assigned) {tot['not_SINE']:,}.",
          "Rules: SINEderella `docs/VERDICTS.md`.", "",
          "Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,",
          "`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,",
          "`pipeline_stages.txt` (stage times)."]
    open(os.path.join(d, "README.md"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
    print(code, tot["assigned"], dict((k, c[k]["n"]) for k in c))
