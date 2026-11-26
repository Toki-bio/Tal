#!/usr/bin/env python3
"""Build ../../plants.html from the files in this folder (run from anywhere).

Inputs: alignments/{mdo_consensuses.aln.fa, mdo_seeds.names.tsv, counts.txt}, genome.fai, ncbi_assembly.json.
The first three come from Biotech2024:/staging/tmp/apple_sines/publish/, genome.fai from
/staging/tmp/apple_sines/genome/apple.fna.fai, ncbi_assembly.json from the NCBI Datasets API.
Also normalises the plate once: the aligner's _R_ prefix is moved from sequence names into the table.
"""
import csv
import json
import urllib.parse
from pathlib import Path

H = Path(__file__).resolve().parent
A = H / "alignments"
REPO = H.parent.parent

plate = (A / "mdo_consensuses.aln.fa").read_text(encoding="utf-8").splitlines()
rev = {l[4:].split()[0] for l in plate if l.startswith(">_R_")}
if rev:
    plate = [">" + l[4:] if l.startswith(">_R_") else l for l in plate]
    (A / "mdo_consensuses.aln.fa").write_bytes(("\n".join(plate) + "\n").encode())
rows = list(csv.DictReader(open(A / "mdo_seeds.names.tsv", encoding="utf-8"), delimiter="\t"))
if "reverse_complemented_in_plate" not in rows[0]:
    for r in rows:
        r["reverse_complemented_in_plate"] = "yes" if r["name"] in rev else "no"
    with open(A / "mdo_seeds.names.tsv", "w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

c = dict(l.split() for l in (A / "counts.txt").read_text(encoding="utf-8").splitlines() if l.strip())
names = [l[1:] for l in plate if l.startswith(">")]
assert len(names) == len(set(names)) == len(rows) == int(c["seeds"])
assert set(names) == {r["name"] for r in rows}
assert c["mismatch"] == "0"
assert sum(r["reverse_complemented_in_plate"] == "yes" for r in rows) == int(c["reverse_complemented"])
fai = [l.split("\t") for l in (H / "genome.fai").read_text(encoding="utf-8").splitlines()]
total = sum(int(x[1]) for x in fai)
nc = json.load(open(H / "ncbi_assembly.json", encoding="utf-8"))
organelle = total - nc["total_sequence_length"]


def n(v):
    return f"{int(v):,}"


raw = "https://raw.githubusercontent.com/Toki-bio/Tal/main/plants/mdo/alignments/mdo_consensuses.aln.fa"
title = f"Apple (mdo): {c['seeds']} AnnoSINE_v2 seed consensuses"
msa = ("https://toki-bio.github.io/MSA-viewer/?url=" + urllib.parse.quote(raw, safe="")
       + "&title=" + urllib.parse.quote(title, safe=""))
tags = [("tRNA", "transfer RNA (tRNA)"), ("5S_rRNA", "5S ribosomal RNA (rRNA)"),
        ("7SL_RNA", "7SL RNA"), ("Unknown", "unknown")]
assert sum(int(c.get("tag_" + k, 0)) for k, _ in tags) == int(c["seeds"])
tagrows = "".join(f"<tr><td>{lab}</td><td class='num'>{n(c['tag_' + k])}</td></tr>" for k, lab in tags)
mm = sorted(((k[6:], int(v)) for k, v in c.items() if k.startswith("match_")), key=lambda x: -x[1])
matchrows = "".join(f"<tr><td>{k}</td><td class='num'>{n(v)}</td></tr>" for k, v in mm)
assert sum(v for _, v in mm) == int(c["matched"])
unmatched = int(c["seeds"]) - int(c["matched"])

peeled = []
for line in (H / "peel" / "peel_round1_owner_consensuses.fa").read_text(encoding="utf-8").splitlines():
    if line.startswith(">"):
        peeled.append([line[1:].strip(), 0])
    elif line.strip():
        peeled[-1][1] += len(line.strip())
peeled_sum = sum(int(k.split("_")[1].replace("seqs", "")) for k, _ in peeled)
peelrows = "".join(f"<tr><td>{k}</td><td class='num'>{v}</td></tr>" for k, v in peeled)

SD = H / "sd"
summ = list(csv.DictReader(open(SD / "summary.by_subfam.tsv", encoding="utf-8"), delimiter="\t"))
verd = {}
for line in (SD / "verdict_counts.tsv").read_text(encoding="utf-8").splitlines():
    fam, v, k = line.split("\t")
    verd.setdefault(fam, {})[v] = int(k)
tot = {v: sum(d.get(v, 0) for d in verd.values()) for v in ("SINE", "grey", "not_SINE")}
copies = sum(tot.values())
base = "https://raw.githubusercontent.com/Toki-bio/Tal/main/plants/mdo/sd/alignments/"


def view(fname, label, cls=""):
    if not (SD / "alignments" / fname).is_file():
        return ""
    u = ("https://toki-bio.github.io/MSA-viewer/?url=" + urllib.parse.quote(base + fname, safe="")
         + "&title=" + urllib.parse.quote("mdo " + fname.replace(".aln.fa", ""), safe=""))
    return f"<a class='btn{cls}' href='{u.replace('&', '&amp;')}' target='_blank'>{label}</a>"


resrows = ""
for r in sorted(summ, key=lambda x: -verd[x["subfam"]].get("SINE", 0)):
    s, d = r["subfam"], verd[r["subfam"]]
    assert sum(d.values()) >= int(r["firm_assigned"])
    plates = (view(f"mdo_{s}_top100.aln.fa", "top copies") + view(f"mdo_{s}_rand100.aln.fa", "random copies", " secondary")
              + view(f"mdo_{s}_subfam.aln.fa", "chunk consensuses", " green"))
    resrows += (f"<tr><td>{s}</td><td class='num'>{n(r['firm_assigned'])}</td><td class='num'>{n(d.get('SINE', 0))}</td>"
                f"<td class='num'>{n(d.get('grey', 0))}</td><td class='num'>{n(d.get('not_SINE', 0))}</td>"
                f"<td class='num'>{r['leak_pct']}</td><td>{plates}</td></tr>")
cont = [x for x in csv.DictReader(open(SD / "continuation_check.tsv", encoding="utf-8"), delimiter="\t")
        if x["consensus"] == "consensus_3seqs" and x["flag"] == "yes"]
selfnote = ""
if cont:
    x = cont[0]
    selfnote = (f"<p><b>consensus_3seqs needs caution.</b> In {n(x['copies'])} of {n(x['tested'])} tested assigned copies the "
                f"sequence continues at the {x['side'].replace(chr(39), '&prime;')} end into the same consensus (dimer or "
                f"tandem copies), and {n(verd['consensus_3seqs'].get('grey', 0))} of its copies are in the grey zone.</p>")

css = (REPO / "echinoderms.html").read_text(encoding="utf-8")
css = css[css.index("<style>"):css.index("</style>")]
old = "main { max-width:1100px; margin:0 auto; padding:28px 24px; }"
assert old in css
css = css.replace(old, "main { margin:0; padding:28px 24px; }\n"
                       "main p, main li { max-width:120ch; }\ncode { overflow-wrap:anywhere; }") + "</style>"

html = f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Plant SINEs — apple</title>
{css}
</head><body>
<header>
  <h1>Plant SINEs &mdash; apple (<i>Malus domestica</i>)</h1>
  <div class="sub">Search for short interspersed nuclear elements (SINEs) with AnnoSINE_v2 and SINEderella. Work in progress, state of 2026-11-26.
  <a href="index.html" style="color:#fff;">&larr; Tal SINE main page</a></div>
</header>
<main>

<section class="card">
  <h2>SINEderella result, whole genome</h2>
  <p>{n(copies)} extracted copies. Verdicts: <b>{n(tot['SINE'])} SINE</b> (assigned to one query, flanks unique), {n(tot['grey'])} grey zone (flanks could not be tested), {n(tot['not_SINE'])} not SINE (not assigned, or flanks shared with other copies).</p>
  <a class="btn" href="plants/mdo/sd/report.html">full SINEderella report</a>
  {view('mdo_consensuses.aln.fa', 'consensus alignment', ' secondary')} {view('mdo_allseeds_chunks.aln.fa', 'all chunk consensuses with the queries', ' secondary')}
  <div style="overflow-x:auto;"><table class="tbl"><thead><tr><th>Query</th><th>Assigned</th><th>SINE</th><th>Grey zone</th><th>Not SINE</th><th>Leak, % of assigned</th><th>Alignments of copies (ViewAlign)</th></tr></thead><tbody>{resrows}</tbody></table></div>
  <p>Assigned: copies that the same query won in all ten voting rounds, above the score threshold. Leak: share of assigned copies whose second-best query scores at least 90&nbsp;% of the best. Grey-zone and not-SINE counts include copies that were not firmly assigned. Repbase consensus rows are not shown in the published alignments.</p>
  {selfnote}
  <p>Sources: <code>Biotech2024:/staging/tmp/apple_sines/sd_peel9/run_excl_20261126_020256/results/summary.by_subfam.tsv</code> and <code>.../step2/step2_output/loci_verdict.tsv</code>.</p>
</section>

<section class="card">
  <h2>Genome</h2>
  <div style="overflow-x:auto;"><table class="tbl">
  <thead><tr><th>Code</th><th>Species</th><th>Accession</th><th>Contig N50</th><th>Sequences</th><th>Size searched, base pairs (bp)</th></tr></thead>
  <tbody><tr><td>mdo</td><td><i>Malus domestica</i> cv. Golden Delicious</td>
  <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/{nc['accession']}/">{nc['accession']}</a> ({nc['assembly_name']}, one haplotype)</td>
  <td class="num">{nc['contig_n50'] / 1e6:.2f} Mb</td><td class="num">{len(fai)}</td><td class="num">{n(total)}</td></tr></tbody></table></div>
  <p>Accession and contig N50: National Center for Biotechnology Information (NCBI). The size searched includes the two organelle genomes ({n(organelle)} bp) on top of the NCBI genome total. Mb: million base pairs. Server copy: <code>Biotech2024:/staging/tmp/apple_sines/genome/apple.fna</code>.</p>
</section>

<section class="card">
  <h2>Consensus plate</h2>
  <p>All {n(c['seeds'])} AnnoSINE_v2 seed consensuses in one alignment ({n(c['columns'])} columns). The seeds are redundant; the consensuses peeled from this plate are listed below. {n(c['reverse_complemented'])} sequences are shown reverse-complemented so that they line up with the rest; the seed table records which. Similar sequences are placed next to each other.</p>
  <a class="btn" href="{msa.replace('&', '&amp;')}" target="_blank">Consensus plate in ViewAlign</a>
  <a class="btn secondary" href="plants/mdo/alignments/mdo_consensuses.aln.fa">alignment (FASTA)</a>
  <a class="btn secondary" href="plants/mdo/alignments/mdo_seeds.names.tsv">seed table: location, length, head type, match</a>
</section>

<section class="card">
  <h2>AnnoSINE_v2 result</h2>
  <p>{n(c['seeds'])} seed consensuses, {c['len_min']}&ndash;{c['len_max']} bp, median {c['len_median']} bp.</p>
  <div class="species-grid">
  <div><table class="tbl"><thead><tr><th>Predicted origin of the head</th><th>Seeds</th></tr></thead><tbody>{tagrows}</tbody></table></div>
  <div><table class="tbl"><thead><tr><th>Apple Repbase entry matched</th><th>Seeds</th></tr></thead><tbody>{matchrows}<tr><td>no query entry matched at this cut</td><td class="num">{n(unmatched)}</td></tr></tbody></table></div>
  </div>
  <p>Match: ssearch36, either strand, at least 80&nbsp;% identity over at least 80&nbsp;% of the seed; the best-scoring such entry is counted. Sources: <code>Biotech2024:/staging/tmp/apple_sines/anno_full/out/Seed_SINE.fa</code> and <code>Biotech2024:/staging/tmp/apple_sines/publish/counts.txt</code>.</p>
</section>

<section class="card">
  <h2>Peeled consensuses</h2>
  <p>{len(peeled)} consensuses taken from the plate by hand; each name carries the number of seeds in its group ({n(peeled_sum)} of {n(c['seeds'])} in total).</p>
  <div style="overflow-x:auto;"><table class="tbl"><thead><tr><th>Consensus</th><th>Length, bp</th></tr></thead><tbody>{peelrows}</tbody></table></div>
  <p></p>
  <a class="btn secondary" href="plants/mdo/peel/peel_round1_owner_consensuses.fa">consensuses (FASTA)</a>
  <a class="btn secondary" href="https://github.com/Toki-bio/Tal/blob/main/plants/mdo/peel/README.md" target="_blank">peel log</a>
</section>

<section class="card">
  <h2>Workflow</h2>
  <ol>
    <li><b>AnnoSINE_v2</b> 2.0.9, mode 3 (structure and homology search), plant profiles, whole genome, RepeatMasker step off.</li>
    <li><b>Comparison set from Repbase</b> (the repeat library of the Genetic Information Research Institute): 84 land-plant SINE entries, 83 of 95&ndash;369 bp plus apple SINE-3_Mad (797 bp), used for the match table above. Repbase sequences are licensed and are not published here.</li>
    <li><b>Peeling</b> of the seed plate into groups, one consensus per group.</li>
    <li><b>SINEderella</b> (commit b5c8a4a): search, copy assignment, flank test, alignments and report, with {len(peeled) + 1} queries: the {len(peeled)} peeled consensuses and Repbase SINE-3_Mad. Repbase SINE-1_Mad was removed in favour of the nearly identical consensus_354seqs, and the copies were assigned again.</li>
  </ol>
</section>

<section class="card">
  <h2>Analysis Status</h2>
  <table class="tbl">
  <tr><td>Genome and tools</td><td>done</td></tr>
  <tr><td>Test of both tools on the first 10 million bp of one chromosome</td><td>done; apple Repbase SINEs found</td></tr>
  <tr><td>AnnoSINE_v2, whole genome</td><td>done, {n(c['seeds'])} seeds</td></tr>
  <tr><td>Peeling of the seed plate</td><td>round 1 done, {len(peeled)} consensuses</td></tr>
  <tr><td>SINEderella, whole genome</td><td>done: {n(tot['SINE'])} SINE loci; report and alignments published</td></tr>
  </table>
</section>

</main>
</body></html>
"""
(REPO / "plants.html").write_bytes(html.encode("utf-8"))
print("plants.html written;", len(names), "sequences")
print(msa)
