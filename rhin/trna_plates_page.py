"""trna_plates_page.py PLATE_DIR REPO_PATH TITLE_HTML - index.html for a trna_plates/ folder (SINEderella tools/trna_plates.py
output: index.tsv + *.trna.aln.fa): short description + one row per consensus with its ViewAlign link. Tied isotypes
(index.tsv tie_gene) are shown as a second gene in the same row."""
import csv, html, os, re, sys, urllib.parse as u
d, repo_path, title = sys.argv[1:4]
repo_path = repo_path.strip("/")
RAW = "https://raw.githubusercontent.com/Toki-bio/Tal/main/" + repo_path + "/"
BLOB = "https://github.com/Toki-bio/Tal/blob/main/" + os.path.dirname(repo_path) + "/README.md"
MSA = "https://toki-bio.github.io/MSA-viewer/"
style = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "chiroptera.html"), encoding="utf-8").read()
style = style[style.index("<style>"):style.index("</style>") + 8]
short = lambda g: g.replace("Homo_sapiens_", "")
rows, ties = [], []
for r in csv.DictReader(open(os.path.join(d, "index.tsv"), encoding="utf-8"), delimiter="\t"):
    tie = r.get("tie_gene", "-") not in ("-", "")
    genes = short(r["tRNA_gene"]) + (" + " + short(r["tie_gene"]) if tie else "")
    link = MSA + "?url=" + u.quote(RAW + r["plate"], safe="") + "&title=" + u.quote("%s vs %s" % (r["consensus"], genes), safe="")
    tag = "<span class='tag tag-ok'>head</span>" if r["call"] == "head" else "<span class='tag tag-no'>below shuffles</span>"
    if tie:
        tag += " <span class='tag tag-partial'>isotype tied</span>"
        ties.append("%s (%s %s bits, %s %s)" % (r["consensus"], short(r["tRNA_gene"]).split("-")[1], r["bits"],
                                                short(r["tie_gene"]).split("-")[1], r["tie_bits"]))
    two = lambda a, b: a + ("<br>" + b if tie else "")
    rows.append("<tr><td><a href='%s' target='_blank'>%s</a></td><td style='white-space:nowrap'>%s</td><td class='num'>%s</td><td class='num'>%s</td>"
                "<td class='num'>%s</td><td class='num'>%s</td><td class='num'>%s</td><td>%s</td>"
                "<td><a class='btn' href='%s' target='_blank'>ViewAlign</a> <a href='%s'>fasta</a></td></tr>" % (
                    html.escape(link), html.escape(r["consensus"]),
                    two(html.escape(short(r["tRNA_gene"])), html.escape(short(r.get("tie_gene", "")))),
                    two(r["ident%"] + " %", r.get("tie_ident%", "") + " %"), two(r["aln_bp"], r.get("tie_aln_bp", "")),
                    r["q_span"], two(r["bits"], r.get("tie_bits", "")), r["null_max_bits"], tag, html.escape(link),
                    html.escape(r["plate"])))
tie_txt = ("Tied here: " + "; ".join(ties) + ".") if ties else "No tied isotypes in this run."
page = """<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>%s</title>
%s
</head><body>
<header><h1>%s</h1><div class="sub">Each subfamily consensus against its closest tRNA gene &middot; %s</div></header>
<main>
<section class="card">
<h2>What these plates are</h2>
<p>Each plate shows one subfamily consensus of the run and the human tRNA gene (GtRNAdb hg19, 419 genes) it matches best. They are aligned by the same ssearch36 local alignment that scored the head in SINEderella's tRNA step (<code>tools/trna_heads.sh</code>, plates by <code>tools/trna_plates.py</code>). The aligned stretch is in upper case and the rest in lower case.</p>
<p><b>Isotype tied</b>: when the best gene of a different amino acid scores within 2 bits of the best, that gene is a third row of the plate and a second line in the table. %s</p>
<p><b>Shuffles</b> is the best score of 20 dinucleotide shuffles of the same consensus against the same genes; a head counts when its score is above that. MEG-RL and MEG-RS are 5S rRNA-derived (Rfam 5S_rRNA, positions 1-119), so their best tRNA hit falls below the shuffles. A composite's plate shows its best-scoring head only.</p>
</section>
<section class="card">
<h2>Plates</h2>
<table class="tbl"><thead><tr><th>subfamily</th><th>best tRNA gene</th><th>identity</th><th>aligned nt</th><th>consensus span</th><th>bits</th><th>shuffles max</th><th>call</th><th>open</th></tr></thead>
<tbody>
%s
</tbody></table>
<p><a href="%s">run README</a> &middot; <a href="index.tsv">index.tsv</a> &middot; <a href="../../chiroptera.html">Chiroptera page</a></p>
</section>
</main></body></html>
""" % (html.escape(html.unescape(re.sub("<[^>]+>", "", title))), style, title, html.escape(os.path.dirname(repo_path)),
       html.escape(tie_txt), "\n".join(rows), BLOB)
open(os.path.join(d, "index.html"), "w", encoding="utf-8", newline="\n").write(page)
print("wrote", os.path.join(d, "index.html"), len(rows), "rows,", len(ties), "tied")
