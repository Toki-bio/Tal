#!/usr/bin/env python3
"""Build chiroptera.html, the Chiroptera SINE project page, grouped by clade.

Inputs (all in the Tal repo):
  chiroptera/genomes.tsv        code, species, family, accession, assembly, suborder, superfamily,
                                prior_evidence; rows in the order of Hao et al. 2023 Fig. 3, top to
                                bottom. code "-" = a family with no genome searched (a note, no card).
  <code>/summary.by_subfam.tsv  SINEderella step3 summary, once a Rhin-1 + VES run is published
  <code>/alignments/*.aln.fa    published alignments (SubFam input, de novo candidates)
  LEGACY below                  the four genomes searched before the one-per-family round, whose
                                published files do not follow that layout

The whole page is regenerated; run this after every publish instead of editing the HTML by hand.
"""
import csv
import glob
import html
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAGE = os.path.join(ROOT, "chiroptera.html")
SINES = ["Rhin-1", "VES"]
DENOVO_PENDING = {"hla", "tbr"}  # de novo chains launched 2026-09-26

# (copies, firm, sim median) per SINE; None = not searched. Earlier runs: Rhin-1 was the
# 177 bp sanitised consensus (IUPAC codes deleted), and rsi/rre/rda had no VES in the bank.
LEGACY = {
    "rsi": {"Rhin-1": (30079, 23545, 0.48), "VES": None, "tag": ("tag-ok", "full report"),
            "note": "Rhin-1-only search; peel bank r1&ndash;r10: 10 subfamilies &middot; 67,162 assigned copies",
            "acts": [("btn", "rsi/report.html", "<i>Rhinolophus sinicus</i> &mdash; full report"),
                     ("msa", "rsi/alignments/rsi_consensuses.aln.fa", "rsi all consensi (r1-r10 + Rhin-1)", "All Consensi in MSA"),
                     ("msa", "rhin/alignments/rsi_subfam_input_30k.aln.fa", "rsi Rhin-1 SubFam input, 600 chunk consensi + anchor", "SubFam input (600 chunks)"),
                     ("sec", "rhin.html", "Peel alignments")]},
    "rre": {"Rhin-1": (32410, 23073, 0.48), "VES": None, "tag": ("tag-partial", "pending manual review"),
            "note": "Rhin-1-only search",
            "acts": [("msa", "rhin/alignments/rre_subfam_input_30k.aln.fa", "rre Rhin-1 SubFam input, 600 chunk consensi + anchor", "SubFam input (600 chunks)")]},
    "rda": {"Rhin-1": (30169, 23850, 0.48), "VES": None, "tag": ("tag-partial", "pending manual review"),
            "note": "Rhin-1-only search",
            "acts": [("msa", "rhin/alignments/rda_subfam_input_30k.aln.fa", "rda (hap1) Rhin-1 SubFam input, 600 chunk consensi + anchor", "SubFam input (600 chunks)")]},
    "lly": {"Rhin-1": (290, 290, 0.16), "VES": (11, 11, 0.14), "tag": ("tag-partial", "Rhin-1 + VES"),
            "note": "fragmented assembly (1.9 M scaffolds); <i>Macroderma</i> (mgi) is the better Megadermatidae genome",
            "acts": [("btn", "chiroptera/lly/report.html", "<i>Lyroderma lyra</i> &mdash; report"),
                     ("msa", "chiroptera/lly/lly_subfam_input.aln.fa", "lly SubFam input, 6 chunk consensi + Rhin-1 + VES", "SubFam input (6 chunks)")]},
}


def fmt(n):
    return f"{int(n):,}"


def esc(s):
    return html.escape(s, quote=True)


def msa_btn(path, title, label, cls="secondary"):
    return (f'<a class="btn {cls}" href="javascript:void(0)" '
            f"onclick=\"openMSA('{path}','{esc(title)}')\">{label}</a>")


def load(g):
    """Return (status, {sine: (copies, firm, sim) or None}, tag, note, actions)."""
    c, sp = g["code"], g["species"]
    if c in LEGACY:
        L = LEGACY[c]
        acts = []
        for a in L["acts"]:
            if a[0] == "btn":
                acts.append(f'<a class="btn" href="{a[1]}">{a[2]}</a>')
            elif a[0] == "sec":
                acts.append(f'<a class="btn secondary" href="{a[1]}">{a[2]}</a>')
            else:
                acts.append(msa_btn(a[1], a[2], a[3]))
        return "done", {s: L[s] for s in SINES}, L["tag"], L["note"], acts
    p = os.path.join(ROOT, c, "summary.by_subfam.tsv")
    if not os.path.exists(p):
        return "running", {}, ("tag-no", "running"), "", []
    with open(p, encoding="utf8") as fh:
        rs = {r["subfam"]: r for r in csv.DictReader(fh, delimiter="\t")}
    vals = {s: ((int(rs[s]["total_assigned"]), int(rs[s]["firm_assigned"]), float(rs[s]["sim_median"]))
                if s in rs else (0, 0, None)) for s in SINES}
    acts = [f'<a class="btn" href="{c}/report.html"><i>{sp}</i> &mdash; full report</a>']
    sub = f"{c}/alignments/{c}_subfam_input_30k.aln.fa"
    if os.path.exists(os.path.join(ROOT, sub)):
        n = sum(1 for l in open(os.path.join(ROOT, sub), encoding="utf8") if l.startswith(">")) - 2
        acts.append(msa_btn(sub, f"{c} ({sp}) SubFam input, {n} chunk consensi + Rhin-1 + VES",
                            f"SubFam input ({n} chunk{'s' if n != 1 else ''})"))
    for dn in sorted(glob.glob(os.path.join(ROOT, c, "alignments", f"{c}_denovo_candidates_*chunks.aln.fa"))):
        k = re.search(r"_(\d+)chunks", dn).group(1)
        acts.append(msa_btn(f"{c}/alignments/{os.path.basename(dn)}",
                            f"{c} de novo candidates (AnnoSINE2 + fragment scan), {k} chunk consensi - for manual review",
                            f"De novo candidates ({k} chunks)", "green"))
    return "done", vals, ("tag-ok", "full report"), "", acts


def tal_style():
    """The Tal main page's own <style> block, copied verbatim so both pages share one design."""
    idx = open(os.path.join(ROOT, "index.html"), encoding="utf8").read()
    return idx[idx.index("<style>"):idx.index("</style>") + len("</style>")]


def copies_tag(v, status):
    if status == "running":
        return '<span class="tag tag-no">running</span>'
    if v is None:
        return '<span class="tag tag-no">not searched</span>'
    if v[0] == 0:
        return '<span class="tag tag-no">0</span>'
    cls = "tag-ok" if v[0] >= 1000 else "tag-partial"
    return f'<span class="tag {cls}">{fmt(v[0])} &middot; sim {v[2]:.2f}</span>'


def main():
    with open(os.path.join(ROOT, "chiroptera", "genomes.tsv"), encoding="utf8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    groups = []  # (suborder, superfamily) in tree order
    for g in rows:
        k = (g["suborder"], g["superfamily"])
        if k not in groups:
            groups.append(k)

    sections, trs = [], []
    n_done = n_gen = 0
    for so, sf in groups:
        cards, notes = [], []
        first = True
        for g in rows:
            if g["superfamily"] != sf:
                continue
            c, sp, fam = g["code"], g["species"], g["family"]
            if c == "-":
                notes.append(f"{fam} is not searched here: VES was described from it, so it was left out.")
                continue
            border = ' style="border-top:3px solid #2d2d8f;"' if first else ""
            first = False
            n_gen += 1
            status, vals, tag, note, acts = load(g)
            n_done += status == "done"
            prior = esc(g["prior_evidence"])

            meta = []
            if status == "done":
                meta.append(" &middot; ".join(f"{s} {fmt(v[0])}" for s, v in vals.items() if v is not None)
                            + " copies")
            meta.append(f'<a href="https://www.ncbi.nlm.nih.gov/datasets/genome/{g["accession"]}/" '
                        f'target="_blank">{g["accession"]}</a>')
            extra = ([f"prior evidence: {prior}"] if prior else []) + ([note] if note else [])
            card = (f'    <div class="sp-card">\n'
                    f'      <h3>{c} <span class="tag {tag[0]}">{tag[1]}</span></h3>\n'
                    f'      <div class="sci">{sp} &mdash; {fam}</div>\n'
                    f'      <div style="font-size:.82rem;color:var(--muted);margin-top:4px;">\n'
                    f'        {" &middot; ".join(meta)}\n'
                    f'      </div>\n')
            if extra:
                card += (f'      <div style="font-size:.82rem;color:var(--muted);margin-top:4px;">'
                         f'{"; ".join(extra)}</div>\n')
            if acts:
                card += '      <div class="actions">\n        ' + "\n        ".join(acts) + '\n      </div>\n'
            cards.append(card + '    </div>\n')

            href = (next((a[1] for a in LEGACY[c]["acts"] if a[0] == "btn"), None) if c in LEGACY
                    else f"{c}/report.html")
            rep = (f'<span class="tag tag-ok"><a href="{href}">HTML</a></span>' if status == "done" and href
                   else '<span class="tag tag-no">&mdash;</span>')
            if status != "done":
                aln = '<span class="tag tag-no">&mdash;</span>'
            elif c in ("rre", "rda", "lly"):
                aln = '<span class="tag tag-partial">SubFam input</span>'
            else:
                aln = '<span class="tag tag-ok">top100 · rand100 · subfam</span>'
            if glob.glob(os.path.join(ROOT, c, "alignments", f"{c}_denovo_candidates_*chunks.aln.fa")):
                dn = '<span class="tag tag-ok">ready</span>'
            elif c in DENOVO_PENDING:
                dn = '<span class="tag tag-partial">running</span>'
            else:
                dn = '<span class="tag tag-no">&mdash;</span>'
            search = ('<span class="tag tag-ok">✓</span>' if status == "done"
                      else '<span class="tag tag-no">running</span>')
            trs.append(
                f'      <tr{border}>\n'
                f'        <td><strong>{c}</strong> <em>{sp}</em> '
                f'<span class="tag tag-other-family" style="margin-left:4px;">{fam}</span></td>\n'
                f'        <td>{prior or "&mdash;"}</td>\n'
                f'        <td>{search}</td>\n'
                f'        <td>{copies_tag(vals.get("Rhin-1"), status)}</td>\n'
                f'        <td>{copies_tag(vals.get("VES"), status)}</td>\n'
                f'        <td>{rep}</td>\n'
                f'        <td>{aln}</td>\n'
                f'        <td>{dn}</td>\n'
                f'      </tr>')
        note_html = "".join(f'  <p style="color:var(--muted);font-size:.9rem;margin-top:4px;">{n}</p>\n'
                            for n in notes)
        sections.append(
            f'<section class="card">\n'
            f'  <h2>{sf} &mdash; {so}</h2>\n'
            f'{note_html}'
            f'  <div class="species-grid">\n\n' + "\n".join(cards) + '\n  </div>\n</section>\n')

    page = f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Chiroptera SINE — bat analysis</title>
{tal_style()}
</head><body>
<header>
  <h1>Chiroptera SINE — bat analysis</h1>
  <div class="sub">Rhin-1 and VES (SINEbase) searched in one genome per bat family, plus three <i>Rhinolophus</i> genomes.
  Species are grouped by clade following the family tree of Hao et al. 2023 (<i>Integrative Zoology</i> 19:989&ndash;998, Fig. 3).
  <a href="index.html" style="color:#fff;">&larr; Tal SINE main page</a></div>
</header>
<main>

<section class="card" style="padding:16px 24px;">
  <h2 style="margin-bottom:10px;">Cross-species resources</h2>
  <p style="margin:0 0 10px;font-size:.9rem;color:var(--muted);">Bank: SINEbase Rhin-1 (182 bp) and VES (220 bp), IUPAC codes resolved to a base. The earlier <i>Rhinolophus</i> runs used Rhin-1 alone (177 bp, codes deleted).</p>
  <a class="btn secondary" href="chiroptera/LOG.md">Analysis Log</a>
  <a class="btn secondary" href="rhin.html">Rhin-1 peel alignments (<i>R. sinicus</i>)</a>
</section>

{chr(10).join(sections)}
<section class="card">
  <h2>Analysis Status</h2>
  <table class="tbl">
    <thead>
      <tr><th>Species</th><th>Prior evidence</th><th>Search</th><th>Rhin-1 &middot; copies</th><th>VES &middot; copies</th><th>Report</th><th>Alignments</th><th>De novo</th></tr>
    </thead>
    <tbody>
{chr(10).join(trs)}
    </tbody>
  </table>
  <p style="color:var(--muted);font-size:.82rem;margin-top:8px;">
    {n_done} of {n_gen} genomes published. Rows are in tree order; a heavy line starts each superfamily.
    Copies = firm + soft assigned; sim = median copy bitscore / consensus self-bitscore; green = 1,000 copies or more.
    Prior evidence is the SINE label on the annotated copy of the Hao et al. figure.
    De novo = AnnoSINE2 + SINEbase-fragment scan candidates, clustered by SubFam for manual review (hla and tbr first).
  </p>
</section>

</main>
<script>
const RAW = 'https://raw.githubusercontent.com/Toki-bio/Tal/main/';
const MSA = 'https://toki-bio.github.io/MSA-viewer/';
function openMSA(relPath, title) {{
  const url = RAW + relPath;
  window.open(MSA + '?url=' + encodeURIComponent(url) + '&title=' + encodeURIComponent(title), '_blank');
}}
</script>
</body></html>
"""
    open(PAGE, "w", encoding="utf8", newline="\n").write(page)
    print(f"{n_done}/{n_gen} published")


if __name__ == "__main__":
    sys.exit(main())
