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


def sine_cells(vals, status):
    if status == "running":
        return '<td class="num muted" colspan="4">running</td>'
    out = []
    for s in SINES:
        v = vals.get(s)
        if v is None:
            out.append('<td class="num muted" colspan="2">not searched</td>')
        elif v[0] == 0:
            out.append('<td class="num">0</td><td class="num">&ndash;</td>')
        else:
            out.append(f'<td class="num">{fmt(v[0])}<div class="sub">{fmt(v[1])} firm</div></td>'
                       f'<td class="num">{v[2]:.2f}</td>')
    return "".join(out)


def main():
    with open(os.path.join(ROOT, "chiroptera", "genomes.tsv"), encoding="utf8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    table, sections = [], []
    n_done = n_gen = 0
    suborders = []
    for g in rows:
        if g["suborder"] not in suborders:
            suborders.append(g["suborder"])

    for so in suborders:
        table.append(f'<tr class="so"><td colspan="9">{so}</td></tr>')
        sfs = []
        for g in rows:
            if g["suborder"] == so and g["superfamily"] not in sfs:
                sfs.append(g["superfamily"])
        sf_blocks = []
        for sf in sfs:
            table.append(f'<tr class="sf"><td colspan="9">{sf}</td></tr>')
            cards, notes = [], []
            for g in rows:
                if g["superfamily"] != sf:
                    continue
                c, sp, fam = g["code"], g["species"], g["family"]
                prior = esc(g["prior_evidence"]) or "&ndash;"
                if c == "-":
                    table.append(f'<tr><td>{fam}</td><td class="muted" colspan="2">not searched here</td>'
                                 f'<td>{prior}</td><td class="muted" colspan="4">&ndash;</td><td></td></tr>')
                    notes.append(f"<i>{fam}</i> is not searched here: VES was described from it, so it was left out.")
                    continue
                n_gen += 1
                status, vals, tag, note, acts = load(g)
                n_done += status == "done"
                href = (next((a[1] for a in LEGACY[c]["acts"] if a[0] == "btn"), None) if c in LEGACY
                        else f"{c}/report.html")
                link = f'<a href="{href}">{c}</a>' if status == "done" and href else c
                dn = ""
                if glob.glob(os.path.join(ROOT, c, "alignments", f"{c}_denovo_candidates_*chunks.aln.fa")):
                    dn = '<span class="tag tag-ok">ready</span>'
                elif c in DENOVO_PENDING:
                    dn = '<span class="tag tag-no">running</span>'
                table.append(f'<tr><td>{fam}</td><td>{link}</td><td><i>{sp}</i></td><td>{prior}</td>'
                             f'{sine_cells(vals, status)}<td>{dn}</td></tr>')
                counts = ""
                if status == "done":
                    parts = [f"{s} {fmt(v[0])}" + (f" ({fmt(v[1])} firm)" if v[0] else "")
                             for s, v in vals.items() if v is not None]
                    counts = f'<div class="meta">{" &middot; ".join(parts)}</div>'
                cards.append(
                    f'    <div class="sp-card">\n'
                    f'      <h3>{c} <span class="tag {tag[0]}">{tag[1]}</span></h3>\n'
                    f'      <div class="sci">{sp} &mdash; {fam}</div>\n'
                    f'      <div class="meta"><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/{g["accession"]}/" '
                    f'target="_blank">{g["accession"]}</a> {g["assembly"]}'
                    f'{" &middot; prior evidence: " + prior if g["prior_evidence"] else ""}</div>\n'
                    f'      {counts}\n'
                    f'      {f"<div class=meta>{note}</div>" if note else ""}\n'
                    f'      <div class="actions">{" ".join(acts)}</div>\n'
                    f'    </div>')
            note_html = "".join(f'<p class="meta">{n}</p>' for n in notes)
            sf_blocks.append(f'  <h3 class="sf">{sf}</h3>\n{note_html}\n  <div class="species-grid">\n'
                             + "\n".join(cards) + "\n  </div>")
        sections.append(f'<section class="card">\n  <h2>{so}</h2>\n' + "\n".join(sf_blocks) + "\n</section>")

    page = f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Chiroptera SINEs</title>
<style>
:root {{ --fg:#222; --bg:#f5f6f8; --card:#fff; --accent:#4C72B0; --muted:#666; --border:#e2e2e2; }}
* {{ box-sizing: border-box; }}
body {{ font-family: -apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; background:var(--bg); color:var(--fg); margin:0; padding:0; }}
header {{ background: linear-gradient(135deg, #2c3e50, #4C72B0); color:#fff; padding:28px 36px; }}
header h1 {{ margin:0 0 6px 0; font-size:1.7rem; }}
header .sub {{ opacity:.85; font-size:.95rem; }}
main {{ max-width:1150px; margin:0 auto; padding:28px 24px; }}
section.card {{ background:var(--card); border:1px solid var(--border); border-radius:8px; padding:20px 24px; margin:18px 0; box-shadow:0 1px 3px rgba(0,0,0,0.05); }}
section.card h2 {{ margin-top:0; font-size:1.2rem; color:#2c3e50; border-bottom:1px solid var(--border); padding-bottom:8px; }}
h3.sf {{ font-size:1rem; color:#4C72B0; margin:18px 0 8px; }}
.tbl {{ border-collapse:collapse; width:100%; font-size:.86rem; margin:8px 0; }}
.tbl th,.tbl td {{ border:1px solid var(--border); padding:5px 9px; text-align:left; vertical-align:top; }}
.tbl th {{ background:#f0f3f7; }}
.tbl td.num {{ font-variant-numeric:tabular-nums; text-align:right; }}
.tbl tr.so td {{ background:#2c3e50; color:#fff; font-weight:600; }}
.tbl tr.sf td {{ background:#e8eef7; color:#2c3e50; font-weight:600; }}
.tbl .sub {{ font-size:.74rem; color:var(--muted); }}
.muted {{ color:var(--muted); }}
.meta {{ font-size:.82rem; color:var(--muted); margin-top:4px; }}
.tag {{ display:inline-block; border-radius:3px; padding:1px 6px; font-size:.78rem; font-weight:600; }}
.tag-ok {{ background:#d4edda; color:#155724; }}
.tag-partial {{ background:#fff3cd; color:#856404; }}
.tag-no {{ background:#f0f0f0; color:#777; }}
a.btn {{ display:inline-block; background:var(--accent); color:#fff; padding:4px 12px; border-radius:4px; font-size:.84rem; text-decoration:none; margin:2px 2px 2px 0; }}
a.btn:hover {{ background:#2c3e50; }}
a.btn.secondary {{ background:#6c757d; }}
a.btn.green {{ background:#28a745; }}
a {{ color:var(--accent); }}
.species-grid {{ display:grid; gap:16px; grid-template-columns: repeat(auto-fit, minmax(320px,1fr)); }}
.sp-card {{ border:1px solid var(--border); border-radius:7px; padding:14px 16px; background:#fff; }}
.sp-card h3 {{ margin:0 0 4px; font-size:1rem; color:#2c3e50; }}
.sp-card .sci {{ font-style:italic; color:var(--muted); font-size:.88rem; }}
.sp-card .actions {{ margin-top:10px; }}
</style>
</head><body>
<header>
  <h1>Chiroptera &mdash; bat SINEs</h1>
  <div class="sub">Rhin-1 and VES across bat families, one genome per family, grouped by clade.
  <a href="index.html" style="color:#fff;">&larr; Tal SINE main page</a></div>
</header>
<main>
<section class="card">
  <h2>Overview</h2>
  <p style="font-size:.9rem;">Clades and family order follow the time-calibrated family tree of Hao et al. 2023
  (<i>Integrative Zoology</i> 19:989&ndash;998, Fig. 3). Every genome was searched with SINEbase <b>Rhin-1</b> and
  <b>VES</b> as the only consensuses (IUPAC codes resolved to a base, so both keep full length); the earlier
  <i>Rhinolophus</i> runs used Rhin-1 alone. <b>Copies</b> = firm + soft assigned (firm count below);
  <b>sim</b> = median copy bitscore / consensus self-bitscore. <b>Prior evidence</b> is the SINE label on the
  annotated copy of that figure. <b>De novo</b> = AnnoSINE2 + SINEbase-fragment scan candidates, clustered by SubFam
  for manual review (run for hla and tbr first). {n_done} of {n_gen} genomes published.</p>
  <div style="overflow-x:auto;">
  <table class="tbl">
    <thead><tr><th>Family</th><th>Code</th><th>Species</th><th>Prior evidence</th>
      <th>Rhin-1 copies</th><th>Rhin-1 sim</th><th>VES copies</th><th>VES sim</th><th>De novo</th></tr></thead>
    <tbody>
{chr(10).join('      ' + t for t in table)}
    </tbody>
  </table>
  </div>
  <p style="margin-top:10px;"><a class="btn secondary" href="chiroptera/LOG.md">Analysis log</a>
  <a class="btn secondary" href="rhin.html">Rhin-1 peel alignments (<i>R. sinicus</i>)</a></p>
</section>
{chr(10).join(sections)}
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
