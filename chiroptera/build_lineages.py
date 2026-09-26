#!/usr/bin/env python3
"""Rebuild the "one genome per family" section of chiroptera.html.

Inputs (all in the Tal repo):
  chiroptera/genomes.tsv             code, species, family, accession, assembly, lineage, prior_evidence
                                     (rows in the order of Hao et al. 2023 Fig. 3, top to bottom)
  <code>/summary.by_subfam.tsv       SINEderella step3 summary, once the run is published
  <code>/alignments/*.aln.fa         published alignments (subfam input, de novo candidates)

The section sits between the LINEAGES START/END markers and is replaced wholesale,
so run this after every publish instead of editing the HTML by hand.
"""
import csv
import glob
import html
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAGE = os.path.join(ROOT, "chiroptera.html")
START, END = "<!-- LINEAGES START -->", "<!-- LINEAGES END -->"
SINES = ["Rhin-1", "VES"]


def summary(code):
    p = os.path.join(ROOT, code, "summary.by_subfam.tsv")
    if not os.path.exists(p):
        return None
    out = {}
    with open(p, encoding="utf8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["subfam"]] = r
    return out


def fmt(n):
    return f"{int(n):,}"


def msa_btn(path, title, label, cls="secondary"):
    return (f'<a class="btn {cls}" href="javascript:void(0)" '
            f"onclick=\"openMSA('{path}','{html.escape(title, quote=True)}')\">{label}</a>")


def main():
    with open(os.path.join(ROOT, "chiroptera", "genomes.tsv"), encoding="utf8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    trs, cards = [], []
    for g in rows:
        c, sp, fam = g["code"], g["species"], g["family"]
        s = summary(c)
        cells = []
        for sine in SINES:
            r = (s or {}).get(sine)
            if s is None:
                cells += ['<td class="num" colspan="2" style="color:var(--muted)">running</td>']
            elif r is None:
                cells += ['<td class="num">0</td><td class="num">&ndash;</td>']
            else:
                cells += [f'<td class="num">{fmt(r["total_assigned"])}'
                          f'<br><span style="font-size:.75rem;color:var(--muted)">{fmt(r["firm_assigned"])} firm</span></td>',
                          f'<td class="num">{float(r["sim_median"]):.2f}</td>']
        link = f'<a href="{c}/report.html">{c}</a>' if s else c
        trs.append(f'<tr><td>{g["lineage"]}</td><td>{fam}</td><td>{link}</td><td><i>{sp}</i></td>'
                   f'<td>{html.escape(g["prior_evidence"]) or "&ndash;"}</td>{"".join(cells)}</tr>')

        acts = []
        if s:
            acts.append(f'<a class="btn" href="{c}/report.html"><i>{sp}</i> &mdash; report</a>')
            sub = f"{c}/alignments/{c}_subfam_input_30k.aln.fa"
            if os.path.exists(os.path.join(ROOT, sub)):
                n = sum(1 for l in open(os.path.join(ROOT, sub), encoding="utf8") if l.startswith(">"))
                acts.append(msa_btn(sub, f"{c} ({sp}) SubFam input, {n - 2} chunk consensi + Rhin-1 + VES",
                                    f"SubFam input ({n - 2} chunks)"))
            for dn in sorted(glob.glob(os.path.join(ROOT, c, "alignments", f"{c}_denovo_candidates_*chunks.aln.fa"))):
                k = re.search(r"_(\d+)chunks", dn).group(1)
                acts.append(msa_btn(f"{c}/alignments/{os.path.basename(dn)}",
                                    f"{c} de novo candidates (AnnoSINE2 + fragment scan), {k} chunk consensi - for manual review",
                                    f"De novo candidates ({k} chunks)", "green"))
        tag = ('<span class="tag tag-ok">published</span>' if s
               else '<span class="tag tag-no">running</span>')
        counts = ""
        if s:
            parts = [f'{sine} {fmt(s[sine]["total_assigned"])}' for sine in SINES if sine in s]
            counts = (f'<div style="font-size:.82rem;color:var(--muted);margin-top:4px;">'
                      f'{" &middot; ".join(parts)} copies assigned</div>')
        cards.append(
            f'    <div class="sp-card">\n'
            f'      <h3>{c} {tag}</h3>\n'
            f'      <div class="sci">{sp} &mdash; {fam}</div>\n'
            f'      <div style="font-size:.82rem;color:var(--muted);margin-top:4px;">'
            f'<a href="https://www.ncbi.nlm.nih.gov/datasets/genome/{g["accession"]}/" target="_blank">{g["accession"]}</a> '
            f'{g["assembly"]}</div>\n'
            f'      {counts}\n'
            f'      <div class="actions">{" ".join(acts)}</div>\n'
            f'    </div>')

    section = f"""{START}
<section class="card">
  <h2>One genome per family &mdash; Rhin-1 and VES</h2>
  <p style="font-size:.9rem;">One representative genome for every family tip on the Hao et al. (2023) family tree,
  except Rhinolophidae (cards above) and Vespertilionidae (left out by request), searched with SINEbase Rhin-1 and VES as the only
  consensuses (IUPAC codes resolved to a base, so both keep full length). Families are in the tree's top-to-bottom order.
  <b>Copies</b> = firm + soft assigned (firm count below); <b>sim</b> = median copy bitscore / consensus self-bitscore.
  <b>Prior evidence</b> is the SINE label on the annotated copy of the figure (<code>Hao_24_IZ fig0003</code>). De novo candidates (AnnoSINE2 + fragment scan)
  are run for hla and tbr first.</p>
  <div style="overflow-x:auto;">
  <table class="tbl">
    <thead><tr><th>Lineage</th><th>Family</th><th>Code</th><th>Species</th><th>Prior evidence</th>
      <th>Rhin-1 copies</th><th>Rhin-1 sim</th><th>VES copies</th><th>VES sim</th></tr></thead>
    <tbody>
    {chr(10).join('    ' + t for t in trs)}
    </tbody>
  </table>
  </div>
  <div class="species-grid" style="margin-top:14px;">
{chr(10).join(cards)}
  </div>
</section>
{END}"""

    page = open(PAGE, encoding="utf8").read()
    if START in page:
        page = re.sub(re.escape(START) + ".*?" + re.escape(END), lambda m: section, page, flags=re.S)
    else:
        page = page.replace("</main>", section + "\n</main>", 1)
    open(PAGE, "w", encoding="utf8", newline="\n").write(page)
    done = sum(1 for g in rows if summary(g["code"]))
    print(f"{done}/{len(rows)} published")


if __name__ == "__main__":
    sys.exit(main())
