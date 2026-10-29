#!/usr/bin/env python3
"""Build chiroptera/orth_rsi-rfe/index.html: the rsi-rfe SINE_orth_loc flank-setting comparison, one line per locus
that differs between runs, each run's alignment opened in the SINE loci browser (sine_loci_browser.html, copied from
the SINE_orth_loc repo, which loads bundles named in its address).

Data (copied beside the page): data/loci.tsv and data/<category>/<run>/aln_rsi-rfe_<PM|MP|SINE|rejected>.aln.gz from
SINE_orth_loc tests/compare_real.py --out (therioserver:/home/toki/sol_bank/view_rsi-rfe); summary.txt = its stdout.
data/by_reference/<dir>/aln_rsi-rfe_<SINE|PM|MP>.aln.gz + index.tsv = orth_split.py output of the current run
(therioserver:/home/toki/sol_bank/rsi-rfe_bwa300both_v2/by_reference), linked from the results table.
Style: the Tal main page's <style>, as for chiroptera.html."""
import csv
import html
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_lineages import ROOT, tal_style

OUT = os.path.join(ROOT, "chiroptera", "orth_rsi-rfe")
RUNS = [("rsi-rfe", "bwa, 300 bp, one flank (baseline)"),
        ("rsi-rfe_bwa300both", "bwa, 300 bp, both flanks"),
        ("rsi-rfe_bwa300both_v2", "bwa, 300 bp, both flanks, junction from the annotated copy (41ccc90)"),
        ("rsi-rfe_bwa500", "bwa, 500 bp, one flank"),
        ("rsi-rfe_bwa500both", "bwa, 500 bp, both flanks"),
        ("rsi-rfe_mm2k11_500both", "minimap2 -x sr -k 11 -w 3, 500 bp, both flanks")]
CATS = [("disagree", "a run's row differs from the copy +-10 kb alignment check"),
        ("partner", "runs pair the locus with different loci in the other genome"),
        ("status", "runs give different calls (SINE / PM / MP)")] + \
       [(f"only_{r}", f"main row in {r} only") for r, _ in RUNS]


def summary_rows():
    t = open(os.path.join(OUT, "data", "summary.txt"), encoding="utf-8-sig").read().replace("\r", "")
    t = t.strip().split("\n\n")[0].splitlines()
    rows = [l.split("\t") for l in t if "\t" in l][1:]
    desc = dict(RUNS)
    return "\n".join(
        f'      <tr><td><strong>{r[0]}</strong><br><span style="color:var(--muted);font-size:.82rem;">'
        f'{html.escape(desc.get(r[0], ""))}</span></td>' + "".join(f'<td class="num">{int(x):,}</td>' for x in r[1:]) +
        "</tr>" for r in rows)


MIN_LOCI = 20   # references with fewer loci are summed into one line


def results_section():
    """Results of the current run (data/summary_rsi-rfe.tsv from orth_summary.py, data/run_facts.tsv)."""
    facts = dict(l.rstrip("\n").split("\t", 1) for l in open(os.path.join(OUT, "data", "run_facts.tsv"), encoding="utf8") if "\t" in l)
    rows = list(csv.DictReader(open(os.path.join(OUT, "data", "summary_rsi-rfe.tsv"), encoding="utf8"), delimiter="\t"))
    total = [r for r in rows if r["reference"] == "total"][0]
    big = [r for r in rows if r["reference"] != "total" and int(r["loci"]) >= MIN_LOCI]
    small = [r for r in rows if r["reference"] != "total" and int(r["loci"]) < MIN_LOCI]
    cols = ["copies_rsi", "copies_rfe", "loci", "SINE", "PM_rsi_only", "MP_rfe_only"]
    # by_reference/<dir>/aln_rsi-rfe_<status>.aln.gz from orth_split.py: the alignments behind the SINE / PM / MP counts
    idx = list(csv.DictReader(open(os.path.join(OUT, "data", "by_reference", "index.tsv"), encoding="utf8"), delimiter="	"))
    dir_of = {i["reference"]: i["dir"] for i in idx}
    status_of = {"SINE": "SINE", "PM_rsi_only": "PM", "MP_rfe_only": "MP"}

    def view(dirs, status, title):
        q = [f"bundle=data/by_reference/{d}/aln_rsi-rfe_{status}.aln.gz" for d in dirs]
        return "sine_loci_browser.html?" + "&".join(q + ["title=" + html.escape(title)])

    def td(r, c, dirs=None):
        v = r.get(c, "")
        if not v.lstrip("-").isdigit():
            return f'<td class="num">{html.escape(v)}</td>'
        if c in status_of and dirs and int(v):
            return f'<td class="num"><a href="{view(dirs, status_of[c], r["reference"] + " " + status_of[c])}" target="_blank">{int(v):,}</a></td>'
        return f'<td class="num">{int(v):,}</td>'

    def tr(r, name=None, strong=False, dirs=None):
        n = html.escape(name or r["reference"])
        n = f"<strong>{n}</strong>" if strong else n
        if dirs is None:
            dirs = [dir_of[r["reference"]]] if r["reference"] in dir_of else None
        return (f"      <tr><td>{n}</td>" + "".join(td(r, c, dirs) for c in cols) +
                f'<td class="num">{html.escape(r["shared_pct"])}</td>' +
                "".join(f'<td class="num">{html.escape(r[c])}</td>' for c in
                        ("identity_to_reference", "identity_between_copies", "identity_left_flanks")) + "</tr>")
    body = [tr(r) for r in big]
    if small:
        agg = {c: str(sum(int(r[c]) for r in small)) for c in cols}
        agg["shared_pct"] = f'{100 * int(agg["SINE"]) / int(agg["loci"]):.1f}' if int(agg["loci"]) else ""
        agg.update(identity_to_reference="", identity_between_copies="", identity_left_flanks="", reference="other")
        body.append(tr(agg, f"other ({len(small)} references with fewer than {MIN_LOCI} loci, mostly partial spans)",
                       dirs=["other"] if "other" in dir_of.values() else None))
    body.append(tr(total, "total", strong=True, dirs=sorted(set(dir_of.values()))))
    return f"""<section class="card">
  <h2>Results &mdash; {html.escape(facts["run"].split(" (")[0])}</h2>
  <p style="margin:0 0 10px;font-size:.9rem;color:var(--muted);">{html.escape(facts["run"])}. Bank: {html.escape(facts["bank"])}.
  Loci more than 300 bp from the next copy, whose flanks were mapped: rsi {int(facts["unique_loci_rsi"]):,}, rfe {int(facts["unique_loci_rfe"]):,}.
  Alignments rejected by ComPair: {html.escape(facts["rejected"])}; {int(facts["widened_references"]):,} clusters got a longer (container) reference;
  {int(facts["superseded"]):,} flank rows inside compound loci were superseded by the nested-locus analysis.
  Copies = SINEderella copies of that label per genome; loci = rows of the orth table (C flank-mapped, M rescued with both flanks, N from nested loci);
  SINE = in both genomes, PM = <i>R. sinicus</i> only, MP = <i>R. ferrumequinum</i> only; shared % = SINE / loci.
  Identities (%) are ComPair's: copy vs its reference consensus (divergence = 100 &minus; identity), copy vs copy at shared loci, and the left flanks;
  median and 10&ndash;90 % range. E1&ndash;E4, HE1 are the elements of nested loci (sine_nest naming), which carry no identities.
  <b>The Shared, rsi-only and rfe-only counts are links</b>: each opens that subfamily's alignments of that status alone in the
  SINE loci browser (data/by_reference/ from orth_split.py; nested-locus rows and the superseded flank rows have no alignment there).
  Inside the browser the subfamily list, the SINE / PM / MP / rejected chips and the name box narrow the loaded alignments further;
  <a href="{html.escape(view(sorted(set(dir_of.values())), "SINE", "rsi-rfe all") + "&bundle=" + "&bundle=".join(f"data/by_reference/{d}/aln_rsi-rfe_{t}.aln.gz" for t in ("PM", "MP") for d in sorted(set(dir_of.values()))))}" target="_blank">all 53,193 alignments at once</a> (30 MB) if you want to filter freely.</p>
  <div style="overflow-x:auto;">
  <table class="tbl">
    <thead><tr><th>Subfamily (reference)</th><th>Copies rsi</th><th>Copies rfe</th><th>Loci</th><th>Shared</th><th>rsi only</th><th>rfe only</th><th>Shared %</th><th>Identity to reference</th><th>Copy vs copy</th><th>Left flanks</th></tr></thead>
    <tbody>
{chr(10).join(body)}
    </tbody>
  </table>
  </div>
  <p style="margin:8px 0 0;"><a class="btn secondary" href="data/summary_rsi-rfe.tsv">summary_rsi-rfe.tsv (all {len(rows) - 1} references)</a></p>
</section>

"""


def main():
    runs_js = "[" + ",".join(f'"{r}"' for r, _ in RUNS) + "]"
    cat_btns = "\n    ".join(
        f'<a class="btn secondary" href="#" data-cat="{c}" title="{html.escape(d)}">{c} <span id="n-{c}"></span></a>'
        for c, d in CATS)
    cat_list = "\n".join(f"    <li><b>{c}</b>: {html.escape(d)}</li>" for c, d in CATS)
    page = f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Orthologous SINE loci rsi-rfe</title>
{tal_style()}
</head><body>
<header>
  <h1>Orthologous SINE loci &mdash; <i>R. sinicus</i> (rsi) vs <i>R. ferrumequinum</i> (rfe)</h1>
  <div class="sub">SINE_orth_loc in bank mode (SINEderella loci, bank21_nor1), five flank settings compared, 2026-10-28/29, plus the best one rerun with the insertion-junction fix.
  <a href="../../chiroptera.html" style="color:#fff;">&larr; Chiroptera page</a></div>
</header>
<main>

{results_section()}<section class="card">
  <h2>Runs</h2>
  <table class="tbl">
    <thead><tr><th>Setting</th><th>Rows</th><th>SINE</th><th>PM</th><th>MP</th><th>Agree</th><th>Disagree</th><th>Untested</th></tr></thead>
    <tbody>
{summary_rows()}
    </tbody>
  </table>
  <p style="color:var(--muted);font-size:.82rem;margin-top:8px;">
    Main flank-mapping rows only (clusters C&hellip;); rescue (M) and nested (N) rows are left out, so a locus a run
    holds as an M row counts as missing from that run. Check: each input copy &plusmn;10 kb was aligned to the other genome
    (minimap2 -x asm20) and its junction lifted; agree = the row's partner is within 50 bp of the lifted junction.
    Disagree is not the same as wrong: at segmental duplications the check itself can pick the wrong copy. Before the
    junction fix (SINE_orth_loc 41ccc90) about 5 % of anchors lay more than 50 bp off the window junction, because a
    reference longer than its copy (partial copy, composite) can start anywhere in the left flank; the comparison then
    uses the window position. With the fix the junction is the annotated copy's own start: 1.6 % remain off, and 198
    more rows become testable.
  </p>
</section>

<section class="card">
  <h2>Loci that differ between runs</h2>
  <p style="margin:0 0 10px;font-size:.9rem;color:var(--muted);">Pick a group, then open a run's alignment for a locus.
  The alignment opens in the SINE loci browser in a new tab, one tab per run, so the same locus can be compared side by side.
  A name starting with ~ means that run has no main row at the locus; it is the alignment the run made there instead
  (an M rescue row, or one rejected, e.g. shortRF / badRF).</p>
  <ul style="font-size:.85rem;color:var(--muted);margin:0 0 10px 18px;">
{cat_list}
  </ul>
  <div class="actions">
    {cat_btns}
  </div>
  <p style="margin:10px 0;"><input id="q" type="text" placeholder="search subfamily, locus, coordinate or alignment name" style="width:340px;max-width:100%;padding:4px 6px;">
  <span id="shown" style="color:var(--muted);font-size:.85rem;margin-left:8px;"></span></p>
  <div style="overflow-x:auto;">
  <table class="tbl" id="loci">
    <thead><tr><th>Locus</th><th>Groups</th><th>rsi junction</th><th>Check</th><th>Subfamily (reference row)</th>{"".join(f"<th>{r}</th>" for r, _ in RUNS)}</tr></thead>
    <tbody></tbody>
  </table>
  </div>
  <p><a class="btn secondary" href="#" id="more">Show 200 more</a>
  <a class="btn secondary" href="data/loci.tsv">loci.tsv</a></p>
</section>

</main>
<script>
const RUNS = {runs_js};
const TYPES = ['PM', 'MP', 'SINE', 'rejected'];
let LOCI = [], cat = 'disagree', shown = 200;
const esc = s => s.replace(/[&<>"]/g, c => ({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}})[c]);
function viewer(group, run, names, title) {{
  const q = TYPES.map(t => 'bundle=' + encodeURIComponent(`data/${{group}}/${{run}}/aln_rsi-rfe_${{t}}.aln.gz`));
  names.forEach(n => q.push('name=' + encodeURIComponent(n)));
  q.push('title=' + encodeURIComponent(title));
  return 'sine_loci_browser.html?' + q.join('&');
}}
function cell(L, run, v) {{
  if (v === '-') return '<span style="color:var(--muted);">nothing at the locus</span>';
  let names, label;
  if (v.startsWith('~')) {{ names = v.slice(1).split(','); label = esc(v); }}
  else {{ const p = v.split(' '); names = [p[2]]; label = `<b>${{esc(p[0])}}</b> ${{esc(p[1])}}<br>${{esc(p[2])}} ${{esc(p[3] || '')}}`; }}
  return `<a href="${{viewer(cat, run, names, L[0] + ' ' + run)}}" target="_blank">${{label}}</a>`;
}}
function render() {{
  const q = document.getElementById('q').value.trim().toLowerCase();
  const rows = LOCI.filter(L => L[1].split(',').includes(cat) && (!q || L.join(' ').toLowerCase().includes(q)));
  document.querySelector('#loci tbody').innerHTML = rows.slice(0, shown).map(L =>
    `<tr><td><strong>${{L[0]}}</strong></td><td>${{esc(L[1]).replace(/,/g, '<br>')}}</td><td>${{esc(L[2])}}</td><td>${{esc(L[3])}}</td><td><b>${{esc(L[4])}}</b></td>` +
    RUNS.map((r, i) => `<td style="font-size:.82rem;">${{cell(L, r, L[5 + i])}}</td>`).join('') + '</tr>').join('');
  document.getElementById('shown').textContent = `${{Math.min(shown, rows.length)}} of ${{rows.length}} shown`;
  document.getElementById('more').style.display = rows.length > shown ? '' : 'none';
  document.querySelectorAll('[data-cat]').forEach(b => b.classList.toggle('secondary', b.dataset.cat !== cat));
}}
fetch('data/loci.tsv').then(r => r.text()).then(t => {{
  LOCI = t.trim().split('\\n').slice(1).map(l => l.split('\\t'));
  const n = {{}};
  LOCI.forEach(L => L[1].split(',').forEach(c => n[c] = (n[c] || 0) + 1));
  for (const c in n) {{ const e = document.getElementById('n-' + c); if (e) e.textContent = '(' + n[c].toLocaleString() + ')'; }}
  render();
}});
document.querySelectorAll('[data-cat]').forEach(b => b.addEventListener('click', e => {{
  e.preventDefault(); cat = b.dataset.cat; shown = 200; render(); }}));
document.getElementById('q').addEventListener('input', () => {{ shown = 200; render(); }});
document.getElementById('more').addEventListener('click', e => {{ e.preventDefault(); shown += 200; render(); }});
</script>
</body></html>
"""
    open(os.path.join(OUT, "index.html"), "w", encoding="utf8", newline="\n").write(page)
    print("wrote", os.path.join(OUT, "index.html"))


if __name__ == "__main__":
    sys.exit(main())
