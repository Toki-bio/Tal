#!/usr/bin/env python3
"""Build chiroptera/orth_rsi-rfe/index.html: the rsi-rfe SINE_orth_loc flank-setting comparison, one line per locus
that differs between runs, each run's alignment opened in the SINE loci browser (sine_loci_browser.html, copied from
the SINE_orth_loc repo, which loads bundles named in its address).

Data (copied beside the page): data/loci.tsv and data/<category>/<run>/aln_rsi-rfe_<PM|MP|SINE|rejected>.aln.gz from
SINE_orth_loc tests/compare_real.py --out (therioserver:/home/toki/sol_bank/view_rsi-rfe); summary.txt = its stdout.
Style: the Tal main page's <style>, as for chiroptera.html."""
import html
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_lineages import ROOT, tal_style

OUT = os.path.join(ROOT, "chiroptera", "orth_rsi-rfe")
RUNS = [("rsi-rfe", "bwa, 300 bp, one flank (baseline)"),
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
  <div class="sub">SINE_orth_loc in bank mode (SINEderella loci, bank21_nor1), four flank settings compared, 2026-10-28.
  <a href="../../chiroptera.html" style="color:#fff;">&larr; Chiroptera page</a></div>
</header>
<main>

<section class="card">
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
    Disagree is not the same as wrong: at segmental duplications the check itself can pick the wrong copy. In about 5 % of
    rows the reported anchor lies more than 50 bp inside the flank; there the window position is used as the junction.
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
  <p style="margin:10px 0;"><input id="q" type="text" placeholder="search locus, coordinate or alignment name" style="width:340px;max-width:100%;padding:4px 6px;">
  <span id="shown" style="color:var(--muted);font-size:.85rem;margin-left:8px;"></span></p>
  <div style="overflow-x:auto;">
  <table class="tbl" id="loci">
    <thead><tr><th>Locus</th><th>Groups</th><th>rsi junction</th><th>Check</th>{"".join(f"<th>{r}</th>" for r, _ in RUNS)}</tr></thead>
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
    `<tr><td><strong>${{L[0]}}</strong></td><td>${{esc(L[1]).replace(/,/g, '<br>')}}</td><td>${{esc(L[2])}}</td><td>${{esc(L[3])}}</td>` +
    RUNS.map((r, i) => `<td style="font-size:.82rem;">${{cell(L, r, L[4 + i])}}</td>`).join('') + '</tr>').join('');
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
