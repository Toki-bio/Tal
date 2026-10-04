#!/usr/bin/env python3
"""Write sicista/orth/section.html (included by build_page.py) from sicista/orth/{summary.json,class_by_subfamily.tsv,sample_loci.tsv}.
Run from the Tal repo root:  python sicista/build_orth_section.py && python sicista/build_page.py"""
import os, json, re, html

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
O = os.path.join(ROOT, 'sicista', 'orth')
S = json.load(open(os.path.join(O, 'summary.json')))
tot = sum(S[c]['n'] for c in S)
def n(x): return '{:,}'.format(x)
def pct(x): return '%.1f%%' % (100.0 * x / tot)

hdr, *rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(O, 'class_by_subfamily.tsv'))]
col = {h: i for i, h in enumerate(hdr)}
NONE = '(no firm Dip locus overlapping)'
label = {'SINE': 'Dip in both species (orthologous, shared)', 'PM': 'Dip in <i>S. betulina</i> only (empty site in <i>S. trizona</i>)',
         'MP': 'Dip in <i>S. trizona</i> only (empty site in <i>S. betulina</i>)'}

out = []
out.append('<section class="card">\n  <h2>Orthologous Dip loci &mdash; <i>S. betulina</i> vs <i>S. trizona</i> (SINE_orth_loc)</h2>\n')
out.append('  <p style="font-size:.9rem;">SINE_orth_loc (GitHub version) on the firmly assigned Dip loci of the two whole-genome runs above '
           '(1,032,760 in <i>S. betulina</i>, 1,036,280 in <i>S. trizona</i>; element: Dip_a1). The 300&nbsp;bp flank on one side of every locus is mapped to the other genome with bwa mem; '
           'each orthologous pair is aligned together with the consensus and classified by ComPair. '
           'Only loci with a unique, mappable flank enter, so these are counts of <b>validated locus pairs</b>, not of copies. '
           'Without an outgroup the direction cannot be read: a locus present in one species only is an insertion in that lineage <i>or</i> a deletion in the other. '
           'The <i>Jaculus jaculus</i> genome is the intended outgroup.</p>\n')
out.append('  <table class="tbl">\n    <thead><tr><th>Class</th><th>Locus pairs</th><th>Share</th><th>Dip_a1</th><th>Dip_a2</th><th>DIP</th><th>Dip_b1 / b2</th><th>No firm locus</th><th>sim_ratio median (n)</th></tr></thead>\n    <tbody>\n')
for r in rows:
    c = r[col['class']]; d = S[c]; get = lambda k: int(r[col[k]])
    out.append('      <tr><td>%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s / %s</td><td class="num">%s</td><td class="num">%.3f (%s)</td></tr>\n'
               % (label[c], n(d['n']), pct(d['n']), n(get('Dip_a1')), n(get('Dip_a2')), n(get('DIP')), n(get('Dip_b1')), n(get('Dip_b2')), n(get(NONE)), d['sim_median'], n(d['sim_ratio_n'])))
out.append('    </tbody>\n  </table>\n')
out.append('  <p style="font-size:.85rem;">Subfamily = the firmly assigned Dip locus that overlaps the aligned element on the side that carries it (<i>S. betulina</i> for the first two classes, <i>S. trizona</i> for the third). '
           '&ldquo;No firm locus&rdquo;: the alignment shows the element but the earlier assignment did not call it firmly. '
           'sim_ratio is the similarity to the assigned consensus (about 1 young, below 0.5 diverged) and exists for only part of the loci (<code>sim_scores.tsv</code>), so the medians are a rough age proxy; '
           'the two genomes have different overall values (Dip_a1 median 0.685 in <i>S. betulina</i>, 0.653 in <i>S. trizona</i>), so compare within a genome, not across. '
           'No subfamily has been curated by hand yet.</p>\n')
out.append('  <p>\n'
           '    <a class="btn" href="sicista/orth/aln_sbe-str_PM.aln.gz">All &ldquo;betulina only&rdquo; alignments (48,110, .gz)</a>\n'
           '    <a class="btn" href="sicista/orth/aln_sbe-str_MP.aln.gz">All &ldquo;trizona only&rdquo; alignments (44,891, .gz)</a>\n'
           '    <a class="btn secondary" href="sicista/orth/PM_sample_200.zip">PM sample 200</a>\n'
           '    <a class="btn secondary" href="sicista/orth/MP_sample_200.zip">MP sample 200</a>\n'
           '    <a class="btn secondary" href="sicista/orth/SINE_sample_200.zip">Shared sample 200</a>\n'
           '    <a class="btn secondary" href="sicista/orth/orth_sbe-str.tsv.gz">All 474,135 locus pairs (.tsv.gz)</a>\n'
           '  </p>\n'
           '  <p style="font-size:.85rem;">A bundle is plain text: a <code>##FILE name</code> line, then that alignment in FASTA (consensus last), repeated; '
           '<code>aln_bundle.sh</code> in SINE_orth_loc lists and extracts single alignments. The bundle of the 381,134 shared loci (163&nbsp;MB) is over the size GitHub accepts, '
           'so only a seeded sample of 200 is here; every locus pair, with coordinates, is in the table. Samples: seed 42.</p>\n')

# sample loci with NCBI fetch of the trizona side
out.append('  <h3 style="margin:16px 0 4px;font-size:1rem;color:#2c3e50;">Sample loci, with the <i>S. trizona</i> sequence fetched from NCBI by coordinates</h3>\n'
           '  <p style="font-size:.85rem;">The <i>S. trizona</i> assembly (GCA_982266845.1) is public, so its side of any locus can be fetched from NCBI by accession and coordinates; '
           'the button does this in your browser (NCBI efetch, start+1 to end, strand applied; tested identical to the local genome on 5 of 5 loci). '
           'The <i>S. betulina</i> assembly is not in a public archive, so its sequences are only in the alignment files.</p>\n')
out.append('  <table class="tbl">\n    <thead><tr><th>Class</th><th>Alignment</th><th><i>S. betulina</i> locus</th><th><i>S. trizona</i> locus (NCBI)</th><th></th></tr></thead>\n    <tbody>\n')
seen = {}
for line in open(os.path.join(O, 'sample_loci.tsv')).read().splitlines()[1:]:
    c, aln, sb, st = line.split('\t')
    seen[c] = seen.get(c, 0) + 1
    if seen[c] > 8:
        continue
    m = re.match(r'str_(.+):(\d+)-(\d+)\(([+-])\)$', st)
    acc, s0, e, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
    url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&seq_start=%d&seq_stop=%d&strand=%d&rettype=fasta&retmode=text'
           % (acc, s0 + 1, e, 1 if strand == '+' else 2))
    out.append('      <tr><td>%s</td><td><code>%s</code></td><td><code>%s</code></td><td><a href="%s" target="_blank"><code>%s:%d-%d(%s)</code></a></td>'
               '<td><a class="btn secondary" href="javascript:void(0)" onclick="fetchNcbi(this,\'%s\')">fetch</a></td></tr>\n'
               % (c, html.escape(aln), html.escape(sb), url, acc, s0, e, strand, url))
out.append('    </tbody>\n  </table>\n</section>\n')
JS = r'''<script>
function fetchNcbi(btn, url) {
  const tr = btn.closest("tr");
  let row = tr.nextElementSibling;
  if (!row || !row.classList.contains("ncbi-row")) {
    row = document.createElement("tr"); row.className = "ncbi-row";
    const td = document.createElement("td"); td.colSpan = 5;
    const pre = document.createElement("pre");
    pre.style.cssText = "margin:0;font-size:.74rem;white-space:pre-wrap;word-break:break-all;";
    td.appendChild(pre); row.appendChild(td); tr.after(row);
  }
  const pre = row.querySelector("pre"); pre.textContent = "fetching from NCBI ...";
  fetch(url).then(r => { if (!r.ok) throw new Error("HTTP " + r.status); return r.text(); })
    .then(t => { const L = t.trim().split("\n"); const seq = L.slice(1).join("");
      pre.textContent = L[0] + "  |  " + seq.length + " bp\n" + seq.slice(0, 400) + (seq.length > 400 ? " ..." : ""); })
    .catch(e => { pre.textContent = "NCBI fetch failed: " + e.message + " (rate limit is about 3 requests per second)"; });
}
</script>
'''
out.append(JS)
open(os.path.join(O, 'section.html'), 'w', encoding='utf-8').write(''.join(out))
print('wrote sicista/orth/section.html')
