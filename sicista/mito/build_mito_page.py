#!/usr/bin/env python3
"""Build sicista_mito.html: Sicista mitogenome and mitochondrial-pseudogene (NUMT) alignments, opened in ViewAlign.

Run from the Tal repo root:  python sicista/mito/build_mito_page.py
The <style> block and the openMSA() script are copied from sicista.html at build time, so the page follows the
other Tal pages. Alignment files are in sicista/mito/alignments/ (built on therioserver, ~/Sicista2026/mito/aln/).
"""
import os, re, html

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ref = open(os.path.join(ROOT, 'sicista.html'), encoding='utf-8').read()
css = re.search(r'<style>.*?</style>', ref, re.S).group(0)
# Same RAW/MSA constants as sicista.html, plus an optional BED (gene track, ViewAlign ?bed=)
js = '''<script>
const RAW = 'https://raw.githubusercontent.com/Toki-bio/Tal/main/';
const MSA = 'https://toki-bio.github.io/MSA-viewer/';
function openMSA(relPath, title, bed) {
  let u = MSA + '?url=' + encodeURIComponent(RAW + relPath) + '&title=' + encodeURIComponent(title);
  if (bed) u += '&bed=' + encodeURIComponent(RAW + bed);
  window.open(u, '_blank');
}
</script>'''

A = 'sicista/mito/alignments/'
B = 'sicista/mito/annotations/'
BED_MITO = B + 'sicista_mitogenomes.NC_069019_genes.bed'      # GenBank annotation of NC_069019.1, on that row
BED_TRI = B + 'strizona_numts.OZ418355_genes_lifted.bed'       # the same, lifted onto OZ418355.1 through units/six.aln
BED_SB1 = B + 'sb1_numts.ptg633_unit2_genes_lifted.bed'        # the same, lifted onto Sb1 ptg633_unit2

def msa(fn, title, label, cls='btn', bed=None):
    return ('<a class="%s" href="javascript:void(0)" onclick="openMSA(\'%s%s\',\'%s\'%s)">%s</a>'
            % (cls, A, fn, html.escape(title, quote=True).replace("'", ''), (",'%s'" % bed) if bed else '', label))

def dl(fn, label):
    return '<a class="btn secondary" href="%s%s">%s</a>' % (A, fn, label)

def dl2(path, label):
    return '<a class="btn secondary" href="%s">%s</a>' % (path, label)

def card(h3, tag, sci, small, buttons):
    return ('    <div class="sp-card">\n      <h3>%s %s</h3>\n      <div class="sci">%s</div>\n'
            '      <div style="font-size:.82rem;color:var(--muted);margin-top:4px;">%s</div>\n'
            '      <div class="actions">\n        %s\n      </div>\n    </div>\n'
            % (h3, tag, sci, small, '\n        '.join(buttons)))

TAG = '<span class="tag tag-partial">unpublished, for review</span>'

mito_cards = [
    card('All Sicista mitogenomes', TAG,
         '61 GenBank mitochondrial records of five species plus the Sb1 mtDNA',
         '62 rows &middot; 16,808 columns &middot; MAFFT, reverse-complemented rows are prefixed <code>_R_</code>',
         [msa('sicista_mitogenomes.aln.fa', 'Sicista mitogenomes: 61 GenBank records + Sb1 own mtDNA', 'Open in MSA', bed=BED_MITO),
          dl('sicista_mitogenomes.aln.fa', 'FASTA'), dl2(BED_MITO, 'genes BED')]),
    card('Control region', TAG,
         'Region after tRNA-Pro to the end of each record, realigned on its own (L-INS-i)',
         '62 rows &middot; 1,294 columns &middot; the 6-bp ATACGC microsatellite is in Sb1 (x52), <i>S. trizona</i> and <i>S. caudata</i> (x56); not in the published <i>S. betulina</i> records',
         [msa('sicista_control_region.aln.fa', 'Sicista control regions (ATACGC repeat)', 'Open in MSA'),
          dl('sicista_control_region.aln.fa', 'FASTA')]),
    card('Sb1 mt-like contigs', TAG,
         'S. betulina ZBS2026-Sb1, hifiasm primary assembly: every segment of the 83 mt-like contigs aligned to Sb1\'s own mtDNA',
         '161 rows (own mtDNA + 160 segments) &middot; 151 segments are &ge;99.9% identical to the own mtDNA &middot; 2 outliers (84% and 86%) are not mtDNA proper &middot; <code>ptg000633l</code> carries five copies',
         [msa('sb1_mtlike_contigs.aln.fa', 'Sb1 mt-like contigs vs own mtDNA', 'Open in MSA', bed=BED_SB1),
          dl('sb1_mtlike_contigs.aln.fa', 'FASTA')]),
]

numt_cards = [
    card('S. trizona', TAG,
         'Sicista trizona &mdash; mSicTri.1, GCA_982266845.1, nuclear copies of OZ418355.1',
         '57 of 141 loci (aligned &ge;500 bp, or &ge;95% identity and &ge;100 bp) &middot; backbone row = the mtDNA, fragments clipped to its coordinates &middot; youngest: ND2 844 bp, 99.2% (OZ418345.1:183,346,458-183,347,301)',
         [msa('strizona_numts.aln.fa', 'S. trizona nuclear mt copies on OZ418355.1', 'Open in MSA', bed=BED_TRI),
          dl('strizona_numts.aln.fa', 'FASTA'), dl2(BED_TRI, 'genes BED')]),
    card('Sb1 nuclear copies, lineage-B family', TAG,
         'Sicista betulina Sb1 &mdash; 30 of the &ge;82 assembled full-length copies (16 kb), two young own-lineage copies, three GenBank references',
         '36 rows &middot; copies are 99.5-99.97% identical to each other, 99.3% to NC_069019 (lineage B) and ~94% to Sb1\'s own mtDNA (lineage A) &middot; about 200 haploid copies by read depth',
         [msa('sb1_nuclear_lineageB_family.aln.fa', 'Sb1 nuclear lineage-B mt family on own mtDNA', 'Open in MSA', bed=BED_SB1),
          dl('sb1_nuclear_lineageB_family.aln.fa', 'FASTA'), dl2(BED_SB1, 'genes BED')]),
    card('Sb1 nuclear copies, other loci', TAG,
         'Sicista betulina Sb1 &mdash; the 84 largest non-family loci (&ge;1 kb, or &ge;95% identity and &ge;100 bp)',
         '85 rows &middot; names carry contig, coordinates, percent identity and the mtDNA interval covered',
         [msa('sb1_nuclear_other_numts.aln.fa', 'Sb1 other nuclear mt copies on own mtDNA', 'Open in MSA', bed=BED_SB1),
          dl('sb1_nuclear_other_numts.aln.fa', 'FASTA')]),
]

body = '''<header>
  <h1>Sicista mitogenomes and mitochondrial pseudogenes</h1>
  <div class="sub">Alignments of mtDNA and its nuclear copies (NUMTs) in <i>Sicista betulina</i> Sb1 and <i>S. trizona</i>, opened in ViewAlign. <a href="sicista.html" style="color:#fff;">&larr; Sicista SINE page</a> &middot; <a href="index.html" style="color:#fff;">Tal SINE main page</a></div>
</header>
<main>
<section class="card">
  <h2>Mitogenomes</h2>
  <p style="font-size:.9rem;">Standard 37-gene mammalian arrangement in all five species; the control region is what differs (938 bp in published
  <i>S. betulina</i>, 982 in <i>S. strandi</i>, 1,044 in <i>S. concolor</i>, 1,206 in <i>S. caudata</i>, 1,217 in <i>S. trizona</i>, 1,199 in Sb1).
  The Sb1 mtDNA was found as repeated copies in the hifiasm primary assembly (no single mt contig) and is 99.2% identical to GenBank MZ570955 (one of two published <i>S. betulina</i> lineages).</p>
  <div class="species-grid">
''' + ''.join(mito_cards) + '''  </div>
</section>
<section class="card">
  <h2>Mitochondrial pseudogenes in the nuclear genome</h2>
  <p style="font-size:.9rem;">Nuclear copies were found by blastn of the species' mtDNA against the nuclear contigs (mt-like contigs excluded), loci merged within 2 kb.
  Each alignment has the mtDNA as its first row; the other rows are the nuclear copies placed on its coordinates, so a column is one mtDNA position.
  Insertions relative to the mtDNA are not shown.</p>
  <div class="species-grid">
''' + ''.join(numt_cards) + '''  </div>
  <table class="tbl" style="margin-top:14px;">
    <thead><tr><th>Genome</th><th>Nuclear mt loci</th><th>mt covered</th><th>Largest / youngest</th></tr></thead>
    <tbody>
      <tr><td><i>S. trizona</i></td><td class="num">141</td><td class="num">95.5%</td><td>4 copies of ~12 kb on OZ418342.1 at ~76%; ND2 844 bp at 99.2%</td></tr>
      <tr><td><i>S. betulina</i> Sb1</td><td class="num">305</td><td class="num">100%</td><td>&ge;82 full-length copies on 43 contigs (lineage B, ~94% to own mtDNA); 16.3 kb at 99.0% on ptg000319l; 12 kb at 99.8% on ptg000880l</td></tr>
    </tbody>
  </table>
</section>
<section class="card">
  <h2>Reading the alignments</h2>
  <ul style="font-size:.9rem;">
    <li>Row names in the NUMT alignments are <code>species_NUMT_contig:start-end_percent-identity[_mtstart-mtend]</code>; the percent is blastn identity over the merged locus.</li>
    <li>A leading <code>_R_</code> means MAFFT reverse-complemented that row.</li>
    <li>Mitogenomes and control regions are full records; Sb1 own mtDNA is the unit of <code>ptg000633l</code> (5 tandem copies, 0-2 differences).</li>
    <li>The gene track above the alignments is the GenBank annotation of NC_069019.1 (RefSeq <i>S. betulina</i>: 13 CDS, 22 tRNA, 2 rRNA, O<sub>L</sub>, control region),
      placed on that row and mapped through its gaps; for the NUMT alignments it was lifted onto OZ418355.1 and Sb1 <code>ptg633_unit2</code> through the six-mitogenome alignment
      (<code>gb2bed.py</code>, <code>lift_bed.py</code> in <code>sicista/mito/</code>). Hover a gene for its coordinates, click it to select its columns, hide it under Annotation in the viewer's settings.</li>
    <li>Sb1 alignments come from ONT reads (104 Gb, 35x) assembled with hifiasm; the Sb1 assembly is not in a public archive yet.</li>
  </ul>
  <p style="margin-top:14px;"><a class="btn secondary" href="sicista/mito/LOG.md">Analysis log</a>
  <a class="btn secondary" href="sicista.html">Sicista SINE page</a></p>
</section>
</main>
'''

page = ('<!doctype html>\n<html lang="en"><head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width,initial-scale=1">\n'
        '<title>Sicista mitogenomes and pseudogenes</title>\n' + css + '\n</head><body>\n' + body + js + '\n</body></html>\n')
open(os.path.join(ROOT, 'sicista_mito.html'), 'w', encoding='utf-8', newline='\n').write(page)
print('wrote sicista_mito.html', len(page))
