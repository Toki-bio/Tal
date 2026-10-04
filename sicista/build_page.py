#!/usr/bin/env python3
"""Build sicista.html in the layout of the other Tal pages (rhin.html / index.html species cards).

Run from the Tal repo root:  python sicista/build_page.py
Reads sicista/<code>/summary.by_subfam.tsv and lists sicista/<code>/alignments/ to decide which links exist,
and keeps the two older sections (step-1 database arm on the Medaka assembly, hit counts) verbatim from sicista/_old_sections.html.
Sections for orthologous loci and satellites are added when sicista/orth/summary.json / sicista/sat/summary.json exist.
"""
import os, re, json, html

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S = os.path.join(ROOT, 'sicista')

SPECIES = [
    dict(code='sbe', sci='Sicista betulina', desc='northern birch mouse, ZBS2026-Sb1 &mdash; hifiasm primary assembly (1,144 contigs, 2.99 Gb)',
         acc='', denovo='sbet', minus='sbe'),
    dict(code='str', sci='Sicista trizona', desc='mSicTri.1 (chromosome level, 2.80 Gb)',
         acc='GCA_982266845.1', denovo='str', minus='str'),
]
ORDER = ['DIP', 'Dip_a1', 'Dip_a2', 'Dip_b1', 'Dip_b2']

def read_summary(code):
    p = os.path.join(S, code, 'summary.by_subfam.tsv')
    rows = []
    if os.path.exists(p):
        with open(p) as f:
            hdr = f.readline().rstrip('\n').split('\t')
            for line in f:
                v = line.rstrip('\n').split('\t')
                rows.append(dict(zip(hdr, v)))
    rows.sort(key=lambda r: ORDER.index(r['subfam']) if r['subfam'] in ORDER else 99)
    return rows

def has(code, name):
    return os.path.exists(os.path.join(S, code, 'alignments', name))

def msa(path, title, label, cls='btn secondary'):
    return ('<a class="%s" href="javascript:void(0)" onclick="openMSA(\'%s\',\'%s\')">%s</a>'
            % (cls, path, html.escape(title, quote=True).replace("'", ''), label))

def num(x):
    try: return '{:,}'.format(int(float(x)))
    except Exception: return x

def head():
    old = open(os.path.join(S, '_old_sections.html'), encoding='utf-8').read()
    css = re.search(r'<style>.*?</style>', old, re.S)
    css = css.group(0) if css else ''
    return css

def old_sections():
    t = open(os.path.join(S, '_old_sections.html'), encoding='utf-8').read()
    a = t.index('<!--OLD_A-->'); b = t.index('<!--/OLD_A-->')
    c = t.index('<!--OLD_D-->'); d = t.index('<!--/OLD_D-->')
    return t[a + len('<!--OLD_A-->'):b], t[c + len('<!--OLD_D-->'):d]

def build():
    sums = {sp['code']: read_summary(sp['code']) for sp in SPECIES}
    out = []
    out.append('<!doctype html>\n<html lang="en"><head>\n<meta charset="utf-8">\n'
               '<meta name="viewport" content="width=device-width,initial-scale=1">\n'
               '<title>Sicista SINEs</title>\n' + head() + '\n</head><body>\n')
    out.append('<header>\n  <h1>Sicista SINE &mdash; Sicista betulina and S. trizona</h1>\n'
               '  <div class="sub">The Dip SINE and de novo SINE candidates in two birch-mouse genomes (Dipodoidea). '
               '<a href="index.html" style="color:#fff;">&larr; Tal SINE main page</a></div>\n</header>\n<main>\n')

    # ---- species cards
    out.append('<section class="card">\n  <h2>Species</h2>\n  <div class="species-grid">\n')
    for sp in SPECIES:
        c = sp['code']; rows = sums[c]
        tot = sum(int(float(r['firm_assigned'])) for r in rows) if rows else 0
        totall = sum(int(float(r['total_assigned'])) for r in rows) if rows else 0
        acc = ('&middot; <a href="https://www.ncbi.nlm.nih.gov/datasets/genome/%s/" target="_blank">%s</a>' % (sp['acc'], sp['acc'])) if sp['acc'] else '&middot; assembly not yet in a public archive'
        rep = os.path.exists(os.path.join(S, c, 'report.html'))
        tag = '<span class="tag tag-ok">full report</span>' if rep else '<span class="tag tag-partial">pending manual review</span>'
        out.append('    <div class="sp-card">\n      <h3>%s %s</h3>\n      <div class="sci">%s &mdash; %s</div>\n' % (c, tag, sp['sci'], sp['desc']))
        if rows:
            out.append('      <div style="font-size:.82rem;color:var(--muted);margin-top:4px;">%d Dip-bank subfamilies &middot; %s firm assigned copies (%s with soft assignments) %s</div>\n'
                       % (len(rows), num(tot), num(totall), acc))
        out.append('      <div class="actions">\n')
        if rep:
            out.append('        <a class="btn" href="sicista/%s/report.html"><i>%s</i> &mdash; full report (Dip bank)</a>\n' % (c, sp['sci']))
        out.append('        ' + msa('sicista/alignments/dip_consensuses.aln.fa', 'Dip bank consensi (as searched)', 'Dip consensi in MSA') + '\n')
        for nm, lab, ttl in ((sp['denovo'] + '_denovo_30k.aln.fa', 'De novo, 10% subsample', sp['code'] + ' de novo (scan + AnnoSINE2, 10% of genome), 30k sample'),
                             (sp['minus'] + '_minus_denovo_30k.aln.fa', 'De novo, Dip-depleted genome', sp['code'] + ' de novo on the Dip-depleted genome, 30k sample')):
            if os.path.exists(os.path.join(S, 'alignments', nm)):
                out.append('        ' + msa('sicista/alignments/' + nm, ttl, lab) + '\n')
        out.append('      </div>\n    </div>\n')
    out.append('  </div>\n  <p style="margin-top:14px;"><a class="btn secondary" href="sicista/LOG.md">Analysis log</a> '
               '<a class="btn secondary" href="sicista/alignments/dip_consensuses.fa">Dip bank (5 consensuses, as searched)</a> '
               '<a class="btn secondary" href="sicista_mito.html">Mitogenomes and mt pseudogenes</a> '
               '<a class="btn secondary" href="sicista_qc.html">Sb1 assembly QC and the Illumina decision</a></p>\n</section>\n')

    # ---- Dip SINEderella, per-subfamily table
    out.append('<section class="card">\n  <h2>Dip SINE &mdash; SINEderella on the whole genome</h2>\n'
               '  <p style="font-size:.9rem;">SINEderella steps 1&ndash;4 and the report on the whole genome, searched with five consensuses: SINEbase <b>DIP</b> '
               'and the literature subfamilies <b>Dip_a1, Dip_a2, Dip_b1, Dip_b2</b> (as searched: IUPAC codes removed, so DIP is 206 bp, Dip_a1 224, Dip_b1 240). '
               'Assignment is by the 10-cycle unanimous vote; <i>firm</i> copies passed it, <i>soft</i> copies did not but have a best-scoring query. '
               'The subfamilies share most of their sequence, so most copies are also hit by the others (the leak and conflict columns in the full report). '
               'Alignments carry 50&nbsp;bp left and 70&nbsp;bp right flanks, the consensus is row&nbsp;1 and the sequence as searched is row&nbsp;2; '
               'a subfamily with fewer than 100 firm copies has its 100-copy plates filled with soft copies (marked <code>[soft]</code>). '
               'No subfamilies have been curated by hand yet.</p>\n')
    for sp in SPECIES:
        c = sp['code']; rows = sums[c]
        if not rows:
            continue
        out.append('  <h3 style="margin:16px 0 4px;font-size:1rem;color:#2c3e50;"><i>%s</i> (%s)</h3>\n' % (sp['sci'], c))
        out.append('  <table class="tbl">\n    <thead><tr><th>Subfamily</th><th>Firm</th><th>Soft</th><th>sim_ratio median</th><th>Top 100</th><th>100 random</th><th>SubFam</th></tr></thead>\n    <tbody>\n')
        for r in rows:
            sf = r['subfam']; cells = []
            for kind, lab in (('top100', 'top 100'), ('rand100', '100 random'), ('subfam', 'SubFam')):
                nm = '%s_%s_%s.aln.fa' % (c, sf, kind)
                cells.append(('<a href="javascript:void(0)" onclick="openMSA(\'sicista/%s/alignments/%s\',\'%s %s %s\')">%s</a>' % (c, nm, c, sf, lab, lab)) if has(c, nm) else '&ndash;')
            out.append('      <tr><td><code>%s</code></td><td class="num">%s</td><td class="num">%s</td><td class="num">%.2f</td><td>%s</td><td>%s</td><td>%s</td></tr>\n'
                       % (sf, num(r['firm_assigned']), num(r['soft_assigned']), float(r.get('sim_median', 0) or 0), cells[0], cells[1], cells[2]))
        out.append('    </tbody>\n  </table>\n')
        out.append('  <p style="font-size:.85rem;margin:4px 0 0;">SubFam alignments exist only for subfamilies with at least 400 copies. '
                   'Same files as the <a href="sicista/%s/report.html">full report</a>.</p>\n' % c)
    out.append('</section>\n')

    # ---- de novo
    out.append('<section class="card">\n  <h2>De novo SINE search</h2>\n'
               '  <p style="font-size:.9rem;">SINE-de-novo-genome-scan (50&nbsp;bp fragments of SINEbase and literature consensuses, ssearch36, identity &ge;65%, query coverage &ge;0.90, hits merged within 500&nbsp;bp) '
               'merged with AnnoSINE2 seeds where run, k-mer thinned (singletons dropped), a random 30,000-copy sample, and SubFam chunks of 50 (rows are chunk consensuses ordered by similarity). '
               'Two searches per genome:</p>\n'
               '  <ul style="font-size:.9rem;">\n'
               '    <li><b>10% subsample, with AnnoSINE2.</b> A seeded random 300 &times; 1&nbsp;Mb windows of each genome (the full-genome scan with Dip present was projected at over 20 h per genome). All AnnoSINE2 seeds are kept in the sample.</li>\n'
               '    <li><b>Dip-depleted whole genome, scan only.</b> Every Dip-bank hit (about 258&nbsp;Mb in <i>S. betulina</i>, 264&nbsp;Mb in <i>S. trizona</i>) is cut out first, pieces &ge;100&nbsp;bp are scanned; about 2 h per genome.</li>\n'
               '  </ul>\n  <table class="tbl">\n    <thead><tr><th>Genome</th><th>Search</th><th>Candidates</th><th>After thinning</th><th>Rows</th><th>Alignment</th></tr></thead>\n    <tbody>\n')
    DN = [('sbe', 'S. betulina', '10% subsample + AnnoSINE2', '279,697 + 1,090 seeds', '269,166', '621', 'sbet_denovo_30k.aln.fa'),
          ('sbe', 'S. betulina', 'Dip-depleted genome', '92,161', '71,856', '600', 'sbe_minus_denovo_30k.aln.fa'),
          ('str', 'S. trizona', '10% subsample + AnnoSINE2', '273,316 + 818 seeds', '262,623', '616', 'str_denovo_30k.aln.fa'),
          ('str', 'S. trizona', 'Dip-depleted genome', '98,540', '78,108', '600', 'str_minus_denovo_30k.aln.fa')]
    for c, g, s, cand, thin, rows, fn in DN:
        link = msa('sicista/alignments/' + fn, '%s %s, 30k sample' % (g, s), 'open in MSA', 'btn') if os.path.exists(os.path.join(S, 'alignments', fn)) else '&ndash;'
        out.append('      <tr><td><i>%s</i></td><td>%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td>%s</td></tr>\n' % (g, s, cand, thin, rows, link))
    out.append('    </tbody>\n  </table>\n  <p style="margin-top:10px;">' +
               ' '.join('<a class="btn secondary" href="sicista/alignments/%s">%s</a>' % (fn, lab) for fn, lab in
                        (('sbet_annosine_seeds.fa', 'AnnoSINE2 seeds, S. betulina (1,090)'), ('str_annosine_seeds.fa', 'AnnoSINE2 seeds, S. trizona (818)'))
                        if os.path.exists(os.path.join(S, 'alignments', fn))) + '</p>\n</section>\n')

    # ---- optional: orth, satellites
    for key, builder in (('orth', None), ('sat', None)):
        p = os.path.join(S, key, 'section.html')
        if os.path.exists(p):
            out.append(open(p, encoding='utf-8').read() + '\n')

    # ---- older material, verbatim
    a, d = old_sections()
    out.append('<section class="card">\n  <h2>Earlier run &mdash; SINEderella step 1 on the Medaka assembly (superseded by the primary assembly)</h2>\n'
               '  <p style="font-size:.9rem;">First pass, 2026-09-26, with the 58 SINEbase Mammalia consensuses on the Medaka assembly (8,737 scaffolds). Kept for the per-consensus hit counts; the runs above use the hifiasm primary assembly.</p>\n</section>\n')
    out.append(a + '\n' + d + '\n')

    out.append('</main>\n<script>\nconst RAW = \'https://raw.githubusercontent.com/Toki-bio/Tal/main/\';\n'
               'const MSA = \'https://toki-bio.github.io/MSA-viewer/\';\nfunction openMSA(relPath, title) {\n  const url = RAW + relPath;\n'
               '  window.open(MSA + \'?url=\' + encodeURIComponent(url) + \'&title=\' + encodeURIComponent(title), \'_blank\');\n}\n</script>\n</body></html>\n')
    open(os.path.join(ROOT, 'sicista.html'), 'w', encoding='utf-8').write(''.join(out))
    print('wrote sicista.html (%d bytes)' % sum(len(x) for x in out))

if __name__ == '__main__':
    build()
