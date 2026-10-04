#!/usr/bin/env python3
"""Build sicista_qc.html: S. betulina Sb1 assembly QC, the ONT-only tests and the decision on additional Illumina.

Run from the Tal repo root:  python sicista/qc/build_qc_page.py
The <style> block is copied from sicista.html at build time, so the page follows the other Tal pages.
Figures and the full notes are in sicista/qc/ (analysis on therioserver ~/Sicista2026/qc/, 2026-10-03/04).
"""
import os, re

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ref = open(os.path.join(ROOT, 'sicista.html'), encoding='utf-8').read()
css = re.search(r'<style>.*?</style>', ref, re.S).group(0)
css = css.replace('</style>', '.fig { max-width:100%; border:1px solid var(--border); border-radius:6px; }\n'
                  '.caption { font-size:.85rem; color:var(--muted); margin:6px 0 14px; }\n'
                  '.verdict li { margin-bottom:6px; }\n</style>')

Q = 'sicista/qc/'
TAG = '<span class="tag tag-partial">unpublished, for review</span>'

body = '''<header>
  <h1>Sicista betulina Sb1: assembly QC and the Illumina decision</h1>
  <div class="sub">Is additional Illumina sequencing needed for PSMC/MSMC, gene phylogeny and SINE insertion polymorphism? Tests on the ONT-only hifiasm assembly, 3&ndash;4 October 2026.
  <a href="sicista.html" style="color:#fff;">&larr; Sicista SINE page</a> &middot; <a href="sicista_mito.html" style="color:#fff;">Mitogenomes and mt pseudogenes</a> &middot; <a href="index.html" style="color:#fff;">Tal SINE main page</a></div>
</header>
<main>
<section class="card">
  <h2>Verdict ''' + TAG + '''</h2>
  <ul class="verdict" style="font-size:.95rem;">
    <li><b>Not needed on this individual</b> for PSMC, gene phylogeny or within-individual SINE indels. The ONT reads (104 Gb, 35&times;, R10.4.1-class) give a het-SNP false-positive rate of about 1% after an allele-fraction filter (measured on the hemizygous X, where a male has no true hets), two independent read halves reproduce 88&ndash;90% of each other's calls and the same PSMC curve, the 20 PSMC bootstraps are tight, and BUSCO is 99.2% complete.</li>
    <li><b>Needed for the population-level goals, whichever platform:</b> MSMC/MSMC2 beyond one genome (2&ndash;4 individuals), SINE insertion frequencies and population genetics need more individuals, not more reads of Sb1. Cheapest route: 15&ndash;20&times; PCR-free Illumina on 4&ndash;8 additional individuals, mapped to this assembly.</li>
    <li><b>Optional on Sb1:</b> 30&times; Illumina (80&ndash;90 Gb) as orthogonal truth for the base-level QV, the false/missed het rate and homopolymer-indel polishing. The ONT-only estimate of the assembly's accuracy is QV&nbsp;~43 (upper bound on the error), with indels in homopolymers dominating, which matters for gene models; self-polishing with the homozygous Clair3 calls is possible without Illumina.</li>
    <li><b>Not verified:</b> the sequencing chemistry and basecaller model (no tags in the read headers), the mutation rate and generation time of <i>Sicista</i> (PSMC time and Ne axes are illustrative), the sex (inferred male from X depth), and the recall of true heterozygotes (no orthogonal call set).</li>
  </ul>
</section>

<section class="card">
  <h2>Assembly QC (hifiasm primary, ONT reads)</h2>
  <table class="tbl">
    <thead><tr><th>Metric</th><th>Value</th></tr></thead>
    <tbody>
      <tr><td>Contigs / total length</td><td>1,144 / 2.992 Gb, no gaps</td></tr>
      <tr><td>N50 / L50; N75 / N90; largest</td><td>18.25 Mb / 48; 8.03 Mb / 1.93 Mb; 84.3 Mb</td></tr>
      <tr><td>Contigs &ge;10 Mb / &ge;1 Mb / &ge;100 kb</td><td>91 (2.11 Gb, 70.5%) / 292 (2.80 Gb, 93.7%) / 750 (2.98 Gb, 99.5%); 394 contigs &lt;100 kb hold 14.3 Mb</td></tr>
      <tr><td>GC</td><td>42.5% (contigs &ge;1 Mb: 36.5&ndash;54.3%, median 42.1%)</td></tr>
      <tr><td>BUSCO mammalia_odb10 (n = 9,226)</td><td>C 99.2% [S 97.1%, D 2.1%], F 0.2%, M 0.6% (the <i>S. trizona</i> reference: C 98.3%)</td></tr>
      <tr><td>Merqury with ONT k-mers</td><td>QV 58.3, completeness 93.4% (relative only: the k-mers come from the same reads)</td></tr>
      <tr><td>Telomeres (TTAGGG, &ge;60 bp in the terminal 10 kb)</td><td>none of the 620 contigs &ge;200 kb has both ends; 30 have one</td></tr>
      <tr><td>GenomeScope2 (ONT 21-mers)</td><td>haploid 2.69 Gb, heterozygosity 0.50&ndash;0.57% (unreliable with ONT k-mers); assembly is 11% larger</td></tr>
      <tr><td>Mitochondrial contamination</td><td>83 mt-like contigs (1.37 Mb) to remove from the nuclear set; 305 nuclear NUMT loci (<a href="sicista_mito.html">mitogenome page</a>)</td></tr>
      <tr><td>Read depth (minimap2 map-ont, all reads)</td><td>median 35.1&times; in 10-kb windows (5&ndash;95%: 15&ndash;46&times;); 75.5% of windows within 0.75&ndash;1.25&times; of the median</td></tr>
      <tr><td>Half-depth sequence</td><td>10.8% of bp at MAPQ &ge;10; 6% (179 Mb, 494 contigs) stays half-depth with MAPQ-0 reads counted, the rest is low-MAPQ repeat</td></tr>
      <tr><td>Sex chromosomes</td><td>81 half-depth contigs &ge;1 Mb have no &ge;90%-identity counterpart in the assembly (not unpurged haplotigs); 89.3% of their alignment to <i>S. trizona</i> goes to one chromosome (OZ418353.1, 148 Mb; 165&times; enriched) and 3.4% to OZ418354.1 (20 Mb): the sample is a male, X (+Y) about 175 Mb, excluded from autosomal PSMC</td></tr>
      <tr><td>Diploid-depth contigs</td><td>292 contigs, 2,671 Mb (89%); 212 contigs &ge;1 Mb = 2.64 Gb. Largest: ptg000057l 84.3, ptg000018l 75.0, ptg000003l 62.5, ptg000566l 58.8, ptg000116l 56.8 Mb</td></tr>
      <tr><td>Reads</td><td>13.5 M reads, 104.25 Gb (35&times;), N50 12.05 kb, mean Q15.8, 21.5% of reads above Q20; read headers carry no basecaller-model tags</td></tr>
    </tbody>
  </table>
</section>

<section class="card">
  <h2>ONT-only variant calling: how reliable are the heterozygotes?</h2>
  <p style="font-size:.9rem;">Clair3 (sup v5.0.0 model for R10.4.1), first on a test set of 5 autosomal contigs (337 Mb, 35&times;), then genome-wide.</p>
  <table class="tbl">
    <thead><tr><th>Test</th><th>Result</th></tr></thead>
    <tbody>
      <tr><td>Heterozygous SNP density, 35&times;</td><td>627,351 PASS het SNPs = 1,860 per Mb (0.19%); 83% have allele fraction 0.35&ndash;0.65 (median 0.46); Ti/Tv 1.64 (true mammalian SNPs: about 2.0&ndash;2.2, so some transversion-rich false calls, of the order of 15%)</td></tr>
      <tr><td>Two disjoint read halves (17.5&times; each)</td><td>659k het SNPs each; 87.8% of one half's calls are in the other without filtering; 99.3% of the shared set is recovered at 35&times;; 99.6% of 35&times; calls are in the union. Calls private to one half (80k each) have Ti/Tv 2.9, so a half overestimates &theta; by about 17%</td></tr>
      <tr><td>Hard filters</td><td>Fixed QUAL/AF thresholds change het density 7-fold (1,950 &rarr; 261 per Mb) and are depth-dependent: an allele-fraction filter is used instead of fixed QUAL across depths</td></tr>
      <tr><td>False-positive rate from the hemizygous X</td><td>On 94.8 Mb of X-class contigs (5 contigs with &gt;600 het/Mb excluded as PAR/duplicated): 163 het SNP/Mb unfiltered = 8.4% of the autosomal 17&times; rate; with AF 0.35&ndash;0.65: 14/Mb, about 1% of the autosomal filtered rate. False hets have median AF 0.29, true ones 0.46</td></tr>
      <tr><td>Assembly base accuracy, ONT-only estimate</td><td>Homozygous-alt calls (reads disagree with the assembly): 5,160 SNPs + 11,505 indels on 337 Mb = 49 per Mb = 4.9 &times; 10<sup>-5</sup>, QV about 43 (upper bound on the error); indels in homopolymers dominate</td></tr>
      <tr><td>Phasing (WhatsHap, ONT reads)</td><td>99.1% of het SNPs phased; blocks span most of each contig (ptg000566l, 58.8 Mb, in 2 blocks; block NG50 48.8 Mb): haplotype-resolved calls are feasible from ONT alone</td></tr>
    </tbody>
  </table>
</section>

<section class="card">
  <h2>PSMC from ONT alone</h2>
  <p style="font-size:.9rem;">Genome-wide: Clair3 on the 351 diploid-depth contigs (2.684 Gb; mt-like and X-class excluded), callable mask MAPQ &ge;10 and depth 18&ndash;62&times; (2.589 Gb callable), 100-bp bins, het-bin fraction 13.5% (AF 0.3&ndash;0.7) or 14.5% (PASS only).
  <code>psmc -N25 -t15 -r5 -p "4+25*2+4+6"</code>, 20 bootstraps on the AF-filtered set. &mu; = 5.4 &times; 10<sup>-9</sup> per generation and 1 generation = 1 year are illustrative, not verified for <i>Sicista</i>.</p>
  <img class="fig" src="''' + Q + '''sicista_psmc_ont_only_2026-10-04.png" alt="ONT-only PSMC of S. betulina Sb1, genome-wide, with 20 bootstraps">
  <div class="caption">Genome-wide ONT-only PSMC: the AF-filtered call set (red) with its 20 bootstraps (grey) and the unfiltered PASS set (blue). A plateau of 150&ndash;190k from before 2 My to about 100 ky, a decline to about 70k at 50&ndash;60 ky, a bump to about 130k at 25&ndash;40 ky, then a collapse to 10&ndash;30k by 15&ndash;20 ky ago. The two call sets agree beyond 60 ky and differ in the timing (up to about 10 ky) and amplitude of the last 40 ky, where PSMC has little power anyway; absolute Ne scales with &theta;, which shifts about 15% between the filters.</div>
  <table class="tbl">
    <thead><tr><th>Time (years ago, illustrative scale)</th><th>Ne, AF-filtered set (bootstrap 2.5&ndash;97.5%)</th></tr></thead>
    <tbody>
      <tr><td>2 My</td><td class="num">~190k (177&ndash;211k)</td></tr>
      <tr><td>1 My</td><td class="num">~148k (144&ndash;157k)</td></tr>
      <tr><td>300 ky</td><td class="num">~166k (160&ndash;168k)</td></tr>
      <tr><td>100 ky</td><td class="num">~117k (109&ndash;126k)</td></tr>
      <tr><td>50 ky</td><td class="num">~71k (69&ndash;85k)</td></tr>
      <tr><td>20 ky</td><td class="num">~80k (45&ndash;96k)</td></tr>
      <tr><td>10 ky</td><td class="num">~12k (9&ndash;12k)</td></tr>
      <tr><td>present-day N<sub>0</sub></td><td class="num">~48k (AF-filtered) / ~55k (PASS only)</td></tr>
    </tbody>
  </table>
  <img class="fig" src="''' + Q + '''sicista_psmc_subset_tests_2026-10-04.png" alt="PSMC on the 337-Mb test set: full depth and two independent read halves, with and without filters">
  <div class="caption">The robustness test on the 337-Mb subset (11% of the genome): the full 35&times; call set (black, red) and the two independent 17.5&times; read halves A and B (blue/green, orange/purple), each unfiltered (pass) or quality- and AF-filtered (q10af, q20af). The halves give nearly the same curve as each other beyond about 30 ky and agree on a 10&ndash;30 ky trough; the full set agrees with them beyond 50 ky and differs 2&ndash;5&times; below 40 ky, where its &theta; is 15% lower. The deep-time history is robust to depth and filtering; the recent part depends on the call set. The genome-wide curve above has the same deep-time Ne (~150k), so the subset was representative.</div>
</section>

<section class="card">
  <h2>SINE insertion polymorphism within Sb1</h2>
  <p style="font-size:.9rem;">Sniffles2 (2.8.1) on the ONT reads against the assembly, autosomal-class contigs (2.68 Gb): 70,853 PASS insertions and deletions &ge;50 bp; 46,423 of 100&ndash;700 bp with sequence, of which 25,313 (54.5%) match the Dip SINE bank (Dip_a2 14,437; Dip_a1 8,255; DIP 2,565): 13,820 deletions and 11,493 insertions relative to the assembly, about 9.4 per Mb, median length 244 bp. All are heterozygous, as expected within one individual. On the test subset, 89&ndash;90% of the 100&ndash;700 bp calls from one read half are found in the other and 98% of the half-calls are in the 35&times; set.
  These are within-individual heterozygous indels; allele frequencies need a population sample, which is the Illumina that is worth buying.</p>
</section>

<section class="card">
  <h2>Files and provenance</h2>
  <ul style="font-size:.9rem;">
    <li>Analysis on therioserver, <code>~/Sicista2026/qc/</code>: <code>sbet_ont.bam</code> (all reads on the primary assembly), <code>calls/</code> (test set), <code>calls_full/</code>, <code>calls_X/</code>, <code>sv_full/</code>, <code>final/</code> (genome-wide PSMC), <code>models/</code>. Assembly: <code>assembly/primary_eval/hifiasm_primary.p_ctg.fa</code>.</li>
    <li>Tools: minimap2 map-ont, Clair3 (r1041_e82_400bps_sup_v500), WhatsHap, psmc, Sniffles2 2.8.1, BUSCO (mammalia_odb10), Merqury, GenomeScope2.</li>
    <li>Full notes with every number: <a href="''' + Q + '''QC_NOTES.md">QC_NOTES.md</a>. Figures: <a href="''' + Q + '''sicista_psmc_ont_only_2026-10-04.png">genome-wide PSMC</a>, <a href="''' + Q + '''sicista_psmc_subset_tests_2026-10-04.png">subset tests</a>.</li>
    <li>Sb1 is not in a public archive yet; the mitogenome, nuclear mt pseudogenes and the SINE scan have their own pages.</li>
  </ul>
  <p style="margin-top:14px;"><a class="btn secondary" href="sicista.html">Sicista SINE page</a>
  <a class="btn secondary" href="sicista_mito.html">Mitogenomes and mt pseudogenes</a></p>
</section>
</main>
'''

page = ('<!doctype html>\n<html lang="en"><head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width,initial-scale=1">\n'
        '<title>Sicista betulina Sb1: assembly QC and the Illumina decision</title>\n' + css + '\n</head><body>\n' + body + '</body></html>\n')
open(os.path.join(ROOT, 'sicista_qc.html'), 'w', encoding='utf-8', newline='\n').write(page)
print('wrote sicista_qc.html', len(page))
