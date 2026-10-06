# Sicista mitogenomes and mitochondrial pseudogenes

## 2026-10-04: alignments for ViewAlign

Built on therioserver, `~/Sicista2026/mito/aln/` (inputs: `genbank/`, `units/`, `bet_numt/`, `trizona_numt/`). Page: `sicista_mito.html`
(generator `sicista/mito/build_mito_page.py`, style copied from `sicista.html`).

- `sicista_mitogenomes.aln.fa`: 61 GenBank Sicista mitochondrial records (5 species) + Sb1 own mtDNA (unit of ptg000633l). MAFFT --auto --adjustdirection.
- `sicista_control_region.aln.fa`: the region after tRNA-Pro (NC_069019 position 15,465) to the record end, extracted from the alignment above and realigned (MAFFT --localpair --maxiterate 100).
- `sb1_mtlike_contigs.aln.fa`: the 83 mt-like contigs of the hifiasm primary assembly (>=60% of length matches Sb1 mtDNA at >=95%), cut into 160 segments by minimap2 -x asm5 (>=800 bp) and added to the Sb1 mtDNA with MAFFT --addfragments --keeplength. 151 segments >=99.9% identical to the own mtDNA.
- `strizona_numts.aln.fa`: S. trizona nuclear loci (blastn of OZ418355.1 vs the other 161 sequences of GCA_982266845.1, e<1e-5, loci merged within 2 kb) with aligned length >=500 bp or >=95% identity and >=100 bp (57 of 141) placed on OZ418355.1.
- `sb1_nuclear_lineageB_family.aln.fa`: 30 evenly spaced copies of the >=82 full-length lineage-B nuclear copies (aligned >=14 kb, 90-95% to Sb1 mtDNA), the two young own-lineage copies (ptg000319l:246942-263266, ptg000880l:1-12009) and GenBank NC_069019, MZ570955, OZ418355, on Sb1 mtDNA.
- `sb1_nuclear_other_numts.aln.fa`: the 84 largest remaining nuclear loci (>=1 kb, or >=95% and >=100 bp) on Sb1 mtDNA.

Fragments are placed with `--keeplength`, so insertions relative to the backbone are dropped. Results and methods:
C:\work\hylomys_ont\SICISTA_MITO_FINDINGS.md.

## 2026-10-04: gene track (BED) for the ViewAlign pages

ViewAlign v233 draws a BED annotation track above the alignment (`?bed=` URL parameter; Settings > Annotation). Files in `sicista/mito/annotations/`:

- `sicista_mitogenomes.NC_069019_genes.bed`: the GenBank feature table of NC_069019.1 (RefSeq *S. betulina* DK_BMNT: 13 CDS, 22 tRNA, 2 rRNA, O_L, control region; `gene` features dropped as duplicates) on the row `Sbet_NC_069019_DK_BMNT`, 0-based half-open, itemRgb by feature class, column 13 = feature type and /product. Made with `gb2bed.py genbank/NC_069019.1.gb Sbet_NC_069019_DK_BMNT` on therioserver.
- `strizona_numts.OZ418355_genes_lifted.bed`, `sb1_numts.ptg633_unit2_genes_lifted.bed`: the same features lifted onto OZ418355.1 and onto Sb1 `ptg633_unit2` through `units/six.aln` (the six-mitogenome MAFFT alignment used earlier for `genes_tri.tsv` / `genes_unit2.tsv`) with `lift_bed.py`; 0 intervals dropped. OZ418355.1 has no feature table in its GenBank record, so the trizona track is a transfer, not an official annotation. These serve the NUMT alignments, whose first row is the mtDNA backbone.

The viewer maps the BED coordinates through the backbone row's gaps on every redraw, so the track is correct in all three alignments that share a backbone.

## 2026-10-06: page revised on the read-polished Sb1 mtDNA; codon-level pseudogene check

- The Sb1 mtDNA used so far (unit of `ptg000633l`, `ptg633_unit2`) has ONT indel errors: COX1 +CT, COX3 +TTTTC, ND5 -ACTC (read-proven, 80-85% of
  ~5,300-6,100 lineage-A reads), a non-coding +AAA, and a probable ND2 C-run C5->C6 (reads split 53/32/15%). Polished sequence:
  `sb1_own_mtDNA_polished.fa` (samtools consensus -X r10.4_sup + the ND2 C6), 16,661 bp, all 13 CDS stop-free.
- `aln_v2.sh` (therioserver `~/Sicista2026/mito/aln_v2/`) rebuilt the five Sb1 alignments with the polished sequence as the Sb1 row / backbone
  (same fragment sets as before); `strizona_numts.aln.fa` unchanged. Gene track for the Sb1 backbone: `annotations/sb1_own_mtDNA_polished_genes.bed`
  (NC_069019 features lifted through the new mitogenome alignment; the 13 CDS match the read-checked coordinates); the old
  `sb1_numts.ptg633_unit2_genes_lifted.bed` was removed.
- `codon/`: codon-level check of the 61 GenBank records against pseudogene controls (82-copy nuclear lineage-B family, young/old Sb1 copies,
  29 S. trizona loci), ViewAlign MACSE v2.07 port, vertebrate mito code, pairwise against the polished Sb1 mtDNA (therioserver
  `~/Sicista2026/mito/pseudo_check/`; scripts in `codon/scripts/`). GenBank S. betulina: 0 frameshifts, 0 premature stops, pN/pS 0.088 (A) / 0.078 (B);
  nuclear family: every copy disabled, pN/pS 1.01, shared stop COX1 codon 179 (82/82). N calls in lineage-B records enriched 4.4x at
  positions where the nuclear family differs (P 8e-9). Details: C:\work\hylomys_ont\SICISTA_MITO_FINDINGS.md sections 5-6.
