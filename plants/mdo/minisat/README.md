# Apple (*Malus domestica*) minisatellite repertoire, satpipe run of 2026-11-26

Genome GCF_042453785.1 (GDT2T hap1, 651,332,481 bp, 17 chromosomes + 2 small sequences). Tool: satpipe (Toki-bio/CLsat_workflow, branch
`check-stage-aware`), published defaults: TRF `2 5 7 80 10 50 2000`, monomer 80-300 bp, at least 10 copies. Run on DRAGEN (Biotech2024):
`/staging/tmp/apple_minisat/run_fam`, config `apple_minisat_fam.yaml` here. Stage 19 (orthology) is meaningless for one genome and is ignored.

## How the families were found (peel)

1. Stages 02-04: TRF on the whole genome, motif extraction, simple-repeat filter: 301 motifs, 259 kept (`mdo_trf_filtered.aln.fa`).
2. Groups were found with a comparison that ignores phase and strand (each motif doubled, blastn both strands against all motifs; pair score =
   identity x covered fraction; groups = connected components at >= 0.70). satpipe's own stage 06-08 clustering is plus-strand only and split
   the largest family into two clusters (92 + 83 motifs), which this comparison joins.
3. For each group: medoid as reference, members oriented and rotated to the medoid phase, MAFFT `--localpair --maxiterate 1000`, consensus =
   majority base in columns with >= 50 % non-gap (`groups/`, `mdo_minisat_consensuses.fa`, `peel_group_stats.tsv`).
4. Control: the same 259 motifs with bases shuffled inside each motif give 0 groups of >= 3 (largest 1). 23 motifs are in no group; none matches a
   group at score >= 0.4.
5. The seven consensuses became the families of a second satpipe run (`anchor.independent: true`, the families are unrelated: no blastn or
   ssearch36 hit between any two). Cluster-to-family table `cluster_family_fam.tsv` from the peel membership of each stage-08 cluster.
6. Scaffolds over 50 Mb (here NC_091675.1, 56.2 Mb) are scanned in chunks by stage 02b and their arrays classified against the family consensuses:
   32 of 33 candidates match a family at >= 70 % identity (14 G01_95, 18 G05_178).

## Result (final catalogue `mdo_minisat_final_loci.tsv`; 239 loci, 592,160 bp = 0.091 % of the genome)

| family | monomer bp | loci | bp | chromosomes | within 1 Mb of a chromosome end (uniform expectation) | median / max locus bp |
|---|---|---|---|---|---|---|
| G01_95 | 95 | 181 | 298,266 | 17 | 10 (6 %; 5.3 %) | 1,429 / 6,341 |
| G05_178 | 178 | 26 | 82,609 | 5 | 25 (96 %; 4.0 %) | 3,050 / 7,406 |
| G02_170 | 170 | 14 | 65,112 | 8 | 0 (0 %; 5.5 %) | 3,346 / 18,924 |
| G03_169 | 169 | 4 | 110,886 | 1 | 0 (0 %; 4.8 %) | 6,819 / 93,558 |
| G04_210 | 210 | 8 | 23,852 | 2 | 0 (0 %; 6.0 %) | 2,454 / 4,165 |
| G07_95 | 95 | 5 | 5,705 | 1 | 5 (100 %; 4.3 %) | 1,141 / 1,141 |
| G06_104 | 104 | 1 | 5,730 | 1 | 0 | 5,730 |

Every locus is assigned to a single family with a large margin (status UNIQUE_SEED_WITH_COMPETITOR in all 239; competitor bits are noise level).

Pattern, from the table:
- **G01_95** is dispersed: 181 short arrays (mean locus length 17 monomers, longest 67; members of the peeled group are on average 99 % identical to the consensus) on all 17 chromosomes, with no excess near the ends.
- **G05_178** is a chromosome-end family: 25 of 26 loci lie within 1 Mb of an end, on five chromosomes (NC_091675.1 end: 18 loci; NC_091665.1 start: 5; NC_091664.1 end, NC_091672.1 start, NC_091677.1 start: 1 each).
- **G07_95** (a second 95 bp family, unrelated to G01_95 by sequence) is five identical copies within 45 kb of the start of NC_091665.1.
- **G03_169** is one cluster on NC_091671.1 at 27.50-27.62 Mb, including one 93.6 kb array.
- **G02_170** lies in the interior of eight chromosomes (11-70 % of chromosome length), single clusters of 1-5 loci, the largest 18.9 kb; **G04_210** on two chromosomes (NC_091661.1: 5 loci, NC_091677.1: 3), **G06_104** one 5.7 kb locus on NC_091662.1.
- All together the minisatellite fraction is small (0.09 %). I could not find a published apple satellite family to compare with (literature search 2026-11-26; apple centromeres are described as retrotransposon-based).

## Position relative to the published centromere ranges

Ranges: Su et al. 2024, Plant Physiol 195:2799 (PMID 38743633), Supplementary Table S6 (`GDT2T_hap1_chrNN`, start/end; defined roughly from transposon abundance, Gypsy and HODOR, so Mb-scale ranges), copied in `gdt2t_hap1_centromeres.tsv`. The paper gives no accession mapping; chromosome NN = NCBI "chromosome NN" in the FASTA header of GCF_042453785.1 (NC_091661.1 = 01 ... NC_091677.1 = 17). Table S6 has no chromosome lengths, so the mapping could not be cross-checked; the paper gives the longest hap1 scaffold as 56,139,794 bp against 56,249,447 bp for NC_091675.1. Per locus in `loci_vs_centromere.tsv` (script `scripts/cen.sh`).

| family | loci | inside a centromere range | within 1 Mb of it | expected if uniform (inside; within 1 Mb) |
|---|---|---|---|---|
| G01_95 | 181 | 0 | 6 (3 %) | 2.4 %; 7.7 % |
| G02_170 | 14 | 2 (14 %) | 6 (43 %) | 2.2 %; 7.8 % |
| G03_169 | 4 | 0 | 0 | 1.8 %; 6.6 % |
| G04_210 | 8 | 0 | 0 | 2.7 %; 8.7 % |
| G05_178 | 26 | 0 | 0 | 1.6 %; 5.6 % |
| G06_104 | 1 | 0 | 0 | 7.4 %; 12.7 % |
| G07_95 | 5 | 0 | 0 | 1.8 %; 6.0 % |

G02_170 loci lie on 8 chromosomes; on 5 of them (chr01, chr08, chr10, chr11, chr14) a locus is within 0.4 Mb of the centromere range (inside it on chr01 and chr11; 60 kb on chr14, 164 kb on chr08, 392 kb on chr10); on chr09, chr12 and chr13 the nearest locus is 1.9, 3.3 and 12.3 Mb away. G02_170 is therefore centromere-proximal on most chromosomes that carry it, but it is a small family (65 kb) and the ranges are rough. No other family is near a centromere; G05_178 and G07_95 are at chromosome ends.

## What the families are (NCBI BLAST, 2026-11-26; raw hit lists in `blast_hits.txt`)

blastn (word size 11) of each consensus against NCBI core nt; every family also matches the apple genome itself (expected). Beyond that:

| family | other hits | reading |
|---|---|---|
| G04_210 | mRNAs of *Pyrus* and *Malus sylvestris* annotated "disease resistance protein", "TMV resistance protein N-like", "protein SUPPRESSOR OF npr1-1, CONSTITUTIVE 1", 87-93 % identity over the full 210 bp | a repeated stretch of a TIR-NBS-LRR-type resistance gene (translation of the consensus has no stop codon and the pattern LNLSGC...ELPSSI) |
| G07_95 | ribosomal RNA gene entries of *Pyracantha*, *Crataegus*, *Aria*, *Sorbus* (SSU/LSU/ITS), 92-96 % over 76-95 bp | derived from rDNA |
| G06_104 | *M. sylvestris* mRNAs of a peptidyl-prolyl cis-trans isomerase (97 %) and an uncharacterised lncRNA of *M. domestica* (96 %), over the whole 104 bp | contains gene or transcribed sequence |
| G05_178 | *M. sylvestris* / *M. domestica* acetyl-CoA synthetase mRNAs, 88-92 % over 133-179 bp (Rosaceae-restricted search) | contains gene-derived sequence |
| G03_169 | *Prunus persica* chromosomes 1-3 at 76-80 % over 91-109 bp, *Argentina anserina* and *P. persica* mRNA fragments at 76-82 % (Rosaceae-restricted search) | weak similarity outside *Malus*; the 166 bp *Prunus* centromeric monomer reported by Wang et al. 2022 (Genes 13:641) is a candidate relation, not tested (monomer sequence not at hand) |
| G01_95, G02_170 | top 12 hits are apple chromosomes only (G01 100 %, G02 85-92 %) | no information from the 12-hit lists |

A match to a gene says that the repeat unit is gene-derived (or that the gene contains the repeat); it does not show that the family is a gene.
The unrestricted G03 and G05 searches and the protein search for G04 stayed queued at NCBI for more than 30 minutes and were not used.

## Alignments of each peeled group (ViewAlign)

- G01 (177 motifs): https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fplants%2Fmdo%2Fminisat%2Fgroups%2Fmdo_minisat_G01.aln.fa
- G02 (16 motifs): https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fplants%2Fmdo%2Fminisat%2Fgroups%2Fmdo_minisat_G02.aln.fa
- G03 (11 motifs): https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fplants%2Fmdo%2Fminisat%2Fgroups%2Fmdo_minisat_G03.aln.fa
- G04 (11 motifs): https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fplants%2Fmdo%2Fminisat%2Fgroups%2Fmdo_minisat_G04.aln.fa
- G05 (10 motifs): https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fplants%2Fmdo%2Fminisat%2Fgroups%2Fmdo_minisat_G05.aln.fa
- G06 (6 motifs): https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fplants%2Fmdo%2Fminisat%2Fgroups%2Fmdo_minisat_G06.aln.fa
- G07 (5 motifs): https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fplants%2Fmdo%2Fminisat%2Fgroups%2Fmdo_minisat_G07.aln.fa

## Not established

- G03_169 and G05_178: the Rosaceae-restricted BLAST gives weak (G03) or gene-fragment (G05) similarity only; whether either is a known satellite is open.
- One haplotype only.
- satpipe's cluster stage is plus-strand only; two clusters of one family (G01_95) were merged by hand-built peel, not by the tool.
