# SVK — Neotetracus sinensis — analysis log

**Status: preliminary, pending manual subfamily inspection. Not a finished result.**

## Provenance

- Specimen: SVK23-170, hybrid ONT + Illumina/MGI *de novo* assembly
  (Flye + Medaka + purge_dups), run on DRAGEN.
- Assembly used: `purge_dups_out/purged.fa` (2.48 GB, 33,740 contigs).
- SINE search: SINEderella step1 (`sear -k`, Slf=0.8, Homology=65, Flank=50)
  against the same combined consensus bank used for AVA (SINEbase nr95 +
  post-2015 literature supplement, 350 consensuses).
- Step1 search produced 1,890,839 candidate loci across the genome.
- **Full-genome subfamily assignment (step2/3/4/6) was deliberately NOT
  run yet** — same reasoning as AVA: subfamilies must be defined by manual
  inspection of a representative sample first.

## What this alignment is

`alignments/svk_subfam_input_30k.aln.fa` — a random 30,000-locus sample from
the 1,890,839 step1 hits, chunked into blocks of 50 by `SubFam`, yielding
600 chunk-consensus rows for manual subfamily grouping. Same method as AVA;
see `../ava/LOG.md` for the full description of what this file is (and
is not).
