# AVA — Neohylomys vietnamensis — analysis log

**Status: preliminary, pending manual subfamily inspection. Not a finished result.**

## Provenance

- Specimen: AVA18-056, hybrid ONT + Illumina/MGI *de novo* assembly
  (Flye + Medaka + purge_dups), run on DRAGEN.
- Assembly used: `purge_dups_out/purged.fa` (2.46 GB, 27,753 contigs).
- SINE search: SINEderella step1 (`sear -k`, Slf=0.8, Homology=65, Flank=50)
  against a combined consensus bank of SINEbase (nr95) + a post-2015
  literature supplement (350 consensuses total, 110 duplicates removed on
  merge).
- Step1 search produced 1,800,089 candidate loci across the genome.
- **Full-genome subfamily assignment (step2/3/4/6) was deliberately NOT
  run yet.** Per standard SINEderella workflow, subfamilies are defined by
  manual inspection of a representative sample before any full assignment —
  running step2 against an uncurated consensus bank first would be premature
  and produce meaningless subfamily calls.

## What this alignment is

`alignments/ava_subfam_input_30k.aln.fa` — a random 30,000-locus sample from
the 1,800,089 step1 hits, chunked into blocks of 50 by `SubFam` (MAFFT
guide-tree ordering + per-chunk consensus), yielding 600 chunk-consensus
rows. This is the standard SINEderella manual-review artifact: each row is a
consensus of ~50 similar copies, ordered by sequence similarity, for
grouping into candidate subfamilies by eye (shared diagnostic
indels/SNPs) before building refined consensuses and re-running assignment.

**This is not itself a subfamily alignment** — it is the input to that
grouping step.

## Run history note

The first attempt at this run crashed partway through the final SubFam
alignment build (missing `rename` utility on the compute host); the run was
resumed and the final alignment rebuilt after fixing the environment issue.
1,800,089 loci and the 30,000-sequence sample are unaffected by that crash —
only the final chunk-consensus alignment build had to be redone.
