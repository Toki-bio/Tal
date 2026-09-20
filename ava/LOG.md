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

## De novo arm (complete, also preliminary)

- **AnnoSINE2** (mode 3, hybrid homology+structural, `-a 2` for both plant
  and animal Dfam HMM families): 206/206 HMM families scanned, full 8-step
  pipeline (HMM search, structural search, TSD search, extension, IRF,
  clustering, RepeatMasker annotation) completed in ~9.5 h. Output:
  `annosine_out/Seed_SINE.fa`, 3,510 seed candidates.
- **SINE-de-novo-genome-scan** (`sine_scan.sh`): the same genome scanned
  against SINEbase's shredded fragment bank (`SINEBase.fragments.nr85.fa`,
  926 fragments, 5'/mid/3' region-tagged, no-library homology search).
  Output: 60,845 candidate loci.
- **Merged**: 3,510 + 60,845 = 64,355 combined candidates from both methods.
- **K-mer thinning** (`kmer_thin_singletons.py`, k=13, min 5 shared k-mers,
  min 1 clusterable partner — tentative first-pass parameters, not yet
  tuned against real review): dropped 4,964 candidates with no k-mer-level
  company (likely spurious singleton hits from either method), kept 59,391.
- **SubFam clustering** (chunks of 50, MAFFT guide-tree ordering + per-chunk
  consensus): 59,391 candidates &rarr; **1,187 chunk consensi**.

`alignments/ava_denovo_candidates_1187chunks.aln.fa` — the 1,187-row
clustered alignment. Same caveat as the database-arm sample: **this is
clustered, not assigned or classified as real SINEs** — de novo candidates
by construction include noise (fragment-of-longer-element hits, low-
complexity false positives, non-SINE repeats the fragment scan happened to
match). This is the input to your manual review, not a result.

One crash during this run: the first `sine_scan.sh` attempt (for SVK,
launched around the same time) filled `/tmp` (a 252 GB RAM-backed tmpfs)
with an orphaned 199 GB `parallel` buffer directory and crashed with a
disk-full error. Cleaned up before this AVA run completed; did not affect
AVA's own output.
