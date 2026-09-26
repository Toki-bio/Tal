# Rhin-1 SINE in Rhinolophidae — SINEderella runs (2026-09-26)

Separate from the Tal (Talpidae) analysis. Server: therioserver, `~/rhin/`.

## Consensus
`Rhin-1` [Borodulina 2005], "bats Rhinolophidae & Hipposideridae", 182 bp, from the SINEbase
download (`~/SINEs_post2015.fas` record 301 — that file is misnamed, it is SINEbase).
File: `alignments/rhin_consensuses.fa`.

**Note:** SINEderella's consensus sanitizing deletes IUPAC ambiguity codes instead of replacing
them. Rhin-1 has five (Y R Y R W), so all three searches used a 177 bp consensus with five
single-base deletions, and the anchor row in the alignments is that 177 bp version.

## Genomes (best Rhinolophidae assemblies by contig N50, one per species)
| code | species | accession | level | contig N50 |
|---|---|---|---|---|
| rsi | *Rhinolophus sinicus* | GCF_057876585.1 (RSI-T2T) | chromosome, T2T | 189 Mb |
| rre | *Rhinolophus rex* | GCA_041825415.1 | chromosome | 75 Mb |
| rda | *Rhinolophus darlingi* | GCA_057427945.1 (mRhiDar1.hap1) | chromosome | 58 Mb |

hap1 was used for *R. darlingi*; hap2 (GCA_057428265.1) was rejected as the secondary haplotype.

## Runs
SINEderella full mode, Rhin-1 as the only consensus, `CHUNK_BP=30000 FLANK=50 THREADS=24`.

| code | run | Rhin-1 hits (step1) | copies extracted | firm assigned | soft | sim_ratio median |
|---|---|---|---|---|---|---|
| rsi | run_20260926_142634 | 30,081 | 30,079 | 23,545 (78.3 %) | 6,534 | 0.480 |
| rre | run_20260926_143118 | 32,414 | 32,410 | 23,073 (71.2 %) | 9,337 | 0.482 |
| rda | run_20260926_143737 | 30,170 | 30,169 | 23,850 (79.1 %) | 6,319 | 0.484 |

`sear` in step1 runs without a hit cap (`-s` not passed), so the ~30,000 counts are real, not a cap.
rsi step2 was killed once by an unrelated watcher and completed with `SINEderella --resume`.

## Alignments (for manual subfamily grouping)
`alignments/<code>_subfam_input_30k.aln.fa` — SubFam of the step1 sample: 600 chunk consensi
(50 copies each) aligned together with the Rhin-1 anchor row (601 rows), from
`genome.clean_step1/subfam_input/input.clw.al`. No subfamilies defined yet.
