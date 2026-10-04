# rsi_sat — the rsi run with the satellite screen on (2026-10-03)

Same genome (GCF_057876585.1) and 21-consensus bank as `rsi_final/`, run with SINEderella `aa178c3`+ (satellite stage inside step 1,
`--exclude-b flagged` mode; docs/SATELLITES.md). Server: therioserver `~/rhin/rsi_sat/run_20261003_202315`.

## What the stage did
* **MEG-RS: tandem array** (SAT_B, 91 % excess over chance): 21 regular runs, **1 600 of 1 755 hits removed** before extraction and assignment.
  128 copies remain assigned (rsi_final: 1 715); the array flag of the remaining copies is 0 %. Those 128 are the material for the question
  "does MEG-RS exist as a dispersed SINE in rsi".
* **MEG-TR: SAT_B** on one run of 72 units (41.6 % of 173 hits); 72 hits removed, 19 assigned (21 before).
* **34 distinct SINE-derived satellite loci (kind A)** across the bank, each attributed to the consensus it derives from (identity within 3 points of the best, then the
  consensus the monomer covers most; `satellites/loci.bed`, column `also_matches` lists the others). The largest: NC_142509.1:30 184 174-30 192 326, a **140 bp monomer
  x 58 copies = the whole r2** (136 of 146 bp at 92.7 %; matches 17 consensuses through the shared head): an r2-derived satellite, with three more small r2 loci.
  Others: 17 P18-derived (61-124 bp monomers of the r10-r8 junction region, up to 38 copies), 3 r3, 3 C11, 3 MEG-T2 (60-89 bp), 1 each r8, r7, P1, r10.
  Excluded hits of these: 2-74 per family (`satellites/excluded_hits.bed`).
* **Kind-B unit check** (`satellites/rsi_sat.kindB.tsv`): of 1 904 regularly spaced runs, 42 have near-identical units (MEG-RS 21, MEG-TR 1, r9 15, r4 2, r10, r7, P26 1 each;
  among them a 26-unit 2.6 kb array at NC_142508.1:42.84 Mb carrying r9/r4/r7/r10 heads) and 1 862 are clusters of ordinary copies (unit identity 25-55 %).
  This run excluded in `flagged` mode (2 066 hits); the current default `verified` would exclude 2 243.
* Everything else unchanged within assignment noise (+-31 copies; r6 +280 is the known run-to-run variation), length-version verdicts identical,
  r8 dimer loci (merged, >= 330 bp) 31 vs 29.

## Files
`report.html` (with the Satellites table in the Similarity section), `alignments/`, `satellites/` (indication.tsv, loci.bed, units.fa, excluded_hits.bed,
per-consensus kindA loci), `array_flag.tsv`, `length_variants/`, `assignment_stats.tsv`, `consensuses.clean.fa`.
Browse the loci: https://toki-bio.github.io/Tal/satellites/viewer.html?loci=rsi_sat/satellites/loci.bed&ind=rsi_sat/satellites/indication.tsv

## Decide
All alignments with MSA-viewer links: [decide.html](https://toki-bio.github.io/Tal/rsi_sat/decide.html).

## Known limits
* The run used the tool copies of its start; the tables in `satellites/` were recomputed on the original hits with the current tools (attribution).
