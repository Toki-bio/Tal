# rsi_sat — the rsi run with the satellite screen on (2026-10-03)

Same genome (GCF_057876585.1) and 21-consensus bank as `rsi_final/`, run with SINEderella `aa178c3`+ (satellite stage inside step 1,
`--exclude-b flagged` mode; docs/SATELLITES.md). Server: therioserver `~/rhin/rsi_sat/run_20261003_202315`.

## What the stage did
* **MEG-RS: tandem array** (SAT_B, 91 % excess over chance): 21 regular runs, **1 600 of 1 755 hits removed** before extraction and assignment.
  128 copies remain assigned (rsi_final: 1 715); the array flag of the remaining copies is 0 %. Those 128 are the material for the question
  "does MEG-RS exist as a dispersed SINE in rsi".
* **MEG-TR: SAT_B** on one run of 72 units (41.6 % of 173 hits); 72 hits removed, 19 assigned (21 before).
* **34 distinct SINE-derived satellite loci (kind A)** across the bank, each attributed to the consensus it aligns to best (`satellites/loci.bed`,
  column `also_matches` lists the others). The largest: NC_142509.1:30 184 174-30 192 326, a **140 bp monomer x 58 copies** that is the head + body
  of the r-families (matches 17 consensuses; attributed to C11 at 84-359): a SINE-derived satellite in a bat. Others: 17 loci of P18-like monomers
  (61-124 bp, up to 38 copies), 11 of C11, 3 of MEG-T2 (65-89 bp). Excluded hits of these: 2-74 per family (`satellites/excluded_hits.bed`).
* Everything else unchanged within assignment noise (+-31 copies; r6 +280 is the known run-to-run variation), length-version verdicts identical,
  r8 dimer loci (merged, >= 330 bp) 31 vs 29.

## Files
`report.html` (with the Satellites table in the Similarity section), `alignments/`, `satellites/` (indication.tsv, loci.bed, units.fa, excluded_hits.bed,
per-consensus kindA loci), `array_flag.tsv`, `length_variants/`, `assignment_stats.tsv`, `consensuses.clean.fa`.
Browse the loci: https://toki-bio.github.io/Tal/satellites/viewer.html?loci=rsi_sat/satellites/loci.bed&ind=rsi_sat/satellites/indication.tsv

## Known limits
* Kind-B runs are geometric (regular spacing); the units are not yet checked for sequence similarity. `--exclude-b long` is therefore not the default.
* Kind-A loci attributed by aligned length x identity favour long composite consensuses (C11, P18) over the unit they really derive from.
* The run used the tool copies of its start; the tables in `satellites/` were recomputed on the original hits with the current tools (attribution).
