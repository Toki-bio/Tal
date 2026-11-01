# rsi_r11_20261101

**Superseded (owner, 2026-11-01):** r11 is not a SINE. Its plates show no A box, a poor B box and a ragged 5' edge with good tails: standalone r11 copies are 5'-truncated copies of the second unit of r10_r8_P18. Without r11 in the bank the partial-copy rule reports them as such (in `rsi_main_20261026`, 411 of these 571 copies are partial copies; stretch Rhin_partial1, 284 copies, = r10_r8_P18 290-426).

*Rhinolophus sinicus* (rsi), rerun with r11 added to the bank. Bank: `bank22_r11` = the rsi bank21_nor1 (as `rsi_main_20261026`) + `r11_4003seqs`, the second unit of r10_r8_P18 / r10_r8_P18del built from 4,003 copies cut at their junctions (`../rsi_main_20261026/unit_pieces_copies/`, section "r11"); families as bank21_nor1, r11 in Rhin. SINEderella main f19fcf1 (main a9d2a51 of the 2026-10-26 run + step 1 end extension + tRNA plates); flank scan stages 1-7 and the report rebuild as for the 2026-10-26 runs. therioserver `/home/toki/rhin/rsi_r11_20261101/` (bank and launcher `/home/toki/rhin/bank22_r11_20261101/`).

- Run: start 2026-11-01 13:56:06 rsi; done 2026-11-01 17:31:50.

## Verdicts (every copy; `loci_verdict.tsv.gz`, rules SINEderella `docs/VERDICTS.md`)

| | copies |
|---|---|
| copies in the list | 68,793 |
| assigned to a subfamily | 64,133 |
| SINE loci (all copies) | 52,624 |
| grey zone (all copies) | 8,178 |
| not SINE (all copies) | 7,991 |

Copies without a subfamily are listed as not SINE in this table, partial copies included (below).

Per family (copies assigned to a subfamily):

| family | subfamilies with copies | assigned | SINE | grey | not SINE |
|---|---|---|---|---|---|
| MEG | 4 | 871 | 124 | 2 | 745 |
| Rhin | 18 | 63,262 | 52,500 | 8,176 | 2,586 |

## Family roll-up and partial copies (`family_view/`)

| family | copies | subfamily-firm | near tie | family only | not firm |
|---|---|---|---|---|---|
| MEG | 892 | 871 | 24 | 21 | 0 |
| Rhin | 67,901 | 63,262 | 1,962 | 4,635 | 4 |

Partial copies (family-firm, covering < 60 % of every competing consensus): 963, plus 361 below their stretch bar; 96 stretches (`partial_stretches.tsv`, `partial_copies.tsv`). Largest stretches:

| stretch | copies | median length | spans |
|---|---|---|---|
| Rhin_partial1 | 91 | 203 | r1_r3_P1:2-204, r1_r2_r3_r3_C11:2-198 |
| Rhin_partial2 | 79 | 73 | r10_r8_P18del:2-78, r10_19seqs:2-73 |
| Rhin_partial3 | 68 | 270 | r1_r2_r3_r3_C11:336-607, r1_r3_P1:200-447 |
| Rhin_partial4 | 67 | 121 | r1_r3_P2:165-286, r1_r3_P1:309-429 |
| Rhin_partial5 | 67 | 88 | r5_r5_P48:1-138 |
| Rhin_partial6 | 39 | 167 | r1_r2_r3_r3_C11:2-168, r1_r3_P1:2-168 |
| Rhin_partial7 | 38 | 86 | r1_r3_P2:2-87, r1_r3_P1:2-85 |
| Rhin_partial8 | 37 | 102 | r1_r3_P1:327-429, r1_r3_P2:184-286 |

## r11 in this run

`composite_units.tsv` now reads both P18 composites with r11:

| composite | r10 | middle | r11 |
|---|---|---|---|
| r10_r8_P18 | 6-150 (100 %) | 151-291 (no unit) | 292-426 (96.8 %; next r8 78.3 %) |
| r10_r8_P18del | 6-150 (98.6 %) | - | 151-255 (92.2 %; next r8 78.0 %) |

r10_r8_P18del was a unit of its own in `rsi_main_20261026`; it is now a composite.

Standalone r11 copies: 571 assigned (430 SINE). In `rsi_main_20261026` they were 456 unassigned (417 rejected for low
bitscore, 39 without a unanimous vote), 96 not extracted, 17 r10_r8_P18del, 2 other. Lengths 84-423 bp, median 134.
r8_83seqs is unchanged (3,388 -> 3,389 assigned).

## Against rsi_main_20261026 (`compare_to_rsi_main_20261026/`)

Copies matched by >= 50 % reciprocal overlap on the same strand (`crosstab.sh`). Assigned 62,650 -> 64,133; SINE
verdicts 51,399 -> 52,624. Moves of more than 50 copies:

| old | new | copies |
|---|---|---|
| unassigned | r6_210seqs | 967 |
| unassigned | r11_4003seqs | 456 |
| not extracted | unassigned | 156 |
| not extracted | r11_4003seqs | 96 |
| r5_r6_P26 | unassigned | 73 |
| unassigned | r1_r3_P1 | 67 |
| r6_210seqs | unassigned | 63 |
| unassigned | not extracted | 55 |

The r6_210seqs gain is not from r11 and not from the end extension (no r6 locus moved). The 450 of these copies with
the same coordinates in both runs had r6 as best but no unanimous vote in the old run (306 of them 9/10); the other votes
went to r5_r6_P26 (658 of 676). In this run they vote 10/10 r6. No step 2 code changed between the runs; the cause is
not determined (r5_r6_P26 contains r6 at 100 %, so the two score almost alike on these copies).

## Contents
| Path | What |
|---|---|
| `report.html` | full report |
| `loci_verdict.tsv.gz`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |
| `family_view/` | family roll-up per copy, partial copies, stretches, the composite stretches searched in step 1 |
| `alignments/` | plates (SINE loci, topped up with grey copies marked [grey]) |
| `alignments/<sp>_all_copies_subfam30k.aln.fa` | step 1 SubFam: 600 chunk consensuses of 30,000 random copies of all subfamilies, aligned with the bank consensuses |
| `end_extension/` | step 1 end extension: per group and side (`ends_extension.tsv`), every moved locus (`loci_extended.tsv`) |
| `compare_to_rsi_main_20261026/` | per-copy cross-tab against rsi_main_20261026 (`crosstab.tsv`, `by_new.tsv`, `by_old.tsv`, `crosstab.sh`) |
| `trna_plates/` | each consensus against its best tRNA gene (SINEderella `tools/trna_plates.py`; `index.html`, `index.tsv`) |
| `continuation/`, `flankscan/` | continuation check; flank scan stages 1-7 |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers |
| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | flank twins, tandem-array share, satellite screen |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its copies |
| `composite_units.tsv`, `consensus_blocks.tsv` | composites read as units; the bank as blocks |
| `consensuses.clean.fa`, `manifest.txt`, `pipeline_stages.txt` | the bank; tool versions and command; stage times |
