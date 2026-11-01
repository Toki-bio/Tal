# rsi_final_20261102

*Rhinolophus sinicus* (rsi), final run with the current SINEderella main (f2457a8). Bank: `bank21_nor1` (21 consensuses, as `rsi_main_20261026`; r1_9seqs removed, r10_r8_P18del added), families as `rsi_families_nor1.tsv`; regular search, composite-piece search, end extension, unit/composite vote rule, 5' flank mask, border loop with the junction design. therioserver `/home/toki/rhin/rsi_final_20261102/` (bank and launcher `/home/toki/rhin/final_bank_20261026/`, code `/home/toki/SINEderella_main_20261102`); flank scan stages 1-7 and the report rebuild as for the 2026-10-26 runs.

- Run: start 2026-11-01 23:56:02 rsi; done 2026-11-02 01:48:13.

## Verdicts (every copy; `loci_verdict.tsv.gz`, rules SINEderella `docs/VERDICTS.md`)

| | copies |
|---|---|
| copies in the list | 68,580 |
| assigned to a subfamily | 64,697 |
| SINE loci (all copies) | 50,790 |
| grey zone (all copies) | 12,404 |
| not SINE (all copies) | 5,386 |

Copies without a subfamily are listed as not SINE in this table, partial copies included (below).

Per family (copies assigned to a subfamily):

| family | subfamilies with copies | assigned | SINE | grey | not SINE |
|---|---|---|---|---|---|
| MEG | 4 | 858 | 117 | 5 | 736 |
| Rhin | 17 | 63,839 | 50,673 | 12,399 | 767 |

## Family roll-up and partial copies (`family_view/`)

| family | copies | subfamily-firm | near tie | family only | not firm |
|---|---|---|---|---|---|
| MEG | 887 | 858 | 4 | 29 | 0 |
| Rhin | 67,693 | 63,839 | 4,500 | 3,836 | 18 |

Partial copies (family-firm, covering < 60 % of every competing consensus): 1,355, plus 334 below their stretch bar; 88 stretches (`partial_stretches.tsv`, `partial_copies.tsv`). Largest stretches:

| stretch | copies | median length | spans |
|---|---|---|---|
| Rhin_partial1 | 407 | 105 | r10_r8_P18:293-426, r10_r8_P18del:147-253 |
| Rhin_partial2 | 99 | 192 | r1_r3_P1:1-192, r1_r2_r3_r3_C11:1-189 |
| Rhin_partial3 | 78 | 73 | r10_r8_P18del:1-77, r10_19seqs:2-72 |
| Rhin_partial4 | 67 | 120 | r1_r3_P1:308-429, r1_r3_P2:165-286 |
| Rhin_partial5 | 66 | 271 | r1_r2_r3_r3_C11:336-607 |
| Rhin_partial6 | 49 | 84 | r1_r3_P1:1-85, r1_r3_P2:1-85 |
| Rhin_partial7 | 37 | 102 | r1_r3_P1:327-429, r1_r3_P2:184-286 |
| Rhin_partial8 | 35 | 139 | r5_r6_P26:1-141, r5_r5_P48:1-140 |

## What is new in this run (SINEderella main f2457a8, against `rsi_main_20261026`)

Same genome and the same bank (bank21_nor1, 21 consensuses). Four changes in the code, each tested on its own before:

- **Unit versus composite in every voting cycle** (`unit_rule_swaps.tsv`, `unit_rule_units.tsv`): a composite beats the unit
  it contains only if the copy covers at least half of the composite outside the unit and at least half of the unit's own
  stretches. 3,416 copies had at least one cycle changed. Largest effect, r6_210seqs against r5_r6_P26: 7,442 copies, of which
  5,285 as before, 2,045 previously unassigned (the split-vote copies of the review `rsi_composites_test/r6_P26_20261101/`)
  and 112 from P26. Copies without a unanimous vote: 3,206 -> 1,123. r5_27seqs 1,981 (was 1,927; the first version of the
  rule had moved 1,522 of them into r5_r5_P48, see `docs/PARKED.md`).
- **5' flank repeat mask** (flank test): 5' flanks were never masked before; 3,000 to 4,000 copies whose 5' flank lies in a
  repeat are now grey (flanks not testable) instead of SINE.
- **Border loop**: the elements are aligned once per side and the flanks packed at the junction; the MEG-RS side verdicts
  are now computed (still shared on both sides) instead of the loop crashing.
- **End extension in step 1** (`end_extension/`), as in the hla run of 2026-10-30; and a plate of every consensus against its
  best tRNA gene (`trna_plates/index.html`; the consensus's second isotype is added when it scores within 2 bits).

Per-copy cross-tab against `rsi_main_20261026`: `compare_to_rsi_main_20261026/` (`crosstab.sh`, `crosstab.tsv`, `by_new.tsv`,
`by_old.tsv`). A judgement zone remains between r10_r8_P18 and r10_r8_P18del: 183 copies moved from P18del to P18 (median
length 313 bp; P18del copies 271, P18 copies 406), i.e. copies that carry part of P18's middle piece.

r11 (the second unit of P18) is not a subfamily here: its standalone copies are 5'-truncated copies and are reported as
partial copies (owner, 2026-11-01).

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
