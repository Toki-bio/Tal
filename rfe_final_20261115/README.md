# rfe_final_20261115

*Rhinolophus ferrumequinum* (rfe, GCF_004115265.2), run identical to `rsi_final_20261102`: SINEderella main f2457a8, bank `bank21_nor1` (the *Rhinolophus sinicus* bank, 21 consensuses; r1_9seqs removed, r10_r8_P18del added), families as `rsi_families_nor1.tsv` (Rhin, MEG); regular search, composite-piece search, end extension, unit/composite vote rule, 5' flank mask, border loop with the junction design. therioserver `/home/toki/rhin/rfe_final_20261115/` (bank and launcher `/home/toki/rhin/final_bank_20261026/launch_final_20261102.sh`, code `/home/toki/SINEderella_main_20261102`); flank scan stages 1-7 and the report rebuild as for rsi_final_20261102.

- Run: start 2026-11-15 19:21:58 rfe; done 2026-11-15 21:03:48.

## Verdicts (every copy; `loci_verdict.tsv.gz`, rules SINEderella `docs/VERDICTS.md`)

| | copies |
|---|---|
| copies in the list | 63,093 |
| assigned to a subfamily | 59,684 |
| SINE loci (all copies) | 48,212 |
| grey zone (all copies) | 10,921 |
| not SINE (all copies) | 3,960 |

Copies without a subfamily are listed as not SINE in this table, partial copies included (below).

Per family (copies assigned to a subfamily):

| family | subfamilies with copies | assigned | SINE | grey | not SINE |
|---|---|---|---|---|---|
| MEG | 4 | 126 | 114 | 8 | 4 |
| Rhin | 17 | 59,558 | 48,098 | 10,913 | 547 |

## Family roll-up and partial copies (`family_view/`)

| family | copies | subfamily-firm | near tie | family only | not firm |
|---|---|---|---|---|---|
| MEG | 144 | 126 | 5 | 18 | 0 |
| Rhin | 62,949 | 59,558 | 4,347 | 3,382 | 9 |

Partial copies (family-firm, covering < 60 % of every competing consensus): 1,204, plus 278 below their stretch bar; 93 stretches (`partial_stretches.tsv`, `partial_copies.tsv`). Largest stretches:

| stretch | copies | median length | spans |
|---|---|---|---|
| MEG_partial1 | 1 | 111 | MEG-RL:1-120, MEG-TR:65-175 |
| Rhin_partial1 | 396 | 105 | r10_r8_P18:292-426, r10_r8_P18del:147-253 |
| Rhin_partial2 | 82 | 121 | r1_r3_P1:308-429, r1_r3_P2:165-286 |
| Rhin_partial3 | 61 | 77 | r10_r8_P18del:1-79, r10_19seqs:1-74 |
| Rhin_partial4 | 49 | 143 | r1_r3_P2:144-286, r1_r3_P1:288-429 |
| Rhin_partial5 | 46 | 119 | r1_r2_r3_r3_C11:1-130, r1_r3_P1:1-121, r1_r3_P2:1-112 |
| Rhin_partial6 | 43 | 99 | r1_r3_P1:332-429, r1_r3_P2:189-286 |
| Rhin_partial7 | 40 | 89 | r1_r3_P1:1-89, r1_r3_P2:1-89 |

## Contents
| Path | What |
|---|---|
| `report.html` | full report |
| `loci_verdict.tsv.gz`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |
| `family_view/` | family roll-up per copy, partial copies, stretches, the composite stretches searched in step 1 |
| `alignments/` | plates (SINE loci, topped up with grey copies marked [grey]) |
| `alignments/<sp>_all_copies_subfam30k.aln.fa` | step 1 SubFam: 600 chunk consensuses of 30,000 random copies of all subfamilies, aligned with the bank consensuses |
| `end_extension/` | step 1 end extension: per group and side (`ends_extension.tsv`), every moved locus (`loci_extended.tsv`) |
| `trna_plates/` | each consensus against its best tRNA gene (SINEderella `tools/trna_plates.py`; `index.html`, `index.tsv`) |
| `continuation/`, `flankscan/` | continuation check; flank scan stages 1-7 |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers |
| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | flank twins, tandem-array share, satellite screen |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its copies |
| `composite_units.tsv`, `consensus_blocks.tsv` | composites read as units; the bank as blocks |
| `consensuses.clean.fa`, `manifest.txt`, `pipeline_stages.txt` | the bank; tool versions and command; stage times |
