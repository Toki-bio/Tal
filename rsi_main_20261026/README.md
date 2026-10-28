# rsi_main_20261026

*Rhinolophus sinicus* (rsi), assembly GCF_057876585.1, final run of the rhinolophid bank. Bank: the *Rhinolophus sinicus* bank without r1_9seqs (21 consensuses: r2-r10, r10_r8_P18del, seven composites, four MEG SINEs; therioserver `/home/toki/rhin/final_bank_20261026/bank21_nor1.fa`), families Rhin and MEG (`families.tsv`). SINEderella main a9d2a51 (2026-10-26): composite stretches of >= 50 bp searched as queries of their own in step 1; partial copies labelled and grouped into stretches in the family roll-up (stretch tables regrouped with main 9d484fa); flank scan stages 1-7 and the report rebuild as for rsi_main_20261022.

- Server: therioserver `/home/toki/rhin/rsi_main_20261026/` (`../final_bank_20261026/launch_final.sh`, `pack_final.sh`).
- Run: start 2026-10-26 20:30:20 rsi; done 2026-10-26 22:34:17.

## Verdicts (every copy; `loci_verdict.tsv.gz`, rules SINEderella `docs/VERDICTS.md`)

| | copies |
|---|---|
| copies in the list | 68,581 |
| assigned to a subfamily | 62,650 |
| SINE loci (all copies) | 51,399 |
| grey zone (all copies) | 8,773 |
| not SINE (all copies) | 8,409 |

Copies without a subfamily are listed as not SINE in this table, partial copies included (below).

Per family (copies assigned to a subfamily):

| family | subfamilies with copies | assigned | SINE | grey | not SINE |
|---|---|---|---|---|---|
| MEG | 4 | 868 | 124 | 731 | 13 |
| Rhin | 17 | 61,782 | 51,275 | 8,042 | 2,465 |

## Family roll-up and partial copies (`family_view/`)

| family | copies | subfamily-firm | near tie | family only | not firm |
|---|---|---|---|---|---|
| MEG | 891 | 868 | 10 | 23 | 0 |
| Rhin | 67,690 | 61,782 | 1,945 | 5,897 | 11 |

Partial copies (family-firm, covering < 60 % of every competing consensus): 1,421, plus 350 below their stretch bar; 95 stretches (`partial_stretches.tsv`, `partial_copies.tsv`). Largest stretches:

| stretch | copies | median length | spans |
|---|---|---|---|
| Rhin_partial1 | 284 | 136 | r10_r8_P18:290-426 |
| Rhin_partial2 | 175 | 187 | r1_r2_r3_r3_C11:1-189, r1_r3_P1:1-189 |
| Rhin_partial3 | 122 | 112 | r10_r8_P18:312-426, r10_r8_P18del:157-253 |
| Rhin_partial4 | 80 | 73 | r10_r8_P18del:1-78, r10_19seqs:2-73 |
| Rhin_partial5 | 69 | 271 | r1_r2_r3_r3_C11:335-607 |
| Rhin_partial6 | 66 | 121 | r1_r3_P2:165-286, r1_r3_P1:309-429 |
| Rhin_partial7 | 49 | 84 | r1_r3_P1:1-85, r1_r3_P2:1-85 |
| Rhin_partial8 | 38 | 102 | r1_r3_P1:327-429, r1_r3_P2:184-286 |

## Against rsi_main_20261022 (bank with r1_9seqs, before the piece search and the partial rule)
Copies matched by coordinate overlap (therioserver `/home/toki/rhin/final_bank_20261026/moved_rsi22_final.txt`):
- the 21 loci that r1_9seqs had: all 21 found; 20 are partial copies of the shared head (family Rhin), 1 is a full
  r1_r3_P1 copy. Without the piece search 3 of them were not extracted (`rsi_nor1_20261026` test);
- copies in the list: 68,581 against 67,192; the composite pieces add loci, most of them not assigned (too short);
- the largest move: 1,342 r6_210seqs copies to "no unanimous vote" (215 the other way). These copies split their ten votes
  between r6 and r5_r6_P26, which contains r6 at 100 % identity (positions 138-376), so the winner of each cycle is chance;
  the same r6 spread was seen between earlier runs (502 copies). At family level they are Rhin;
- every other subfamily changes by 90 copies or fewer (r9_15seqs +90, r5_r6_P26 +61, r10_r8_P18 +60, r10_19seqs -42);
  the small r3_58seqs went from 52 to 23 assigned copies.
Why r1_9seqs left the bank: `../rsi_main_20261022/README.md` ("Which r1 module is r2's?", "Is r1 more than r10 + r4?",
"Rerun without r1_9seqs") and SINEderella `docs/PARTIAL_COPIES.md`.

## Contents
| Path | What |
|---|---|
| `report.html` | full report |
| `loci_verdict.tsv.gz`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |
| `family_view/` | family roll-up per copy, partial copies, stretches, the composite stretches searched in step 1 |
| `alignments/` | plates (SINE loci, topped up with grey copies marked [grey]) |
| `alignments/rsi_all_copies_subfam30k.aln.fa` | step 1 SubFam: 600 chunk consensuses of 30,000 random copies of all subfamilies, aligned with the 21 bank consensuses |
| `continuation/`, `flankscan/` | continuation check; flank scan stages 1-7 |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers |
| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | flank twins, tandem-array share, satellite screen |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its copies |
| `composite_units.tsv`, `consensus_blocks.tsv` | composites read as units; the bank as blocks |
| `consensuses.clean.fa`, `manifest.txt`, `pipeline_stages.txt` | the bank; tool versions and command; stage times |
