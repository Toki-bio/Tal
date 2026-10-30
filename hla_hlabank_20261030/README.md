# hla_hlabank_20261030

*Hipposideros larvatus* (hla), rerun with the hla-adjusted bank. Bank: `bank30_hla` (30 consensuses), the rsi bank21_nor1 adjusted from the peel of the hla step 1 SubFam plate: r8_83seqs, r6_210seqs, r5_r6_P26 and r10_r8_P18 replaced by r8_hla, r5_r6_short_hla, r5_r6_long_hla and consensus_386seqs; added r9_ext, consensus_28seqs, consensus_47seqs, consensus_10seqs_noterm, consensus_12seqs, consensus_10seqs_b, consensus_29seqs (owner's peel rounds 1-2) and r1_r3_P1_hla020, r5_r5_P48_hla007 (plate rows named in the colleague's notes). Build and provenance: `../hla_rsibank_20261026/peel/README.md` and `../hla_rsibank_20261026/peel/bank30_hla/`; therioserver `/home/toki/rhin/hla_bank_20261030/bank30_hla.fa`, families `families30_hla.tsv` (all new consensuses in Rhin). SINEderella branch `end-extension` bf2480a (main a9d2a51 + the step 1 end extension, `docs/END_EXTENSION.md`; not merged into main); flank scan stages 1-7 and the report rebuild as for the 2026-10-26 runs.

- Server: therioserver `/home/toki/rhin/hla_hlabank_20261030/` (`../hla_bank_20261030/launch_hla30.sh`, `build_bank.sh`, `pack_final.sh` with `FAM=families30_hla.tsv`).

- Run: start 2026-10-29 22:04:42 hla; done 2026-10-30 03:23:49.

## Verdicts (every copy; `loci_verdict.tsv.gz`, rules SINEderella `docs/VERDICTS.md`)

| | copies |
|---|---|
| copies in the list | 261,557 |
| assigned to a subfamily | 247,503 |
| SINE loci (all copies) | 188,923 |
| grey zone (all copies) | 28,608 |
| not SINE (all copies) | 44,026 |

Copies without a subfamily are listed as not SINE in this table, partial copies included (below).

Per family (copies assigned to a subfamily):

| family | subfamilies with copies | assigned | SINE | grey | not SINE |
|---|---|---|---|---|---|
| MEG | 4 | 583 | 115 | 9 | 459 |
| Rhin | 26 | 246,920 | 188,808 | 28,599 | 29,513 |

## Family roll-up and partial copies (`family_view/`)

| family | copies | subfamily-firm | near tie | family only | not firm |
|---|---|---|---|---|---|
| MEG | 630 | 583 | 55 | 47 | 0 |
| Rhin | 260,927 | 246,920 | 4,633 | 13,976 | 31 |

Partial copies (family-firm, covering < 60 % of every competing consensus): 8,159, plus 139 below their stretch bar; 225 stretches (`partial_stretches.tsv`, `partial_copies.tsv`). Largest stretches:

| stretch | copies | median length | spans |
|---|---|---|---|
| Rhin_partial1 | 1,072 | 99 | consensus_386seqs:1-99 |
| Rhin_partial2 | 750 | 70 | r5_r5_P48_hla007:251-323 |
| Rhin_partial3 | 584 | 170 | consensus_386seqs:1-170 |
| Rhin_partial4 | 553 | 111 | consensus_386seqs:328-438 |
| Rhin_partial5 | 525 | 191 | consensus_386seqs:248-438 |
| Rhin_partial6 | 468 | 154 | consensus_386seqs:285-438 |
| Rhin_partial7 | 266 | 171 | consensus_386seqs:267-438 |
| Rhin_partial8 | 241 | 147 | consensus_386seqs:1-147 |

## Against hla_rsibank_20261026 (same genome, rsi bank21_nor1)

| | hla_rsibank_20261026 | hla_hlabank_20261030 |
|---|---|---|
| copies in the list | 233,325 | 261,557 |
| assigned to a subfamily | 191,014 (81.9 %) | 247,503 (94.6 %) |
| SINE loci | 144,006 | 188,923 |
| grey zone | 21,012 | 28,608 |
| not SINE | 68,307 | 44,026 |
| largest group: mean similarity to its consensus | r10_r8_P18 0.365 | consensus_386seqs 0.876 |

Copies matched between the runs by >= 50 % reciprocal overlap on the same strand (`compare_to_hla_rsibank/`); 229,812 matched, 31,765 only in the new run, 5,743 only in the old run. Where the copies went (largest flows):

| old consensus | copies | new consensus (copies) |
|---|---|---|
| r10_r8_P18 | 140,029 | consensus_386seqs 138,267 |
| unassigned | 42,311 | consensus_47seqs 18,764; consensus_28seqs 7,983; still unassigned 5,608 |
| (not in the old list) | 31,765 | consensus_47seqs 22,931; unassigned 4,429; consensus_28seqs 3,473 |
| r5_r6_P26 | 14,868 | consensus_29seqs 9,542; r5_r6_long_hla 3,107; r1_r3_P1_hla020 515 |
| r8_83seqs | 11,625 | r8_hla 10,986 |
| r1_r2_r3_r3_C11 | 9,727 | consensus_10seqs_noterm 3,508; consensus_10seqs_b 2,665; consensus_12seqs 2,620 |
| r9_15seqs | 5,130 | r9_ext 4,178; r9_15seqs 796 |
| r6_210seqs | 2,938 | r5_r6_short_hla 2,060; r5_r6_long_hla 511 |
| r1_r3_P1 | 1,834 | consensus_12seqs 936; r1_r3_P1_hla020 652 |
| r10_r8_P18del | 1,423 | unassigned 895; consensus_386seqs 390; r10_r8_P18del 115 |

rsi consensuses kept in the bank that now collect few firm copies: r1_r2_r3_r3_C11 14, r1_r3_P1 8, r7_133seqs 67, r10_19seqs 36, r5_27seqs 185, r10_r8_P18del 144.

## Step 1 end extension in this run (`end_extension/`)

2,299 of 261,803 loci moved (1,764 at the 5' end, 535 at the 3' end). The shortie groups (consensus_28seqs, consensus_47seqs) did not move: their consensuses now carry the tail. Groups that moved:

| group | side | windows | extended | median bp | flag |
|---|---|---|---|---|---|
| consensus_10seqs_noterm | 5' | 7,830 | 1,233 | 54 | |
| r10_19seqs | 3' | 1,913 | 508 | 45 | >= 20 % of the group grew |
| r1_r3_P1_hla020 | 5' | 3,830 | 454 | 37 | |
| r5_r6_short_hla | 5' | 2,885 | 58 | 39 | |
| MEG-RS | 3' / 5' | 515 | 25 / 19 | 60 / 38 | |
| r10_r8_P18del | 3' | 732 | 2 | 32 | |

## Contents
| Path | What |
|---|---|
| `report.html` | full report |
| `loci_verdict.tsv.gz`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |
| `family_view/` | family roll-up per copy, partial copies, stretches, the composite stretches searched in step 1 |
| `alignments/` | plates (SINE loci, topped up with grey copies marked [grey]) |
| `alignments/<sp>_all_copies_subfam30k.aln.fa` | step 1 SubFam: 600 chunk consensuses of 30,000 random copies of all subfamilies, aligned with the bank consensuses |
| `end_extension/` | step 1 end extension: per group and side (`ends_extension.tsv`), every moved locus (`loci_extended.tsv`) |
| `compare_to_hla_rsibank/` | per-copy cross-tab against hla_rsibank_20261026 (`crosstab.tsv`, `by_new.tsv`, `by_old.tsv`, `crosstab.sh`) |
| `continuation/`, `flankscan/` | continuation check; flank scan stages 1-7 |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers |
| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | flank twins, tandem-array share, satellite screen |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its copies |
| `composite_units.tsv`, `consensus_blocks.tsv` | composites read as units; the bank as blocks |
| `consensuses.clean.fa`, `manifest.txt`, `pipeline_stages.txt` | the bank; tool versions and command; stage times |
