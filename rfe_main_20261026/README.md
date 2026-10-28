# rfe_main_20261026

*Rhinolophus ferrumequinum* (rfe), assembly GCF_004115265.2, searched with the *R. sinicus* bank. Bank: the *Rhinolophus sinicus* bank without r1_9seqs (21 consensuses: r2-r10, r10_r8_P18del, seven composites, four MEG SINEs; therioserver `/home/toki/rhin/final_bank_20261026/bank21_nor1.fa`), families Rhin and MEG (`families.tsv`). SINEderella main a9d2a51 (2026-10-26): composite stretches of >= 50 bp searched as queries of their own in step 1; partial copies labelled and grouped into stretches in the family roll-up (stretch tables regrouped with main 9d484fa); flank scan stages 1-7 and the report rebuild as for rsi_main_20261022.

- Server: therioserver `/home/toki/rhin/rfe_main_20261026/` (`../final_bank_20261026/launch_final.sh`, `pack_final.sh`).
- Run: start 2026-10-26 20:30:20 rfe; done 2026-10-26 22:16:14.

## Verdicts (every copy; `loci_verdict.tsv.gz`, rules SINEderella `docs/VERDICTS.md`)

| | copies |
|---|---|
| copies in the list | 63,089 |
| assigned to a subfamily | 58,105 |
| SINE loci (all copies) | 49,125 |
| grey zone (all copies) | 7,049 |
| not SINE (all copies) | 6,915 |

Copies without a subfamily are listed as not SINE in this table, partial copies included (below).

Per family (copies assigned to a subfamily):

| family | subfamilies with copies | assigned | SINE | grey | not SINE |
|---|---|---|---|---|---|
| MEG | 4 | 123 | 117 | 2 | 4 |
| Rhin | 17 | 57,982 | 49,008 | 7,047 | 1,927 |

## Family roll-up and partial copies (`family_view/`)

| family | copies | subfamily-firm | near tie | family only | not firm |
|---|---|---|---|---|---|
| MEG | 143 | 123 | 3 | 20 | 0 |
| Rhin | 62,946 | 57,982 | 2,506 | 4,950 | 14 |

Partial copies (family-firm, covering < 60 % of every competing consensus): 1,291, plus 286 below their stretch bar; 93 stretches (`partial_stretches.tsv`, `partial_copies.tsv`). Largest stretches:

| stretch | copies | median length | spans |
|---|---|---|---|
| Rhin_partial1 | 371 | 134 | r10_r8_P18:293-426 |
| Rhin_partial2 | 97 | 73 | r10_r8_P18del:1-79, r10_19seqs:2-74 |
| Rhin_partial3 | 83 | 121 | r1_r3_P1:308-429, r1_r3_P2:165-286 |
| Rhin_partial4 | 64 | 265 | r1_r2_r3_r3_C11:335-607, r1_r3_P1:201-447 |
| Rhin_partial5 | 60 | 143 | r1_r3_P1:288-429, r1_r3_P2:145-286 |
| Rhin_partial6 | 54 | 199 | r1_r3_P1:1-199, r1_r2_r3_r3_C11:1-192 |
| Rhin_partial7 | 50 | 123 | r1_r2_r3_r3_C11:1-130, r1_r3_P1:1-121, r1_r3_P2:1-114 |
| Rhin_partial8 | 42 | 99 | r1_r3_P1:332-429, r1_r3_P2:189-286 |

## Genome
*Rhinolophus ferrumequinum* RefSeq GCF_004115265.2 (mRhiFer1_v1.p, chromosome level, 2,075,768,562 bp in 134 sequences),
downloaded 2026-10-26 and checked against NCBI's MD5 (therioserver `/home/toki/rhin/genomes/rfe.fna`). First SINEderella
run on this species; searched with the *R. sinicus* bank as it stands, so families of its own that the bank lacks show up
only as weak hits, partial copies or continuations.

## Contents
| Path | What |
|---|---|
| `report.html` | full report |
| `loci_verdict.tsv.gz`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |
| `family_view/` | family roll-up per copy, partial copies, stretches, the composite stretches searched in step 1 |
| `alignments/` | plates (SINE loci, topped up with grey copies marked [grey]) |
| `alignments/rfe_all_copies_subfam30k.aln.fa` | step 1 SubFam: 600 chunk consensuses of 30,000 random copies of all subfamilies, aligned with the 21 bank consensuses |
| `continuation/`, `flankscan/` | continuation check; flank scan stages 1-7 |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers |
| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | flank twins, tandem-array share, satellite screen |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its copies |
| `composite_units.tsv`, `consensus_blocks.tsv` | composites read as units; the bank as blocks |
| `consensuses.clean.fa`, `manifest.txt`, `pipeline_stages.txt` | the bank; tool versions and command; stage times |
