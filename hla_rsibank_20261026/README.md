# hla_rsibank_20261026

*Hipposideros larvatus* (hla), assembly GCA_031876335.1, searched with the *R. sinicus* bank. Bank: the *Rhinolophus sinicus* bank without r1_9seqs (21 consensuses: r2-r10, r10_r8_P18del, seven composites, four MEG SINEs; therioserver `/home/toki/rhin/final_bank_20261026/bank21_nor1.fa`), families Rhin and MEG (`families.tsv`). SINEderella main a9d2a51 (2026-10-26): composite stretches of >= 50 bp searched as queries of their own in step 1; partial copies labelled and grouped into stretches in the family roll-up (stretch tables regrouped with main 9d484fa); flank scan stages 1-7 and the report rebuild as for rsi_main_20261022.

- Server: therioserver `/home/toki/rhin/hla_rsibank_20261026/` (`../final_bank_20261026/launch_final.sh`, `pack_final.sh`).
- Run: start 2026-10-26 20:30:20 hla; done 2026-10-27 01:11:09.

## Verdicts (every copy; `loci_verdict.tsv.gz`, rules SINEderella `docs/VERDICTS.md`)

| | copies |
|---|---|
| copies in the list | 233,325 |
| assigned to a subfamily | 191,014 |
| SINE loci (all copies) | 144,006 |
| grey zone (all copies) | 21,012 |
| not SINE (all copies) | 68,307 |

Copies without a subfamily are listed as not SINE in this table, partial copies included (below).

Per family (copies assigned to a subfamily):

| family | subfamilies with copies | assigned | SINE | grey | not SINE |
|---|---|---|---|---|---|
| MEG | 4 | 562 | 121 | 10 | 431 |
| Rhin | 17 | 190,452 | 143,885 | 21,002 | 25,565 |

## Family roll-up and partial copies (`family_view/`)

| family | copies | subfamily-firm | near tie | family only | not firm |
|---|---|---|---|---|---|
| MEG | 639 | 562 | 20 | 77 | 0 |
| Rhin | 232,686 | 190,452 | 6,848 | 42,203 | 31 |

Partial copies (family-firm, covering < 60 % of every competing consensus): 18,454, plus 2,055 below their stretch bar; 118 stretches (`partial_stretches.tsv`, `partial_copies.tsv`). Largest stretches:

| stretch | copies | median length | spans |
|---|---|---|---|
| Rhin_partial1 | 15,022 | 86 | r5_27seqs:1-86, r5_r6_P26:1-81 |
| Rhin_partial2 | 449 | 100 | r5_27seqs:1-100 |
| Rhin_partial3 | 384 | 112 | r10_r8_P18:290-401 |
| Rhin_partial4 | 309 | 127 | r10_r8_P18:158-284 |
| Rhin_partial5 | 157 | 234 | r10_r8_P18:1-237 |
| Rhin_partial6 | 146 | 75 | r5_27seqs:12-86, r5_r5_P48:12-86 |
| Rhin_partial7 | 145 | 218 | r10_r8_P18:172-384 |
| Rhin_partial8 | 129 | 93 | r10_r8_P18:292-384 |

## Genome and the earlier run
*Hipposideros larvatus* GCA_031876335.1 (mHipLar1.pri), therioserver `/home/toki/chiro/genomes/hla.fna`. The earlier run
(`../hla_main_20261024/`) searched the Rhin-1 + VES bank of the Chiroptera page: 71,135 copies, 13,597 SINE loci. This run
uses the *R. sinicus* bank: 233,325 copies, most of them assigned to r10_r8_P18 (140,029). The two runs answer different
questions and are not a before/after pair.
The largest partial stretch (r5_27seqs 1-86 = r5_r6_P26 1-81, about 15,000 copies) is a head that matches the first 86 bp
of the r5 head and stops: a family of this genome whose head resembles r5 and whose body the bank does not hold. Its
copies are family-firm Rhin but have no subfamily; a consensus built from them would answer what they are.

## Contents
| Path | What |
|---|---|
| `report.html` | full report |
| `loci_verdict.tsv.gz`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |
| `family_view/` | family roll-up per copy, partial copies, stretches, the composite stretches searched in step 1 |
| `alignments/` | plates (SINE loci, topped up with grey copies marked [grey]) |
| `alignments/hla_all_copies_subfam30k.aln.fa` | step 1 SubFam: 600 chunk consensuses of 30,000 random copies of all subfamilies, aligned with the 21 bank consensuses |
| `continuation/`, `flankscan/` | continuation check; flank scan stages 1-7 |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers |
| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | flank twins, tandem-array share, satellite screen |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its copies |
| `composite_units.tsv`, `consensus_blocks.tsv` | composites read as units; the bank as blocks |
| `consensuses.clean.fa`, `manifest.txt`, `pipeline_stages.txt` | the bank; tool versions and command; stage times |

## r10_r8_P18 in this genome is not rsi's P18 (2026-10-27)
The consensus audit (`consensus_audit/summary.tsv`) calls r10_r8_P18 DIVERGED here: rebuilt from its 140,029 hla copies
it is 440-442 bp against 426, with 66-70 mismatches and 106-110 gap columns, the same with both seeds. On the subfamily
(chunk) plate every one of the 201 chunk consensuses carries an insertion of 18-53 bp after P18 position ~386-401, so this
is the hla element's own 3' part, not a subfamily mark. Rarer: 16-20 bp inside the 10-bp-unit middle piece (after ~209,
2 of 201 chunks; 9 of the top 100 copies) and 107-140 bp before the start (2 chunks).
`alignments/hla_r10_r8_P18_peel.aln.fa`: the subfamily plate (SubFam chunk order kept, bases unchanged) with the two
rebuilt hla consensuses added under rsi's P18 (MAFFT --add, therioserver `/home/toki/tmp/hla_p18_peel/`), for peeling.
