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

## Shorties are complete short elements (`shortie_check/`)

consensus_47seqs (128 bp) and consensus_28seqs (134 bp) are the head of the long elements (consensus_10seqs_b, consensus_29seqs, consensus_12seqs: first ~150 bp at 96-99 %, with an ~18 bp deletion in the shorties). Their genome copies stop where the consensus stops:

| test | result |
|---|---|
| shared 10-mers past either end, all copies (step 1 end extension) | 0 (42,402 windows for 47seqs, 11,191 for 28seqs) |
| 150 bp 3' flank of 100 random SINE copies each, ssearch36 against the long-element bodies (consensus position 150 to end), E < 1e-3 | 47seqs 2 of 100; 28seqs 0 of 100 |
| control: 150 bp after the head in 100 random SINE copies each of consensus_10seqs_b and consensus_29seqs, same search | 99 of 100; 100 of 100 |
| verdicts of assigned copies | 47seqs SINE 30,934, grey 6,661, not SINE 4,141; 28seqs SINE 8,998, grey 1,880, not SINE 626 |

Plate: `shortie_check/hla_shorties_flank150.aln.fa`: the 5 consensuses and 100 random SINE copies of each shortie with 150 bp genome flanks (copy upper case, flanks lower case, element orientation), mafft localpair; script `shortie_plate.sh`.

## Heads against tRNA genes (`trna_plates/`)

One plate per consensus: the consensus against its best-scoring human tRNA gene (GtRNAdb hg19), aligned by the same ssearch36
that scored it in the workflow's tRNA step (SINEderella `tools/trna_heads.sh` + `tools/trna_plates.py`). Aligned part in upper case,
the rest in lower case. "Shuffles" = best score of 20 dinucleotide shuffles of the consensus. Consensus names link to the plate in ViewAlign; the chiroptera page has a picker too.

| consensus | best tRNA gene | identity | aligned | consensus span | bits | shuffles max |
|---|---|---|---|---|---|---|
| [MEG-RL](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2FMEG-RL.trna.aln.fa&title=hla%20MEG-RL%20vs%20tRNA-Arg-TCT-1-1) | tRNA-Arg-TCT-1-1 | 62.5 % | 40 nt | 86-125 | 20.4 | 24.9 |
| [MEG-RS](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2FMEG-RS.trna.aln.fa&title=hla%20MEG-RS%20vs%20tRNA-Arg-TCT-1-1) | tRNA-Arg-TCT-1-1 | 62.5 % | 40 nt | 86-125 | 19.7 | 25.3 |
| [MEG-T2](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2FMEG-T2.trna.aln.fa&title=hla%20MEG-T2%20vs%20tRNA-Val-TAC-1-1) | tRNA-Val-TAC-1-1 | 94.5 % | 74 nt | 6-79 | 80.3 | 24.3 |
| [MEG-TR](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2FMEG-TR.trna.aln.fa&title=hla%20MEG-TR%20vs%20tRNA-Val-TAC-1-1) | tRNA-Val-TAC-1-1 | 94.5 % | 55 nt | 1-55 | 62.5 | 25.7 |
| [consensus_10seqs_b](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fconsensus_10seqs_b.trna.aln.fa&title=hla%20consensus_10seqs_b%20vs%20tRNA-Ile-GAT-1-2) | tRNA-Ile-GAT-1-2 | 76.1 % | 71 nt | 255-321 | 39.9 | 27.0 |
| [consensus_10seqs_noterm](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fconsensus_10seqs_noterm.trna.aln.fa&title=hla%20consensus_10seqs_noterm%20vs%20tRNA-Ile-GAT-1-2) | tRNA-Ile-GAT-1-2 (tied: tRNA-Lys-CTT-4-1, 41.6 bits) | 77.6 % | 69 nt | 255-321 | 43.0 | 24.8 |
| [consensus_12seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fconsensus_12seqs.trna.aln.fa&title=hla%20consensus_12seqs%20vs%20tRNA-Ile-GAT-1-2) | tRNA-Ile-GAT-1-2 | 78.3 % | 71 nt | 239-307 | 46.9 | 26.1 |
| [consensus_28seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fconsensus_28seqs.trna.aln.fa&title=hla%20consensus_28seqs%20vs%20tRNA-Ile-GAT-1-2) | tRNA-Ile-GAT-1-2 | 78.6 % | 82 nt | 6-87 | 34.1 | 24.2 |
| [consensus_29seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fconsensus_29seqs.trna.aln.fa&title=hla%20consensus_29seqs%20vs%20tRNA-Ile-GAT-1-2) | tRNA-Ile-GAT-1-2 | 78.6 % | 71 nt | 124-193 | 50.5 | 25.2 |
| [consensus_386seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fconsensus_386seqs.trna.aln.fa&title=hla%20consensus_386seqs%20vs%20tRNA-Pro-AGG-1-1) | tRNA-Pro-AGG-1-1 (tied: tRNA-Ala-AGC-4-1, 31.3 bits) | 77.8 % | 62 nt | 14-67 | 32.0 | 24.6 |
| [consensus_47seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fconsensus_47seqs.trna.aln.fa&title=hla%20consensus_47seqs%20vs%20tRNA-Lys-TTT-5-1) | tRNA-Lys-TTT-5-1 (tied: tRNA-Ser-GCT-5-1, 28.6 bits) | 72.2 % | 79 nt | 6-84 | 28.9 | 23.0 |
| [r10_19seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr10_19seqs.trna.aln.fa&title=hla%20r10_19seqs%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 79.7 % | 74 nt | 1-74 | 62.8 | 24.3 |
| [r10_r8_P18del](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr10_r8_P18del.trna.aln.fa&title=hla%20r10_r8_P18del%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 79.7 % | 74 nt | 6-79 | 60.7 | 26.1 |
| [r1_r2_r3_r3_C11](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr1_r2_r3_r3_C11.trna.aln.fa&title=hla%20r1_r2_r3_r3_C11%20vs%20tRNA-Ile-AAT-8-1) | tRNA-Ile-AAT-8-1 | 75.7 % | 75 nt | 6-80 | 50.8 | 25.4 |
| [r1_r3_P1](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr1_r3_P1.trna.aln.fa&title=hla%20r1_r3_P1%20vs%20tRNA-Ile-AAT-8-1) | tRNA-Ile-AAT-8-1 | 74.3 % | 75 nt | 6-80 | 48.7 | 27.9 |
| [r1_r3_P1_hla020](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr1_r3_P1_hla020.trna.aln.fa&title=hla%20r1_r3_P1_hla020%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 76.7 % | 75 nt | 6-79 | 49.0 | 28.5 |
| [r1_r3_P2](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr1_r3_P2.trna.aln.fa&title=hla%20r1_r3_P2%20vs%20tRNA-Ile-AAT-8-1) | tRNA-Ile-AAT-8-1 | 74.3 % | 75 nt | 6-80 | 49.4 | 27.7 |
| [r2_3seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr2_3seqs.trna.aln.fa&title=hla%20r2_3seqs%20vs%20tRNA-Ile-AAT-2-1) | tRNA-Ile-AAT-2-1 | 74.2 % | 67 nt | 1-66 | 45.9 | 26.8 |
| [r3_58seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr3_58seqs.trna.aln.fa&title=hla%20r3_58seqs%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 74.6 % | 68 nt | 1-67 | 46.1 | 26.5 |
| [r4_32seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr4_32seqs.trna.aln.fa&title=hla%20r4_32seqs%20vs%20tRNA-Ile-AAT-2-1) | tRNA-Ile-AAT-2-1 | 75.8 % | 67 nt | 1-66 | 47.6 | 22.6 |
| [r5_27seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr5_27seqs.trna.aln.fa&title=hla%20r5_27seqs%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 78.1 % | 74 nt | 5-77 | 54.6 | 27.5 |
| [r5_r3_P34](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr5_r3_P34.trna.aln.fa&title=hla%20r5_r3_P34%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 75.7 % | 71 nt | 138-207 | 48.9 | 26.5 |
| [r5_r5_P48](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr5_r5_P48.trna.aln.fa&title=hla%20r5_r5_P48%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 78.1 % | 74 nt | 5-77 | 55.0 | 25.5 |
| [r5_r5_P48_hla007](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr5_r5_P48_hla007.trna.aln.fa&title=hla%20r5_r5_P48_hla007%20vs%20tRNA-Ile-GAT-1-2) | tRNA-Ile-GAT-1-2 | 78.6 % | 71 nt | 124-193 | 51.4 | 24.1 |
| [r5_r6_long_hla](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr5_r6_long_hla.trna.aln.fa&title=hla%20r5_r6_long_hla%20vs%20tRNA-Ile-GAT-1-2) | tRNA-Ile-GAT-1-2 | 81.2 % | 71 nt | 122-190 | 53.7 | 25.1 |
| [r5_r6_short_hla](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr5_r6_short_hla.trna.aln.fa&title=hla%20r5_r6_short_hla%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 81.4 % | 71 nt | 6-75 | 58.5 | 23.8 |
| [r7_133seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr7_133seqs.trna.aln.fa&title=hla%20r7_133seqs%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 86.5 % | 74 nt | 5-78 | 79.6 | 23.7 |
| [r8_hla](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr8_hla.trna.aln.fa&title=hla%20r8_hla%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 87.8 % | 74 nt | 6-79 | 76.5 | 24.5 |
| [r9_15seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr9_15seqs.trna.aln.fa&title=hla%20r9_15seqs%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 93.2 % | 74 nt | 5-78 | 86.0 | 24.0 |
| [r9_ext](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fhla_hlabank_20261030%2Ftrna_plates%2Fr9_ext.trna.aln.fa&title=hla%20r9_ext%20vs%20tRNA-Ile-AAT-5-3) | tRNA-Ile-AAT-5-3 | 90.5 % | 74 nt | 6-79 | 83.9 | 26.9 |

MEG-RL and MEG-RS are 5S-derived (Rfam 5S_rRNA 1-119); their best tRNA hit is below the shuffles. When another amino acid's best gene scores within 2 bits
of the best, the isotype is tied and both genes are on the plate (consensus_10seqs_noterm Ile/Lys, consensus_386seqs Pro/Ala,
consensus_47seqs Lys/Ser). Page: https://toki-bio.github.io/Tal/hla_hlabank_20261030/trna_plates/

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
| `trna_plates/` | each consensus against its best tRNA gene (ViewAlign plates) + `index.tsv` |
