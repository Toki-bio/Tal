# hla peel of the step 1 SubFam plate (round 1, 2026-10-29)

Plate: `../alignments/hla_all_copies_subfam30k.aln.fa` (600 chunk consensuses of 30,000 random copies, aligned with
the 21 rsi bank consensuses; SubFam's order). A number such as 335 means chunk `input_335`; "row" is the plate row.

## Owner's calls (verbatim in substance)

| rows | chunks | call |
|---|---|---|
| top 6 | MEG-RL, MEG-RS, MEG-TR, MEG-T2, 001, 002 | MEG, deleted |
| bottom 15 | 003, 600 | noise |
| | 463 | hard to decide, skipped |
| | r10_19seqs, r10_r8_P18, r10_r8_P18del | no good matches here |
| | 511, 462 | a bit similar to r10, singletons (not enough to support a subfamily) |
| 607-613 | r1_r3_P2 ... r1_r2_r3_r3_C11 (6 consensuses) | no matches, removed |
| 590-606 | 513 ... 478, with r9_15seqs | two forms of r9: short as in rsi (513, 514; 525 maybe) and longer (514-519); 501, 478 unclear (deletion). Keep r9 as in rsi, add **r9_ext** |
| 558-589 | 471 ... 512, with r8_83seqs, r7_133seqs | r8 with a difference in ttcccc near the tail; 473 outlier (insertion). Record **r8_hla** |
| 533-557 | 008 ... r6_210seqs | closest r5_r6; more than half carry the insertion ATGTTCCCCTTCCCCAATTCCCCAACACTGTTCCCTTCCC: **r5_r6_short_hla**, **r5_r6_long_hla** (014, 015, 017: outliers with an extra insertion) |
| 7-532 | the remaining 526 | no rsi consensus; mosaic, hard to peel |
| 454-532 | 321 ... 079 (79) | shorties; 081, 056, 007, 079 outliers, deleted |
| 501-528 | 526 ... 598 (28) | vary by presence of the tail CgCTGAGCTGCCAAcagGtGGcAGGATGGCTCAGCT; "need to clarify from copies": **consensus_28seqs** |
| 454-500 | 321 ... 347 (47) | same tail variation: **consensus_47seqs** |
| 393-453 | 006 ... 068 | next batch (owner continues) |

The seven consensuses: `peel_round1_owner_consensuses.fa`.

## Tail of the shorties, from the copies

`tail_check.py`: every copy of the 75 shortie chunks (3,750 copies), extended 50 bp 5' and 300 bp 3' in the genome
(strand from the copy name), aligned with ssearch36 -z 11 to consensus_28seqs, consensus_47seqs and r5_r6_long_hla
(`tail_check.tsv`).

| group | copies | full tail (end >= 128 of 134) | partial tail | no tail | genome copy runs > 10 bp past the step 1 hit | body (r8_hla / r5_r6_short_hla / r9_ext) downstream |
|---|---|---|---|---|---|---|
| 28 chunks (526-598) | 1,400 | 1,315 (94 %) | 62 | 23 | 764 | 60 |
| 47 chunks (321-347) | 2,350 | 2,308 (98 %) | 32 | 10 | 708 | 199 |

The tail is in nearly every genome copy. Where a chunk consensus lacks it, the step 1 hit stopped before the tail (the
bank has no consensus for this element), so the absence on the plate is a hit-boundary effect, not a variant.

What follows the tail: 200 random copies per group without a body downstream, from the start of the alignment to
120 bp past its end (`hla_shortie_G28_copies_tail.aln.fa`, `hla_shortie_G47_copies_tail.aln.fa`; copies whose shown
flank shares a 20-mer with another copy dropped: 4 of 200 and 0 of 200; the owner's consensus on top). Column
consensus (columns occupied >= 50 %):

- 28: `...CCAACAGGTGGCAGGATGGCTCAGCT` + 28 A + then AT-rich columns with plurality share 0.25-0.41 (flank)
- 47: `...CCAACAGGTAGCAGGATGGCTCAGCT` + 31 A + then 0.24-0.34 (flank)

So both shorties end right after GGCTCAGCT with an A tail; the sequence after it is unique flank.

## Chunk consensuses against the bank plus the seven new consensuses

`chunks_best.py` -> `chunks_best.tsv`: best consensus per chunk (ssearch36 -z 11), identity, spans, second best.

- rows 7-392 (and other P18-like rows; 384 chunks in all): best r10_r8_P18, mean 80.9 % identity, 369 cover >= 90 % of it
  (the P18 that the consensus audit calls DIVERGED in hla).
- next batch, rows 394-419 (022 ... 021, 26 chunks of 412-420 bp): r5_r6_long_hla, 94-96 % over all 393 bp; the chunk
  aligns over 1-410, ~17 bp longer. Row 393 (006, 312 bp) is the 3' part (cons 119-393).
- rows 420-453 except a few: 510-565 bp chunks, best r1_r2_r3_r3_C11 at 75-82 % over its 81-608 (058: r1_r3_P1
  81.6 % over all 448); second best r5_r6_long_hla at 0.84-0.97 of the score. Rows 422, 428, 430, 451 (048, 080, 082,
  077) match only the first ~250-265 bp of r5_r6_long_hla.

## Round 2 (owner, 2026-10-29): rows 393-453 and the big block

Consensuses: `peel_round2_owner_consensuses.fa`. The owner named two of them consensus_10seqs; here they are
consensus_10seqs_noterm (069-068) and consensus_10seqs_b (048-076 group).

| chunks | call |
|---|---|
| 069 ... 068 (10) | good match, **without terminator**; the two bottom ones odd (large insertion and deletion): **consensus_10seqs_noterm** |
| 048 ... 076 | with terminator, except the odd 080, 597, 082. AGAGCACGTGCTCC splits them into groups, and CACGAGCTGCCGACGGGAAGCCAGATGGCGAC marks the lower group; all else matches: **consensus_12seqs**, **consensus_10seqs_b** |
| 006 ... 027 (029 rows) | **consensus_29seqs** (026, 027 carry the long insertion of the previous group) |
| 004 ... 005 (386) | could be subdivided, maybe artefacts; proposed as one unresolved group: **consensus_386seqs** |

consensus_386seqs against the bank and the consensus audit (ssearch36): r10_r8_P18_seed1 (the hla P18 rebuilt by
`../consensus_audit/` from its copies) 99.1 % over 1-438 / 1-437; seed2 99.1 %; rsi r10_r8_P18 80.3 % over 1-354 / 1-384.
So the big block is hla's own P18, the one the audit already rebuilt and called DIVERGED.
