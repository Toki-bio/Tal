# rsi: genome copies of every elementary unit, standalone and cut out of composites (2026-11-01)

Run `rsi_main_20261026` (bank21_nor1). Groups: 50 random SINE-verdict copies (unique flanks) per standalone consensus, used whole
(`<unit>@alone`); and 50 random SINE-verdict copies per composite, each aligned to the composite consensus (ssearch36, BTOP) and
cut at its own junctions, i.e. the composite stretch of `composite_units.tsv` projected onto the copy (`<unit>@<composite>:a-b`,
kept when the alignment covers >= 80 % of the stretch). r10_r8_P18del is cut like P18: r10 head 1-151, second unit 152-255.
P18's second unit (290-426) is called `P18_2nd` here. Scripts: `unit_pieces_copies.py`, `plates.sh` (therioserver
`/home/toki/tmp/r8units/`). Plates: MAFFT --localpair --maxiterate 1000, input order: consensus rows first, then the copies
group by group; copy names end with the group.

## r8 alone vs the second unit of P18 / P18del

Median identity of each group's pieces to each consensus (best local hit), and how many pieces score best against it:

| group | pieces | r8_83seqs | r7_133seqs | r10_19seqs | P18 290-426 | best r8 / r7 / r10 / P18_2nd |
|---|---|---|---|---|---|---|
| r8_83seqs@alone | 50 | 91.1 | 87.5 | 79.5 | 75.9 | 49 / 0 / 1 / 0 |
| r8_83seqs@r1_r3_P2:2-85 | 50 | 80.1 | 77.1 | 86.1 | 65.7 | 0 / 0 / 50 / 0 |
| P18_2nd@r10_r8_P18:290-426 | 48 | 78.2 | 77.7 | 71.9 | 90.3 | 0 / 0 / 0 / 48 |
| P18_2nd@r10_r8_P18del:152-255 | 48 | 78.1 | 77.0 | 71.8 | 91.5 | 0 / 0 / 0 / 48 |

Plate: [r8 question](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr8_question.aln.fa&title=rsi%20r8%20alone%20vs%20P18%20second%20unit%20%28copies%20cut%20at%20junctions%29) (r7 and r10 consensuses on top, then r8, the P2 "r8" piece and the two P18_2nd groups).

## Every group against every standalone unit

Each piece scored against the standalone consensuses (r2-r10) and P18 290-426; the unit that scores best, and the median identity to it:

| group | pieces | best unit (pieces) | median identity to best |
|---|---|---|---|
| r10_19seqs@alone | 50 | r10_19seqs 50 | 90.0 |
| r10_r8_P18del@alone | 50 | r10_19seqs 35, P18_2nd 15 | 90.7 |
| r2_3seqs@alone | 19 | r2_3seqs 19 | 69.8 |
| r3_58seqs@alone | 16 | r3_58seqs 16 | 92.8 |
| r4_32seqs@alone | 11 | r4_32seqs 10, r3_58seqs 1 | 70.5 |
| r5_27seqs@alone | 50 | r5_27seqs 50 | 92.3 |
| r6_210seqs@alone | 50 | r6_210seqs 49, r5_27seqs 1 | 92.3 |
| r7_133seqs@alone | 50 | r7_133seqs 50 | 92.9 |
| r8_83seqs@alone | 50 | r8_83seqs 49, r10_19seqs 1 | 91.1 |
| r9_15seqs@alone | 50 | r9_15seqs 50 | 90.4 |
| r10_19seqs@r10_r8_P18:6-151 | 49 | r10_19seqs 49 | 91.0 |
| unassigned@r10_r8_P18:152-289 | 48 | r10_19seqs 30, none 8, r7_133seqs 5, r9_15seqs 3, r4_32seqs 2 | 83.8 |
| P18_2nd@r10_r8_P18:290-426 | 48 | P18_2nd 48 | 90.3 |
| r10_r8_P18del@r1_r2_r3_r3_C11:1-80 | 42 | r10_19seqs 42 | 88.4 |
| r2_3seqs@r1_r2_r3_r3_C11:81-338 | 49 | r2_3seqs 48, r3_58seqs 1 | 94.2 |
| r4_32seqs@r1_r2_r3_r3_C11:339-471 | 50 | r3_58seqs 48, r4_32seqs 2 | 95.1 |
| r3_58seqs@r1_r2_r3_r3_C11:472-589 | 46 | r3_58seqs 43, r6_210seqs 2, r8_83seqs 1 | 90.1 |
| r10_r8_P18del@r1_r3_P1:1-80 | 47 | r10_19seqs 47 | 87.8 |
| r2_3seqs@r1_r3_P1:81-200 | 50 | r2_3seqs 49, r4_32seqs 1 | 92.2 |
| r3_58seqs@r1_r3_P1:201-401 | 50 | r3_58seqs 50 | 93.8 |
| unassigned@r1_r3_P1:402-448 | 45 | r5_27seqs 31, r8_83seqs 9, r6_210seqs 4, r3_58seqs 1 | 95.8 |
| r8_83seqs@r1_r3_P2:2-85 | 50 | r10_19seqs 50 | 86.1 |
| r3_58seqs@r1_r3_P2:86-258 | 50 | r3_58seqs 50 | 94.7 |
| r5_27seqs@r1_r3_P2:259-294 | 47 | r5_27seqs 34, r8_83seqs 9, r6_210seqs 2, r9_15seqs 1, r3_58seqs 1 | 96.0 |
| r5_27seqs@r5_r3_P34:2-140 | 47 | r5_27seqs 46, r8_83seqs 1 | 82.1 |
| r3_58seqs@r5_r3_P34:141-339 | 49 | r3_58seqs 45, r6_210seqs 3, r8_83seqs 1 | 91.8 |
| r5_27seqs@r5_r3_P34:340-370 | 41 | r5_27seqs 36, r3_58seqs 4, r4_32seqs 1 | 95.7 |
| r5_27seqs@r5_r5_P48:1-136 | 46 | r5_27seqs 31, r8_83seqs 14, r6_210seqs 1 | 89.0 |
| r6_210seqs@r5_r5_P48:137-308 | 49 | r8_83seqs 27, r5_27seqs 17, r6_210seqs 4, r7_133seqs 1 | 88.9 |
| r5_27seqs@r5_r6_P26:2-137 | 47 | r5_27seqs 45, r8_83seqs 2 | 88.1 |
| r6_210seqs@r5_r6_P26:138-376 | 49 | r6_210seqs 49 | 92.0 |
| r10_19seqs@r10_r8_P18del:1-151 | 48 | r10_19seqs 48 | 89.4 |
| P18_2nd@r10_r8_P18del:152-255 | 48 | P18_2nd 48 | 91.5 |

## Plates per unit label (all groups carrying that label in `composite_units.tsv`)

- [P18_2nd](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2FP18_2nd.aln.fa&title=rsi%20unit%20pieces%3A%20P18_2nd%20%28alone%20%2B%20cut%20from%20composites%29)
- [r10_19seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr10_19seqs.aln.fa&title=rsi%20unit%20pieces%3A%20r10_19seqs%20%28alone%20%2B%20cut%20from%20composites%29)
- [r10_r8_P18del](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr10_r8_P18del.aln.fa&title=rsi%20unit%20pieces%3A%20r10_r8_P18del%20%28alone%20%2B%20cut%20from%20composites%29)
- [r2_3seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr2_3seqs.aln.fa&title=rsi%20unit%20pieces%3A%20r2_3seqs%20%28alone%20%2B%20cut%20from%20composites%29)
- [r3_58seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr3_58seqs.aln.fa&title=rsi%20unit%20pieces%3A%20r3_58seqs%20%28alone%20%2B%20cut%20from%20composites%29)
- [r4_32seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr4_32seqs.aln.fa&title=rsi%20unit%20pieces%3A%20r4_32seqs%20%28alone%20%2B%20cut%20from%20composites%29)
- [r5_27seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr5_27seqs.aln.fa&title=rsi%20unit%20pieces%3A%20r5_27seqs%20%28alone%20%2B%20cut%20from%20composites%29)
- [r6_210seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr6_210seqs.aln.fa&title=rsi%20unit%20pieces%3A%20r6_210seqs%20%28alone%20%2B%20cut%20from%20composites%29)
- [r7_133seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr7_133seqs.aln.fa&title=rsi%20unit%20pieces%3A%20r7_133seqs%20%28alone%20%2B%20cut%20from%20composites%29)
- [r8_83seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr8_83seqs.aln.fa&title=rsi%20unit%20pieces%3A%20r8_83seqs%20%28alone%20%2B%20cut%20from%20composites%29)
- [r9_15seqs](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Fr9_15seqs.aln.fa&title=rsi%20unit%20pieces%3A%20r9_15seqs%20%28alone%20%2B%20cut%20from%20composites%29)
- [unassigned](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_main_20261026%2Funit_pieces_copies%2Fplates%2Funassigned.aln.fa&title=rsi%20unit%20pieces%3A%20unassigned%20%28alone%20%2B%20cut%20from%20composites%29)
