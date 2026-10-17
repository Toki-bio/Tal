# rsi: unknown partners from one-sided shared flanks (flankscan stage 9b, draft), 2026-10-17

SINEderella branch `fs9b-partners` (cbf1348), `flankscan/fs9b_partners.sh`, run on this rsi run with the bat partner
library (therioserver `~/refs/partnerlib/partnerlib_9397.fa`). Stage 9 links copies whose flanks are similar from the
junction outward; its body check marks pairs with independent bodies "flank_only" (hot-spot signature) and drops them.
9b keeps flank_only groups shared on one side, removes members that are duplicated loci (other side shared too, or a
stage-9 twin / array), keeps groups with >= 10 copies on >= 3 contigs after that, builds the
partner, names it if the bank or the partner library knows it, counts copies of every family carrying it at the junction.

`unknown_partners.summary.tsv`: one row per distinct partner. kind = known partner (bank / library) | family-specific:
composite candidate (>= 50 % of the carriers in one family and >= 2x that family's share of all copies) | all families:
background. `unknown_partners.tsv`: one row per group.

Review alignments of the family-specific candidates (`<id>.review.aln.fa`): first row = the candidate composite (partner +
family consensus); members laid on its positions (MAFFT --addfragments --keeplength), outer flanks (100 bp) lower case,
degapped, packed against each member's own first / last aligned base; `ins=` = bases a member has beyond the composite.
`r5_r6_P26_5p1` is the X group found by eye: 21 copies after removing 9 duplicated members (among them the block of 6 on
NC_142512 that shares both flanks and the body); 11 of the 21 are hand-found copies.

## X identified (2026-10-17)
Searched against Dfam curated families (current release, 30,646 families; the online Dfam search queue did not move for
3 h, so locally: therioserver `/home/toki/refs/dfam_current/`, ssearch36 -z 11, E <= 1e-3). The X partner (168 bp,
`r5_r6_P26_5p1.partner.fa`) = the 3' part of **RhiSin-1.14** (DF000280481, *R. sinicus*, 237 bp, **DNA/TcMar-Tc2**, Zoonomia):
95.5 % over 161 bp, X 8-168 = family 82-237 (to its end), same orientation as P26. Also RhiSin-2.11 (TcMar-Tc2) 89.6 % over
146 bp on the opposite strand, and Mariner1_rFer / nMar4_rFer (Bat1K, *R. ferrumequinum*, TcMar-Tc2) 88 % over 77 bp.
A 237 bp TcMar-Tc2 family is a non-autonomous (MITE-like) element; its terminal inverted repeats fit the inverted copy of
X 91-115 found just upstream of X in 6 copies. The (TATTA)n tract at the start of X lies inside the family (positions
~82-90 on), not at an element's end. The 700 bp upstream of one member (NC_142518.1:59695071-59695800) is an L1 ORF2
fragment (L1M4_orf2 82 % over 410 bp): local context only.

## Why the candidates could not be decided: the partners are known repeats, and the other side was untested (2026-10-17)
All 60 partners (>= 30 bp) searched against Dfam curated (ssearch36 -z 11, E <= 1e-3): 56 are named, mostly bat L1s and
mariner / Tc2 MITEs at 92-100 % (`partners_vs_dfam.best.tsv`). For every group, the members' OTHER side (300 bp after the
element) was searched against the partner's own Dfam family (`partners_dfam_otherside.tsv`: hits >= 70 % over >= 50 bp, and
the median start of the hit on the family). In most groups the other side continues the same repeat where the partner
stops (within ~20 bp, the target-site duplication): the copies sit INSIDE an older repeat, i.e. an insertion hot spot,
shared on both sides. Stage 9 did not see the other side because those flanks were untestable (copy_status "masked": the
repeat filter blanks bases of 20-mers frequent in the genome), and fs9b read "no other-side pair" as "other side unique".
Family-specific candidates: r1_r3_P1_5p34 other side inside RhiSin-2.11 (MITE) in 38/41 (masked 38/41); 5p15 15/16 (L1,
masked 12/16); 5p13 9/10; 5p28 20/24; 5p25 11/13; r1_r3_P2_5p10 11/14. MEG-RS_5p1: array family (array_flag.tsv: 85.9 % of
MEG-RS in arrays), members regularly spaced (3.6, 6.2, 10.2 kb; 4.8 + 1.4 kb), 39/41 other sides masked; stage 9's array
class needs neighbours < 2 kb apart. Only the X group (r5_r6_P26_5p1) has a tested-unique other side (21/21 "unique",
other sides unrelated): its partner ends at RhiSin-1.14's own 3' end (237 of 237).
