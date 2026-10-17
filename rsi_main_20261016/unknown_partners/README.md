# rsi: unknown partners from one-sided shared flanks (flankscan stage 9b, draft), 2026-10-17

SINEderella branch `fs9b-partners` (3c51893), `flankscan/fs9b_partners.sh`, run on this rsi run with the bat partner
library (therioserver `~/refs/partnerlib/partnerlib_9397.fa`). Stage 9 links copies whose flanks are similar from the
junction outward; its body check marks pairs with independent bodies "flank_only" (hot-spot signature) and drops them.
9b keeps flank_only groups shared on one side only (>= 10 copies, >= 3 contigs, other side shared by <= 25 %), builds the
partner, names it if the bank or the partner library knows it, counts copies of every family carrying it at the junction.

`unknown_partners.summary.tsv`: one row per distinct partner. kind = known partner (bank / library) | family-specific:
composite candidate (>= 50 % of the carriers in one family and >= 2x that family's share of all copies) | all families:
background. `unknown_partners.tsv`: one row per group.

Review alignments of the family-specific candidates (`<id>.review.aln.fa`): first row = the candidate composite (partner +
family consensus); members laid on its positions (MAFFT --addfragments --keeplength), outer flanks (100 bp) lower case,
degapped, packed against each member's own first / last aligned base; `ins=` = bases a member has beyond the composite.
`r5_r6_P26_5p1` is the X + r6 group found by eye (14 of its 30 copies are the hand-found ones).
