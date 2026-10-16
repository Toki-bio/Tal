# Review alignments for the composite decisions (2026-10-17)

**Composites table** (`<composite>.table.aln.fa`): the composite consensus, then each declared unit part laid at its
placed interval (MAFFT `--addfragments --keeplength` against that interval only, so a unit cannot drift onto a shared
head). Row names give the declared part, the placed interval and its identity. `r10_r8_P18` also shows r8 where it
aligns (2-133, on r10's head): nothing of r8 aligns to P18's second half.

**r6 vs P26** (`r6_vs_P26_review.aln.fa`): copies called r6 or P26 by main (296fa33) or the branch (coverage-tiebreak
6e8eacf) when the two disagree, plus 15 + 15 copies both call firm r6 / firm P26 (references). Each copy cut from the
genome with 150 bp flanks (lower case, degapped, packed against the copy's own element ends); the element is the part
that aligns to r6 or P26, laid on P26's own positions (MAFFT `--addfragments --keeplength`: 376 element columns =
P26 1-376; r5 rows = r5 1-127, the part inside P26). Bases a copy has in addition to P26 are left out of the columns and
counted in `ins=`. `other=NNbp@PP%` = the copy's sequence upstream of its r6 match aligned against P26's r5 +
linker part (P26 1-137) alone, so locus overhang and internal insertions do not count. Rows sorted by `other`.
