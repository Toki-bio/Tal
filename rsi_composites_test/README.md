# rsi: unit vs declared composite, main vs branch coverage-tiebreak (2026-10-16)

Copies whose step-2 call differs between SINEderella main (296fa33) and branch `coverage-tiebreak` (d7c09b9) with the
declared composites table (`therioserver:/home/toki/rhin/rsi_tiebreak_20261016/rsi_composites_explicit.tsv`), on the same
step-1 copies (`/home/toki/rhin/rsi_main_20261016/run_20261016_194647`). Up to 60 copies per file, evenly sampled.
Rows: consensus rows `CONS_*`, then copies `<locus>|main=<call>|fix=<call>|reach=<bp>`; `~` = not firm;
reach = bp of composite sequence the copy's alignment covers upstream (5') of the unit's interval in the composite
(for r10 vs P18 the extra lies 3', so reach is 0 there). Rows are sorted by reach.
