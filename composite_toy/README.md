# Composite / dimer toy (SINEderella coverage tie-break), 2026-10-16

Two families, each with unit A (120 bp), unit B (150 bp) and dimer AB = A + 25 bp linker + B:
X = dimer halves identical to the units; Y = dimer halves 4 % different from the units.
40 planted copies per class (A, B, AB, ABt = 5'-truncated dimer) per family, 5-15 % diverged, 2 Mb genome
(therioserver `/home/toki/tmp/composite_toy/small`, generator `make_composite_toy.py`, seed 11).

Each file: the family's three consensus rows (`CONS_*`) + up to 40 planted copies of one true class (planted element,
oriented). Row names: `<true class>[_eNNN]_<locus>|main=<call>|fix=<call>`. `eNNN` = bp of A/linker kept before B
(truncated dimers only). `~` = not firm (no_unanimous / rejected). main = SINEderella 296fa33,
fix = branch `coverage-tiebreak` 1281386. Runs: `/home/toki/tmp/composite_toy/s3_main`, `s3_branch`.
