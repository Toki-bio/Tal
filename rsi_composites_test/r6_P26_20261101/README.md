# rsi: r6_210seqs vs r5_r6_P26 - who wins the vote, and should it (2026-11-01)

Replaces the 2026-10-16/17 review files one folder up (those compared main with the unmerged `coverage-tiebreak`
branch, on an older run). This one is built from the current run `rsi_r11_20261101` (step 2 of SINEderella main;
r11 does not touch r6 or P26), with copies of all four kinds labelled in the row names.

**P26** = r5-type head (positions 2-137) + the whole of r6 (138-376, 100 % identical to r6_210seqs). A copy that has
only r6 sequence therefore aligns to r6 and to P26 with the same alignment, and the 10 voting cycles split by score
noise (`-z 11`). The question: which copies are standalone r6, which are P26?

## Plate: `r6_P26_review.aln.fa` ([ViewAlign](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi_composites_test%2Fr6_P26_20261101%2Fr6_P26_review.aln.fa&title=rsi%20r6%20vs%20P26%3A%20A%20firm%20r6%2C%20B%20split%20votes%2C%20C%20P26%20with%20r6%20runner-up%2C%20D%20firm%20P26))

Rows, top to bottom:

| rows | group | what it is |
|---|---|---|
| 1-3 | consensus | `CONS_r5_r6_P26`, then r6 and r5_27seqs laid on P26's columns (r6 starts at P26 138) |
| A (40) | firm r6 | assigned r6, 10/10 votes |
| B (40) | split r6/P26 | votes for both in the 10 cycles; **not assigned** (no unanimous vote) |
| C (40) | P26, runner-up r6 | assigned P26 10/10, r6 second |
| D (40) | firm P26 | assigned P26 10/10, runner-up another composite (reference for real P26 copies) |

Only copies whose flanks fs9 calls unique are shown (fixed 5' mask; flanks shared or masked are left out), 40 random per
group, sorted within the group by the first P26 position they reach.

**How to read a row.** The middle block of 376 columns is P26: each copy's bases sit at the P26 position they align to,
in **upper case** (insertions in the copy are left out and counted in `ins=`). Everything of the copy's genome window
(locus + 300 bp upstream + 100 bp downstream) that does not align is in **lower case**, packed against the aligned part:
the upstream sequence ends right before the copy's first aligned P26 column. So:
- a **standalone r6** copy: lower-case (its own flank) in the head columns, upper case from about column 138 on;
- a **P26** copy: upper case from column 1 (or from its truncation point), its flank before that.

**Look at the head columns (P26 1-137) and at column 138**, group by group.

Row names: `group | votes r6:N_P26:N | locus | P26span=first-last P26 position aligned | head=bp@identity%/bits_null<bits> | ins=`.
`head` = the copy's sequence upstream of its r6 part aligned to P26 1-137 alone; `null` = the best score of the same
sequence against 20 dinucleotide shuffles of P26 1-137. A head counts when bits > null.

## Measured (`summary.tsv`, `r6_P26_review.tsv`; all copies of the run for the counts)

| group | copies | flanks unique | first P26 position: <= 20 / 21-119 / >= 120 (of 40 shown) | real r5 head (bits > null, >= 50 bp) |
|---|---|---|---|---|
| A firm r6 | 6,243 | 4,973 | 0 / 1 / 39 | 0 of 40 |
| B split r6/P26 | 1,188 | 958 | 0 / 6 / 34 | 1 of 40 |
| C P26, runner r6 | 497 | 358 | 25 / 11 / 4 | 31 of 40 |
| D firm P26 | 5,165 | 3,991 | 40 / 0 / 0 | 38 of 40 |

- B looks exactly like A: alignments start at r6's start (median P26 137); the six rows that start a little earlier
  creep a few bases into the head by chance (heads of 11-43 bp, all at or below the shuffle score).
- C is mostly real P26: 25 full heads, 10 5'-truncated P26 with 58-137 bp of real head (bits 64-127 vs null 27-41);
  but 5 of 40 have no head (4 start at 127-156, one with a head below the shuffles): r6 copies P26 won 10 times out of 10.
- D: real P26 throughout (the two "short head" rows have long diverged heads aligned in the main alignment).

## Analysis and proposed verdict

**B (1,188 copies, 958 with unique flanks) are standalone r6.** They carry no part of P26's r5 head; P26 gets some of
their votes only because it contains r6 entire and its score on the same alignment differs by `-z 11` noise. Today
they are unassigned, so r6 loses them: it is the bug you describe (the standalone loses the vote to the composite that
contains it). About 12 % of C (about 60 copies) are the same case where P26 happened to win all 10 cycles.

**Proposed fix, in the workflow (step 2):** when a copy's best and runner-up are a composite and a unit it contains, the
composite wins only if the copy's alignment covers at least half of the composite's sequence outside that unit; otherwise
the unit wins. The unit's interval inside the composite is taken from `composite_units.tsv`, which the run already
computes from the bank (P26: r6 = 138-376), so nothing is declared by hand. This is the rule of the parked
`coverage-tiebreak` branch with its hand-made table replaced by the automatic one. Expected here: B -> r6 (~1,000
copies), the headless part of C -> r6, real P26 (C with head, D) unchanged. If you agree, I will implement it on current
main, toy-test it (unit / composite / truncated-composite copies) and rerun rsi.

Script: `r6p26_review.py` (therioserver `/home/toki/tmp/r6p26/`; fs9 on all 13,575 r6/P26 copies: `fs9/copy_status.tsv`).
