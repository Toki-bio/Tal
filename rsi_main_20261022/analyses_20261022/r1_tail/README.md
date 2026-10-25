# What follows the 3' end of r1_9seqs copies (rsi_main_20261022, 2026-10-26)

r1's consensus (155 bp) ends at `...CTGCAACTAACAA`, where r4 and r3 go on. Checked on the 21 copies assigned to r1
(`r1_tail.py`, therioserver `/home/toki/tmp/r1tail/`; copies padded 60 bp 5', 1 kb 3'); 7 have a complete r1 match.

| copy | verdict | after r1's end |
|---|---|---|
| NC@U@142514.1:146321124-146321306(+) | SINE | A-rich tail at once (`AAAAAAAACAAAACCAAAA`) |
| NC@U@142516.1:153233188-153233343(-) | SINE | A-rich tail at once; 10-bp TSD 13 bp past r1's end |
| NC@U@142514.1:17801016-17801172(+) | SINE | unrelated sequence |
| NC@U@142505.1:27424682-27424832(-) | SINE | unrelated sequence |
| NC@U@142504.1:13868426-13868586(-) | grey | r2's own tail (`TTGAACACGGCACCTTGAG...`, r2 106-138), then A-rich |
| NC@U@142518.1:103814724-103814878(+) | grey | the r3-type terminal piece (`TCCTGTTCCCC...AATAAAA`, P2 252-304) |
| NC@U@142515.1:55908440-55908594(-) | grey | r2's tail and then a second module: an r1_r3_P1 copy (P1 162-319, 93 %) |

TSD within 300 bp past r1's end: 1 of 7 (shuffled pairs 0 of 7).

Copies with r1's head are mostly assigned to the composites that start with it: r1_r3_P1 14,471, r1_r3_P2 6,910,
r1_r2_r3_r3_C11 2,907, r1 alone 21. In those consensuses r1 goes on with r2's tail (`TTGAACACGGCACC...`, P1 157-200)
and then a second element. The 21 stand-alone r1 copies are the few that stop at r1's end or soon after (2 of 7 with an
A-rich tail right there).

`rsi_r1_copies_3prime.aln.fa`: r1, r4, r3 and r1_r3_P1 consensuses + the 7 complete copies, each from r1's start to 250 bp
past its end (to the TSD for the one copy with one), MAFFT --localpair. `r1_22.report.txt`: the script's table.
