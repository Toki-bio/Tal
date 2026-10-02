# rsi_fresh

Full SINEderella run on *Rhinolophus sinicus* started **from the 21 established consensuses** (the v7 bank: r1-r10, MEG-RL/RS/T2/TR, 7 composites), as one from-scratch run (not add-mode on top of an earlier run), for comparison with the semi-agnostic chain (rsi_peel10, flank scan, rsi_v7 add-mode).

- `report.html` (alignment links, Similarity section, length-version table), `alignments/`, `assignment_stats.tsv`, `consensuses.clean.fa`, `length_variants/`.

Comparison with rsi_v7 (assigned copies, fresh / v7): r10 779/851, P18 4239/4240, r5 1907/2090, P26 5659/5641, r6 7068/7062, r7 6847/6859, r8 3347/3395, r9 6850/6851, P1 14486/14389, P2 6877/6726, C11 2919/2924. Differences: r1 28/372, r3 50/284, MEG-RS 1695/818 (these are copies that sit between their bank neighbours; the composites take them).
Length-version verdicts are the same as v7 (r9/r7, r9/r8, r7/r8 TWO_VERSIONS; r7/r5 UNLINKED_ENDS; r5_r5_P48/r5_r6_P26 SINGLE_MODE; MEG-RS/MEG-RL and r4/r3 NOT_TESTED).
asSINEment check against the flank-scan classes: single copies 84.7 % same unit, 10.4 % (3,064) of them take P18 (r10 + r8 piece): the same weak point as in v7.
