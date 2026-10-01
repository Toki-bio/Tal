# rsi_v7

Full SINEderella report for *Rhinolophus sinicus* with the **bank-cleaning bug fixed** (SINEderella commit 2452058). In rsi_v3-v6 the cleaning merged r9 into r7, r4 into r2 and MEG-RS into MEG-RL by position-by-position identity, and the merged name took the other's sequence. Here the 21-sequence bank (`bank.pre_canon.fa`, `consensuses.clean.fa`) is unmerged, so r7 plates start from the true 154 bp, MEG-RL from 207 bp, and r9, r4, r2 and MEG-RS have their own plates.

- `report.html`, `alignments/` (60 files: top100, rand100, subfamily plates), `assignment_stats.tsv`
- `length_variants/`: SINEderella `tools/length_variants.py` results (end modes, linkage, TSD, residual) for the candidate pairs and single-family controls; `summary.tsv` is the pair table. Method and calibration: SINEderella `docs/LENGTH_VARIANTS.md`. The positive control (rle MEG-RS / MEG-RL) is in `chiroptera/length_variants/`.

Caveats: the MEG-RS top100 plate row 1 (proposed extension) is 1,086 bp because the border loop walked the 5' side to its cap (as in v6 on other genomes); not read plate by plate.
