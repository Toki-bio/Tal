# rsi_fresh2

Second full run from the same 21-consensus bank, after the mafft `--threadit 0` fix (SINEderella c34d052): the alignments are reproducible. Same files as `rsi_fresh/`, plus `stage8_r5/` (the r5 single-copy / control / in-composite plates from flankscan stage 8, deterministic, same on 16 and 32 threads).

Against rsi_fresh (run 1): the final consensus bank is identical (21 sequences), length-version verdicts identical, assignment counts differ slightly because ssearch36 `-z 11` shuffles (by design): most families within 1 %, r6 7,068 -> 6,585, P48 512 -> 470, MEG-RL 15 -> 12. Alignments differ in 52 of 56 files (run 1 was not reproducible).
