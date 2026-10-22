# r10_19seqs, all copies

- **`rsi_r10_19seqs_all_plate.aln.fa` (use this one).** All 593 SINE copies plus the 168 grey copies (marked `[grey]`) of
  r10_19seqs, built by SINEderella's own plate chain with `PLATE_N=1000` (step8a sorted by score, then
  boundary_justify, rebuild_consensus_row, trim_display_flanks, correct_published_aln, add_seed_row; SINEderella main
  4305f5c; MAFFT L-INS-i on 32 threads). Copies are cut at the element edges and the flanks are packed against them, so the
  conserved 3' tail (G-rich, then A-rich) lines up in the columns right after the element. Row 1 = consensus rebuilt from
  the copies, row 2 = the consensus as searched. The 18 copies with verdict not SINE (shared flanks) are not shown.
  Server work dir: therioserver `/home/toki/tmp/r10_allplate32`.
- `rsi_r10_19seqs_all_SINE_grey.aln.fa`, `..._packed_after853.aln.fa`: superseded. Made by a separate script that aligned the
  raw loci (84-519 bp) whole and appended the flanks, which smeared the 3' tail; kept only because they were linked.
