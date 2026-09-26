# Sicista — SINEderella database arm

## 2026-09-26 — step 1, Mammalia-only bank

- Genome: therioserver `~/Sicista2026/assembly/consensus.fasta` (Medaka assembly).
- Bank: the 58 SINEbase consensuses classed Mammalia (`alignments/sicista_consensuses.fa`, therioserver
  `~/sine_runs/Sicista/db_based/mammal_consensus.fa`). The first attempt with the full combined bank was too slow and was
  stopped (partial kept as `partial_fullbank_run_20260926_143745`).
- Run: `~/sine_runs/Sicista/db_based/run_20260926_150828`, THREADS=32, CHUNK_BP=30000, FLANK=50; stopped on purpose
  when step 2 started (15:08 -> 16:25).
- Result: 2,170,177 merged hits. `alignments/sicista_subfam_input_30k.aln.fa` = SubFam `input.clw.al`:
  600 chunk consensi of a 30,000-copy sample + the 58 bank consensuses (658 rows).
- Per-consensus raw hit counts are on sicista.html (overlapping; not copy numbers).
