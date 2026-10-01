# rsi_v6 — Rhinolophus sinicus, flank-scan elements (2026-10-01/02)

Pages
- [report.html](report.html) — the SINEderella report with the seven accepted flank-scan elements in the bank.
- [decisions.html](decisions.html) — the open decisions, each with its data.
- rsi composition page (Claude artifact): https://claude.ai/artifact/PfHq2MnKWxT4NcA76qJRkM

Files
- `alignments/` — plates of every unit and element: top 100, random 100, subfamily consensi (open in [MSA viewer](https://toki-bio.github.io/MSA-viewer/)).
- `stage8_plates/<family>_<S|B|A>.plate.aln.fa` — stage 8: singles (S), the same family inside composites (B), control singles (A); consensus row 1, 100 bp raw flank each side, element in upper case.
- `stage8_tsd/<family>/` — `summary.tsv`, `tsd_calibration.tsv` (real vs shuffled by minimum length), `ends.tsv`, `copies.tsv`.
- `consensus_pairs/` — pairwise alignments (original bank sequences) of the consensuses that the bank cleaning had merged by mistake (r9/r7, r4/r2, MEG-RS/MEG-RL). **In this v6 report the row 1 of the r7 plates is the r9 sequence and that of the MEG-RL plates the MEG-RS sequence** (bug in `canonicalize_consensus_bank.py`, fixed 2026-10-02); the rsi_v7 report replaces it.
- `summary.by_subfam.tsv`, `assignment_stats.tsv` — report tables.

Method notes: SINEderella `flankscan/HANDOFF.md` (stage 8: TSD search with a 45 bp 3′ allowance, flanks cut from the raw copy).
